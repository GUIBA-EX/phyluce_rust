//! Small bounded worker pool for deterministic, ordered parallel work.
//!
//! Backed by `rayon`'s work-stealing scheduler rather than a hand-rolled
//! mutex+`VecDeque`+channel pool (which this module used to be): benchmarked
//! ~1.8-4.4x faster across repeated runs on a synthetic many-small-CPU-tasks
//! workload shaped like what these CLI commands actually do (see
//! `tests::bench_parallel_dispatch_hand_rolled_vs_rayon`), likely because
//! rayon's deques avoid a mutex lock/unlock per item dequeue. The public
//! `try_map_ordered` signature is unchanged, so none of its 14 call sites
//! needed to change.

use std::sync::atomic::{AtomicBool, Ordering};

pub(crate) fn catch_operation<R>(
    operation: impl FnOnce() -> anyhow::Result<R>,
) -> anyhow::Result<R> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(operation)).unwrap_or_else(|payload| {
        let message = payload
            .downcast_ref::<&str>()
            .copied()
            .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
            .unwrap_or("unknown panic payload");
        Err(anyhow::anyhow!("parallel worker panicked: {message}"))
    })
}

/// One item's outcome. Kept distinct from `Result` so that `Cancelled`
/// items (skipped after an earlier item failed) never compete with a real
/// `Failed` item for "earliest-index error" selection below -- otherwise a
/// skip that happens to land at a lower index than the failure that caused
/// it could shadow the real error message.
enum Outcome<R> {
    Done(R),
    Failed(anyhow::Error),
    Cancelled,
}

/// Apply `operation` with at most `workers` threads while returning results
/// in the same order as `items`. On the first error, already-dispatched
/// work finishes but not-yet-started items are skipped (best effort, same
/// as before); the earliest-input-index error is returned.
pub fn try_map_ordered<T, R, F>(
    items: Vec<T>,
    workers: usize,
    operation: F,
) -> anyhow::Result<Vec<R>>
where
    T: Send,
    R: Send,
    F: Fn(T) -> anyhow::Result<R> + Sync,
{
    anyhow::ensure!(workers > 0, "worker count must be greater than zero");
    if items.len() <= 1 || workers == 1 {
        return items
            .into_iter()
            .map(|item| catch_operation(|| operation(item)))
            .collect();
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .map_err(|e| anyhow::anyhow!("failed to build worker pool: {e}"))?;

    let cancelled = AtomicBool::new(false);
    let outcomes: Vec<Outcome<R>> = pool.install(|| {
        use rayon::prelude::*;
        items
            .into_par_iter()
            .map(|item| {
                if cancelled.load(Ordering::Acquire) {
                    return Outcome::Cancelled;
                }
                match catch_operation(|| operation(item)) {
                    Ok(value) => Outcome::Done(value),
                    Err(error) => {
                        cancelled.store(true, Ordering::Release);
                        Outcome::Failed(error)
                    }
                }
            })
            .collect()
    });

    let mut ordered: Vec<Option<R>> = Vec::with_capacity(outcomes.len());
    let mut first_error: Option<(usize, anyhow::Error)> = None;
    for (index, outcome) in outcomes.into_iter().enumerate() {
        match outcome {
            Outcome::Done(value) => ordered.push(Some(value)),
            Outcome::Cancelled => ordered.push(None),
            Outcome::Failed(error) => {
                ordered.push(None);
                if first_error
                    .as_ref()
                    .is_none_or(|(error_index, _)| index < *error_index)
                {
                    first_error = Some((index, error));
                }
            }
        }
    }
    if let Some((_, error)) = first_error {
        return Err(error);
    }

    ordered
        .into_iter()
        .map(|result| {
            result.ok_or_else(|| anyhow::anyhow!("parallel worker did not return a result"))
        })
        .collect()
}

pub fn ensure_unique_output_names(names: impl IntoIterator<Item = String>) -> anyhow::Result<()> {
    let mut seen = std::collections::HashSet::new();
    for name in names {
        anyhow::ensure!(
            seen.insert(name.clone()),
            "multiple inputs map to the same output file {name:?}"
        );
    }
    Ok(())
}

/// Long external jobs benefit from one shared priority queue rather than
/// indexed range splitting: workers take the largest remaining job first.
/// Results and error selection retain original input order. The queue lock
/// is held only while popping a job, never during the expensive operation.
pub fn try_map_ordered_by_cost<T, R, F>(
    items: Vec<T>,
    workers: usize,
    cost: impl Fn(&T) -> u128,
    operation: F,
) -> anyhow::Result<Vec<R>>
where
    T: Send,
    R: Send,
    F: Fn(T) -> anyhow::Result<R> + Sync,
{
    anyhow::ensure!(workers > 0, "worker count must be greater than zero");
    if workers == 1 || items.len() <= 1 {
        return try_map_ordered(items, workers, operation);
    }
    let count = items.len();
    let workers = workers.min(count);
    let mut ranked: Vec<_> = items
        .into_iter()
        .enumerate()
        .map(|(index, item)| (cost(&item), index, item))
        .collect();
    // pop() must select highest cost, then smallest original index.
    ranked.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| b.1.cmp(&a.1)));
    let queue = std::sync::Mutex::new(ranked);
    let cancelled = AtomicBool::new(false);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()?;
    let batches: Vec<Vec<_>> = pool.install(|| {
        use rayon::prelude::*;
        (0..workers)
            .into_par_iter()
            .map(|_| {
                let mut results = Vec::new();
                loop {
                    if cancelled.load(Ordering::Acquire) {
                        break;
                    }
                    let job = queue.lock().unwrap().pop();
                    let Some((_, index, item)) = job else {
                        break;
                    };
                    let result = catch_operation(|| operation(item));
                    if result.is_err() {
                        cancelled.store(true, Ordering::Release);
                    }
                    results.push((index, result));
                }
                results
            })
            .collect()
    });
    let mut ordered: Vec<_> = batches.into_iter().flatten().collect();
    ordered.sort_by_key(|(index, _)| *index);
    let results: Vec<R> = ordered
        .into_iter()
        .map(|(_, result)| result)
        .collect::<anyhow::Result<_>>()?;
    anyhow::ensure!(
        results.len() == count,
        "priority worker did not return a result"
    );
    Ok(results)
}

/// Process bounded batches on one pool, consuming results in input order on
/// the calling thread (e.g. the SQLite connection owner). At most 2 * workers
/// results are retained; a slow database cannot accumulate the entire run.
pub fn try_for_each_ordered<T, R, F, C>(
    items: &[T],
    workers: usize,
    operation: F,
    mut consume: C,
) -> anyhow::Result<()>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> anyhow::Result<R> + Sync,
    C: FnMut(R) -> anyhow::Result<()>,
{
    anyhow::ensure!(workers > 0, "worker count must be greater than zero");
    if workers == 1 || items.len() <= 1 {
        for item in items {
            consume(catch_operation(|| operation(item))?)?;
        }
        return Ok(());
    }
    let workers = workers.min(items.len());
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()?;
    for batch in items.chunks(workers.saturating_mul(2)) {
        let cancelled = AtomicBool::new(false);
        let results: Vec<Outcome<R>> = pool.install(|| {
            use rayon::prelude::*;
            batch
                .par_iter()
                .map(|item| {
                    if cancelled.load(Ordering::Acquire) {
                        return Outcome::Cancelled;
                    }
                    match catch_operation(|| operation(item)) {
                        Ok(result) => Outcome::Done(result),
                        Err(error) => {
                            cancelled.store(true, Ordering::Release);
                            Outcome::Failed(error)
                        }
                    }
                })
                .collect()
        });
        for result in results {
            match result {
                Outcome::Done(result) => consume(result)?,
                Outcome::Failed(error) => return Err(error),
                Outcome::Cancelled => (), // A real error in this batch follows.
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::time::Duration;

    // Ad hoc benchmark: is rayon's work-stealing scheduler actually faster
    // than this module's hand-rolled mutex+VecDeque+channel pool for the
    // shape of work these CLI commands do (many short-to-medium CPU tasks,
    // bounded worker count), or is swapping it in mostly a code-size win?
    // `rayon` is a dev-dependency only here -- not pulled into the
    // production binary unless this comparison justifies it. Run with:
    //   cargo +stable test --release -p phyluce-cli --bin phyluce -- --ignored --nocapture bench_parallel_dispatch
    fn cpu_work(n: u64) -> u64 {
        // A deterministic, allocation-free stand-in for "parse+process one
        // alignment file": no sleeps, no syscalls, just enough scalar work
        // to be measurable and non-trivial for the optimizer to elide.
        let mut acc = n;
        for _ in 0..2000 {
            acc = acc.wrapping_mul(6364136223846793005).wrapping_add(1);
        }
        acc
    }

    fn rayon_map_ordered(items: Vec<u64>, workers: usize) -> Vec<u64> {
        use rayon::prelude::*;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| items.into_par_iter().map(cpu_work).collect())
    }

    #[test]
    #[ignore]
    fn bench_parallel_dispatch_hand_rolled_vs_rayon() {
        for (n_items, workers) in [(2_000usize, 4usize), (2_000, 8)] {
            let items: Vec<u64> = (0..n_items as u64).collect();

            let start = std::time::Instant::now();
            let hand_rolled = try_map_ordered(items.clone(), workers, |n| Ok(cpu_work(n))).unwrap();
            let hand_rolled_elapsed = start.elapsed();

            let start = std::time::Instant::now();
            let rayon_result = rayon_map_ordered(items, workers);
            let rayon_elapsed = start.elapsed();

            assert_eq!(hand_rolled, rayon_result);
            eprintln!(
                "[bench] {n_items} items, {workers} workers: hand-rolled {:?} vs rayon {:?} ({:.2}x)",
                hand_rolled_elapsed,
                rayon_elapsed,
                hand_rolled_elapsed.as_secs_f64() / rayon_elapsed.as_secs_f64()
            );
        }
    }

    #[test]
    fn preserves_order_and_respects_worker_limit() {
        let active = AtomicUsize::new(0);
        let peak = AtomicUsize::new(0);
        let output = try_map_ordered((0..12).collect(), 3, |value| {
            let now = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
            std::thread::sleep(Duration::from_millis(5));
            active.fetch_sub(1, Ordering::SeqCst);
            Ok(value * 2)
        })
        .unwrap();

        assert_eq!(output, (0..12).map(|value| value * 2).collect::<Vec<_>>());
        assert!(peak.load(Ordering::SeqCst) > 1);
        assert!(peak.load(Ordering::SeqCst) <= 3);
    }

    #[test]
    fn propagates_worker_errors_in_input_order() {
        let error = try_map_ordered((0..4).collect(), 2, |value| {
            anyhow::ensure!(value != 2, "failed at {value}");
            Ok(value)
        })
        .unwrap_err();
        assert_eq!(error.to_string(), "failed at 2");
    }

    #[test]
    fn cancels_queued_work_after_an_error() {
        let completed = AtomicUsize::new(0);
        let error = try_map_ordered((0..100).collect(), 2, |value| {
            if value == 0 {
                anyhow::bail!("stop");
            }
            std::thread::sleep(Duration::from_millis(5));
            completed.fetch_add(1, Ordering::SeqCst);
            Ok(value)
        })
        .unwrap_err();
        assert_eq!(error.to_string(), "stop");
        assert!(completed.load(Ordering::SeqCst) < 100);
    }

    #[test]
    fn converts_worker_panics_to_errors() {
        let error = try_map_ordered(vec![1], 1, |_| -> anyhow::Result<()> {
            panic!("broken worker")
        })
        .unwrap_err();
        assert!(error.to_string().contains("broken worker"));
    }

    #[test]
    fn rejects_output_name_collisions_before_parallel_work() {
        let error =
            ensure_unique_output_names(["a.nexus", "a.nexus"].map(str::to_string)).unwrap_err();
        assert!(error.to_string().contains("a.nexus"));
    }
}

#[cfg(test)]
mod bounded_tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;

    #[test]
    fn ordered_consumer_bounds_live_results_and_threads() {
        let active = AtomicUsize::new(0);
        let peak_active = AtomicUsize::new(0);
        let live = AtomicUsize::new(0);
        let peak_live = AtomicUsize::new(0);
        let mut consumed = Vec::new();
        try_for_each_ordered(
            &(0..40).collect::<Vec<_>>(),
            3,
            |&index| {
                let n = active.fetch_add(1, Ordering::SeqCst) + 1;
                peak_active.fetch_max(n, Ordering::SeqCst);
                std::thread::sleep(std::time::Duration::from_millis((index % 3) + 1));
                active.fetch_sub(1, Ordering::SeqCst);
                let n = live.fetch_add(1, Ordering::SeqCst) + 1;
                peak_live.fetch_max(n, Ordering::SeqCst);
                Ok(index)
            },
            |index| {
                consumed.push(index);
                live.fetch_sub(1, Ordering::SeqCst);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(consumed, (0..40).collect::<Vec<_>>());
        assert!(peak_active.load(Ordering::SeqCst) <= 3);
        assert!(peak_live.load(Ordering::SeqCst) <= 6);
    }

    #[test]
    fn errors_and_panics_stop_later_batches() {
        for panic in [false, true] {
            let started = AtomicUsize::new(0);
            let result = try_for_each_ordered(
                &(0..100).collect::<Vec<_>>(),
                2,
                |&index| {
                    started.fetch_add(1, Ordering::SeqCst);
                    if index == 0 {
                        assert!(!panic, "worker failure");
                        anyhow::bail!("worker failure");
                    }
                    Ok(index)
                },
                |_| Ok(()),
            );
            assert!(result.unwrap_err().to_string().contains("worker failure"));
            assert!(started.load(Ordering::SeqCst) <= 4);
        }
        assert!(try_for_each_ordered(&[1], 0, |&n| Ok(n), |_| Ok(())).is_err());
    }
}

#[cfg(test)]
mod priority_tests {
    use super::*;
    #[test]
    fn priority_preserves_results_and_bounds_workers() {
        let active = std::sync::atomic::AtomicUsize::new(0);
        let peak = std::sync::atomic::AtomicUsize::new(0);
        let result = try_map_ordered_by_cost(
            (0..100).collect(),
            4,
            |x| *x as u128,
            |x| {
                let n = active.fetch_add(1, Ordering::SeqCst) + 1;
                peak.fetch_max(n, Ordering::SeqCst);
                std::thread::yield_now();
                active.fetch_sub(1, Ordering::SeqCst);
                Ok(x * 2)
            },
        )
        .unwrap();
        assert_eq!(result, (0..100).map(|x| x * 2).collect::<Vec<_>>());
        assert!(peak.load(Ordering::SeqCst) <= 4);
    }
    #[test]
    fn priority_starts_heaviest_jobs_first() {
        let barrier = std::sync::Barrier::new(2);
        let started = std::sync::Mutex::new(Vec::new());
        try_map_ordered_by_cost(
            vec![1, 2, 3, 4],
            2,
            |x| *x as u128,
            |x| {
                started.lock().unwrap().push(x);
                if x >= 3 {
                    barrier.wait();
                }
                Ok(x)
            },
        )
        .unwrap();
        let mut first = started.lock().unwrap()[..2].to_vec();
        first.sort();
        assert_eq!(first, vec![3, 4]);
    }
    #[test]
    fn priority_propagates_failure_and_panic() {
        assert!(try_map_ordered_by_cost(vec![1, 2], 0, |_| 1, Ok).is_err());
        let err = try_map_ordered_by_cost(
            vec![1, 2],
            2,
            |_| 1,
            |_| -> anyhow::Result<()> { panic!("priority failure") },
        )
        .unwrap_err();
        assert!(err.to_string().contains("priority failure"));
        assert!(try_map_ordered_by_cost(
            vec![1, 2],
            2,
            |_| 1,
            |_| -> anyhow::Result<()> { anyhow::bail!("failed") }
        )
        .is_err());
    }
}
