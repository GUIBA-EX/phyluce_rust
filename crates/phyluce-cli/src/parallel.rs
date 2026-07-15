//! Small bounded worker pool for deterministic, ordered parallel work.

use std::collections::VecDeque;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

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

/// Apply `operation` with at most `workers` threads while returning results
/// in the same order as `items`.
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

    let item_count = items.len();
    let queue = Arc::new(Mutex::new(
        items.into_iter().enumerate().collect::<VecDeque<_>>(),
    ));
    let (sender, receiver) = std::sync::mpsc::channel();
    let cancelled = AtomicBool::new(false);

    std::thread::scope(|scope| {
        for _ in 0..workers.min(item_count) {
            let queue = Arc::clone(&queue);
            let sender = sender.clone();
            let operation = &operation;
            let cancelled = &cancelled;
            scope.spawn(move || loop {
                if cancelled.load(Ordering::Acquire) {
                    break;
                }
                let item = {
                    let mut queue = queue
                        .lock()
                        .unwrap_or_else(|poisoned| poisoned.into_inner());
                    if cancelled.load(Ordering::Acquire) {
                        return;
                    }
                    queue.pop_front()
                };
                let Some((index, item)) = item else {
                    break;
                };
                let result = catch_operation(|| operation(item));
                if result.is_err() {
                    cancelled.store(true, Ordering::Release);
                    queue
                        .lock()
                        .unwrap_or_else(|poisoned| poisoned.into_inner())
                        .clear();
                }
                if sender.send((index, result)).is_err() {
                    break;
                }
            });
        }
        drop(sender);

        let mut ordered: Vec<Option<R>> = (0..item_count).map(|_| None).collect();
        let mut first_error: Option<(usize, anyhow::Error)> = None;
        for (index, result) in receiver {
            match result {
                Ok(value) => ordered[index] = Some(value),
                Err(error) => {
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
    })
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::time::Duration;

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
