//! CLI wiring for `phyluce probe run-multiple-lastzs-sqlite`, mirroring
//! `phyluce_probe_run_multiple_lastzs_sqlite`.
//!
//! Like `phyluce.many_lastz.multi_lastz_runner`, chromosome genomes are
//! split by `.2bit` sequence and scaffold genomes are decoded into roughly
//! 10 Mbp FASTA chunks. LASTZ runs concurrently, then outputs are joined in
//! target order before one main-thread SQLite transaction per species.

use std::collections::{BTreeMap, HashMap};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use anyhow::Context;
use phyluce_io::sql::ident;
use phyluce_io::twobit::TwoBitFile;
use rusqlite::Connection;

const SCAFFOLD_CHUNK_SIZE: usize = 10_000_000;
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

pub struct RunMultipleLastzsArgs {
    pub chromolist: Vec<String>,
    pub scaffoldlist: Vec<String>,
    pub append: bool,
    pub no_dir: bool,
    pub genome_base_path: String,
    pub coverage: f64,
    pub identity: f64,
    pub cores: usize,
}

#[derive(Debug)]
struct TempPath {
    path: PathBuf,
}

impl Drop for TempPath {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

struct TargetUnit {
    target_spec: String,
    _input: Option<TempPath>,
}

#[cfg(test)]
struct GenomeUnits {
    genome: String,
    units: Vec<TargetUnit>,
}

struct UnitTask {
    genome_index: usize,
    unit_index: usize,
    unit: TargetUnit,
}

struct UnitResult {
    genome_index: usize,
    unit_index: usize,
    output: Option<TempPath>,
}

enum QueueEvent {
    GenomeStarted {
        genome_index: usize,
        genome: String,
    },
    GenomeFinished {
        genome_index: usize,
        unit_count: usize,
    },
    UnitFinished(UnitResult),
    Failed(anyhow::Error),
}

struct PendingGenome {
    genome: String,
    expected_units: Option<usize>,
    results: BTreeMap<usize, Option<TempPath>>,
}

fn create_temp_file(suffix: &str) -> anyhow::Result<(TempPath, File)> {
    for _ in 0..1000 {
        let serial = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "phyluce-lastz-{}-{serial}{suffix}",
            std::process::id()
        ));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => return Ok((TempPath { path }, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error.into()),
        }
    }
    anyhow::bail!("could not allocate a unique LASTZ temporary file")
}

fn genome_path(base_path: &str, no_dir: bool, name: &str) -> PathBuf {
    if no_dir {
        Path::new(base_path).join(format!("{name}.2bit"))
    } else {
        Path::new(base_path).join(name).join(format!("{name}.2bit"))
    }
}

fn create_species_lastz_table(conn: &Connection, g: &str) -> anyhow::Result<()> {
    let table = ident(g);
    let index = ident(&format!("{g}_name2_idx"));
    conn.execute_batch(&format!(
        "CREATE TABLE {table} (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            score INTEGER NOT NULL,
            name1 TEXT NOT NULL,
            strand1 TEXT NOT NULL,
            zstart1 INTEGER NOT NULL,
            end1 INTEGER NOT NULL,
            length1 INTEGER NOT NULL,
            name2 TEXT NOT NULL,
            strand2 TEXT NOT NULL,
            zstart2 INTEGER NOT NULL,
            end2 INTEGER NOT NULL,
            length2 INTEGER NOT NULL,
            diff TEXT NOT NULL,
            cigar TEXT NOT NULL,
            identity TEXT NOT NULL,
            percent_identity FLOAT NOT NULL,
            continuity TEXT NOT NULL,
            percent_continuity FLOAT NOT NULL,
            coverage TEXT NOT NULL,
            percent_coverage FLOAT NOT NULL);
         CREATE INDEX {index} on {table}(name2);"
    ))?;
    Ok(())
}

/// Mirrors `clean_lastz_data`: strip literal `%` characters (lastz writes
/// percentages like `100.0%` in several columns; the DB stores them as
/// bare numbers).
#[cfg(test)]
fn clean_lastz_data(text: &str) -> String {
    text.replace('%', "")
}

#[cfg(test)]
fn insert_species_rows(conn: &Connection, g: &str, cleaned: &str) -> anyhow::Result<()> {
    let table = ident(g);
    let mut stmt = conn.prepare(&format!(
        "INSERT INTO {table} (score, name1, strand1, zstart1, end1,
            length1, name2, strand2, zstart2, end2, length2, diff, cigar,
            identity, percent_identity, continuity, percent_continuity,
            coverage, percent_coverage) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    ))?;
    for line in cleaned.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        anyhow::ensure!(fields.len() == 19, "malformed lastz row: {line:?}");
        stmt.execute(rusqlite::params_from_iter(fields.iter()))?;
    }
    Ok(())
}

fn store_species_file(conn: &Connection, genome: &str, cleaned: &Path) -> anyhow::Result<()> {
    let tx = conn.unchecked_transaction()?;
    create_species_lastz_table(&tx, genome)?;
    {
        let table = ident(genome);
        let mut stmt = tx.prepare(&format!(
            "INSERT INTO {table} (score, name1, strand1, zstart1, end1,
                length1, name2, strand2, zstart2, end2, length2, diff, cigar,
                identity, percent_identity, continuity, percent_continuity,
                coverage, percent_coverage) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
        ))?;
        for line in BufReader::new(
            File::open(cleaned)
                .with_context(|| format!("opening cleaned lastz output {}", cleaned.display()))?,
        )
        .lines()
        {
            let line = line?;
            let fields: Vec<&str> = line.trim().split('\t').collect();
            anyhow::ensure!(fields.len() == 19, "malformed lastz row: {line:?}");
            stmt.execute(rusqlite::params_from_iter(fields.iter()))?;
        }
    }
    tx.execute("INSERT INTO species (name) VALUES (?1)", [genome])?;
    tx.commit()?;
    Ok(())
}

#[cfg(test)]
fn store_species_rows(conn: &Connection, genome: &str, cleaned: &str) -> anyhow::Result<()> {
    let tx = conn.unchecked_transaction()?;
    create_species_lastz_table(&tx, genome)?;
    insert_species_rows(&tx, genome, cleaned)?;
    tx.execute("INSERT INTO species (name) VALUES (?1)", [genome])?;
    tx.commit()?;
    Ok(())
}

#[cfg(test)]
fn chromosome_units(target: &Path) -> anyhow::Result<Vec<TargetUnit>> {
    let twobit = TwoBitFile::open(target)?;
    Ok(twobit
        .names()
        .into_iter()
        .map(|name| TargetUnit {
            target_spec: format!("{}/{name}", target.display()),
            _input: None,
        })
        .collect())
}

fn finish_scaffold_chunk(temporary: TempPath, mut file: File) -> anyhow::Result<TargetUnit> {
    file.flush()?;
    Ok(TargetUnit {
        target_spec: temporary.path.to_string_lossy().to_string(),
        _input: Some(temporary),
    })
}

#[cfg(test)]
fn scaffold_units(target: &Path, chunk_size: usize) -> anyhow::Result<Vec<TargetUnit>> {
    anyhow::ensure!(
        chunk_size > 0,
        "scaffold chunk size must be greater than zero"
    );
    let twobit = TwoBitFile::open(target)?;
    let names: Vec<String> = twobit.names().into_iter().map(str::to_string).collect();
    let mut chunks = Vec::new();
    let mut current: Option<(TempPath, File, usize)> = None;

    for name in names {
        if current.is_none() {
            let (temporary, file) = create_temp_file(".fasta")?;
            current = Some((temporary, file, 0));
        }
        let sequence = twobit.read_full(&name)?;
        let (_, file, length) = current.as_mut().expect("scaffold chunk was initialized");
        writeln!(file, ">{name}")?;
        file.write_all(&sequence)?;
        writeln!(file)?;
        *length += sequence.len();

        if *length > chunk_size {
            let (temporary, file, _) = current.take().expect("scaffold chunk exists");
            chunks.push(finish_scaffold_chunk(temporary, file)?);
        }
    }
    if let Some((temporary, file, _)) = current {
        chunks.push(finish_scaffold_chunk(temporary, file)?);
    }
    anyhow::ensure!(
        !chunks.is_empty(),
        "{} contains no sequences",
        target.display()
    );
    Ok(chunks)
}

#[cfg(test)]
fn execute_units<F>(
    units: Vec<TargetUnit>,
    cores: usize,
    execute: F,
) -> anyhow::Result<Vec<Option<TempPath>>>
where
    F: Fn(&str, &Path) -> anyhow::Result<()> + Sync,
{
    crate::parallel::try_map_ordered(units, cores, |unit| {
        let (temporary, file) = create_temp_file(".lastz")?;
        drop(file);
        execute(&unit.target_spec, &temporary.path)?;
        let has_results = temporary.path.metadata()?.len() > 0;
        Ok(has_results.then_some(temporary))
    })
}

#[cfg(test)]
fn execute_genome_queue<F>(
    genomes: Vec<GenomeUnits>,
    cores: usize,
    execute: F,
) -> anyhow::Result<Vec<(String, Vec<Option<TempPath>>)>>
where
    F: Fn(&str, &Path) -> anyhow::Result<()> + Sync,
{
    let mut boundaries = Vec::with_capacity(genomes.len());
    let mut units = Vec::new();
    for genome in genomes {
        boundaries.push((genome.genome, genome.units.len()));
        units.extend(genome.units);
    }

    let results = execute_units(units, cores, execute)?;
    let mut results = results.into_iter();
    let mut grouped = Vec::with_capacity(boundaries.len());
    for (genome, unit_count) in boundaries {
        grouped.push((genome, results.by_ref().take(unit_count).collect()));
    }
    debug_assert!(results.next().is_none());
    Ok(grouped)
}

#[cfg(test)]
fn concatenate_results(output: &Path, results: &[Option<TempPath>]) -> anyhow::Result<()> {
    let mut destination = File::create(output)?;
    for temporary in results.iter().flatten() {
        let mut source = File::open(&temporary.path)?;
        std::io::copy(&mut source, &mut destination)?;
    }
    destination.flush()?;
    Ok(())
}

#[cfg(test)]
fn align_genome(
    lastz_bin: &str,
    genome: &str,
    scaffolds: bool,
    probefile: &Path,
    output_dir: &Path,
    args: &RunMultipleLastzsArgs,
) -> anyhow::Result<String> {
    let target = genome_path(&args.genome_base_path, args.no_dir, genome);
    let units = if scaffolds {
        scaffold_units(&target, SCAFFOLD_CHUNK_SIZE)?
    } else {
        chromosome_units(&target)?
    };
    let query = probefile.to_string_lossy().to_string();
    let mut grouped = execute_genome_queue(
        vec![GenomeUnits {
            genome: genome.to_string(),
            units,
        }],
        args.cores,
        |target_spec, output| {
            crate::lastz_align::run_many_lastz(
                lastz_bin,
                target_spec,
                &query,
                args.coverage,
                args.identity,
                &output.to_string_lossy(),
            )
        },
    )?;
    let (_, results) = grouped.pop().expect("one genome was queued");

    finalize_genome(genome, probefile, output_dir, &results)
}

fn finalize_genome_file(
    genome: &str,
    probefile: &Path,
    output_dir: &Path,
    results: &[Option<TempPath>],
) -> anyhow::Result<PathBuf> {
    let probe_name = probefile
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("probes");
    let output =
        crate::output_path::output_file(output_dir, &format!("{probe_name}_v_{genome}.lastz"))?;
    let clean_path = PathBuf::from(format!("{}.clean", output.display()));
    let mut destination = File::create(&clean_path)
        .with_context(|| format!("creating cleaned lastz output {}", clean_path.display()))?;
    let mut input_buffer = [0u8; 64 * 1024];
    let mut clean_buffer = Vec::with_capacity(input_buffer.len());
    for temporary in results.iter().flatten() {
        let mut source = File::open(&temporary.path)
            .with_context(|| format!("opening lastz unit output {}", temporary.path.display()))?;
        loop {
            let count = source.read(&mut input_buffer)?;
            if count == 0 {
                break;
            }
            clean_buffer.clear();
            clean_buffer.extend(
                input_buffer[..count]
                    .iter()
                    .copied()
                    .filter(|byte| *byte != b'%'),
            );
            destination.write_all(&clean_buffer)?;
        }
    }
    destination.flush()?;
    Ok(clean_path)
}

#[cfg(test)]
fn finalize_genome(
    genome: &str,
    probefile: &Path,
    output_dir: &Path,
    results: &[Option<TempPath>],
) -> anyhow::Result<String> {
    Ok(std::fs::read_to_string(finalize_genome_file(
        genome, probefile, output_dir, results,
    )?)?)
}

fn send_task(
    sender: &std::sync::mpsc::SyncSender<UnitTask>,
    cancelled: &AtomicBool,
    mut task: UnitTask,
) -> anyhow::Result<()> {
    loop {
        anyhow::ensure!(
            !cancelled.load(Ordering::Acquire),
            "LASTZ task production was cancelled"
        );
        match sender.try_send(task) {
            Ok(()) => return Ok(()),
            Err(std::sync::mpsc::TrySendError::Full(returned)) => {
                task = returned;
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
            Err(std::sync::mpsc::TrySendError::Disconnected(_)) => {
                anyhow::bail!("LASTZ worker queue disconnected")
            }
        }
    }
}

fn produce_genome_tasks(
    sender: &std::sync::mpsc::SyncSender<UnitTask>,
    cancelled: &AtomicBool,
    genome_index: usize,
    genome: &str,
    scaffolds: bool,
    args: &RunMultipleLastzsArgs,
) -> anyhow::Result<usize> {
    let target = genome_path(&args.genome_base_path, args.no_dir, genome);
    let twobit = TwoBitFile::open(&target)
        .with_context(|| format!("opening 2bit genome file {}", target.display()))?;
    let names: Vec<String> = twobit.names().into_iter().map(str::to_string).collect();

    if !scaffolds {
        let unit_count = names.len();
        for (unit_index, name) in names.into_iter().enumerate() {
            send_task(
                sender,
                cancelled,
                UnitTask {
                    genome_index,
                    unit_index,
                    unit: TargetUnit {
                        target_spec: format!("{}/{name}", target.display()),
                        _input: None,
                    },
                },
            )?;
        }
        return Ok(unit_count);
    }

    let mut current: Option<(TempPath, File, usize)> = None;
    let mut unit_count = 0usize;
    for name in names {
        if current.is_none() {
            let (temporary, file) = create_temp_file(".fasta")?;
            current = Some((temporary, file, 0));
        }
        let sequence = twobit.read_full(&name)?;
        let (_, file, length) = current.as_mut().expect("scaffold chunk was initialized");
        writeln!(file, ">{name}")?;
        file.write_all(&sequence)?;
        writeln!(file)?;
        *length += sequence.len();
        if *length > SCAFFOLD_CHUNK_SIZE {
            let (temporary, file, _) = current.take().expect("scaffold chunk exists");
            send_task(
                sender,
                cancelled,
                UnitTask {
                    genome_index,
                    unit_index: unit_count,
                    unit: finish_scaffold_chunk(temporary, file)?,
                },
            )?;
            unit_count += 1;
        }
    }
    if let Some((temporary, file, _)) = current {
        send_task(
            sender,
            cancelled,
            UnitTask {
                genome_index,
                unit_index: unit_count,
                unit: finish_scaffold_chunk(temporary, file)?,
            },
        )?;
        unit_count += 1;
    }
    anyhow::ensure!(
        !scaffolds || unit_count > 0,
        "{} contains no sequences",
        target.display()
    );
    Ok(unit_count)
}

fn finalize_completed_genome(
    pending: &mut HashMap<usize, PendingGenome>,
    genome_index: usize,
    conn: &Connection,
    probefile: &Path,
    output_dir: &Path,
) -> anyhow::Result<bool> {
    let Some(state) = pending.get(&genome_index) else {
        return Ok(false);
    };
    let Some(expected_units) = state.expected_units else {
        return Ok(false);
    };
    if state.results.len() != expected_units {
        return Ok(false);
    }

    let mut state = pending
        .remove(&genome_index)
        .ok_or_else(|| anyhow::anyhow!("completed LASTZ genome disappeared"))?;
    let results = (0..expected_units)
        .map(|unit_index| {
            state.results.remove(&unit_index).ok_or_else(|| {
                anyhow::anyhow!(
                    "LASTZ results for genome {:?} are missing unit {unit_index}",
                    state.genome
                )
            })
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let cleaned = finalize_genome_file(&state.genome, probefile, output_dir, &results)?;
    store_species_file(conn, &state.genome, &cleaned)?;
    Ok(true)
}

fn align_and_store_all(
    conn: &Connection,
    lastz_bin: &str,
    probefile: &Path,
    output_dir: &Path,
    args: &RunMultipleLastzsArgs,
) -> anyhow::Result<()> {
    let query = probefile.to_string_lossy().to_string();
    let queue_capacity = args.cores.saturating_mul(2).max(1);
    let (task_sender, task_receiver) = std::sync::mpsc::sync_channel(queue_capacity);
    let task_receiver = Arc::new(Mutex::new(task_receiver));
    let (event_sender, event_receiver) = std::sync::mpsc::sync_channel(queue_capacity);
    let cancelled = Arc::new(AtomicBool::new(false));

    std::thread::scope(|scope| -> anyhow::Result<()> {
        {
            let event_sender = event_sender.clone();
            let cancelled = Arc::clone(&cancelled);
            scope.spawn(move || {
                let result = crate::parallel::catch_operation(|| {
                    let genomes = args
                        .scaffoldlist
                        .iter()
                        .map(|genome| (genome.as_str(), true))
                        .chain(
                            args.chromolist
                                .iter()
                                .map(|genome| (genome.as_str(), false)),
                        );
                    for (genome_index, (genome, scaffolds)) in genomes.enumerate() {
                        anyhow::ensure!(
                            !cancelled.load(Ordering::Acquire),
                            "LASTZ task production was cancelled"
                        );
                        event_sender.send(QueueEvent::GenomeStarted {
                            genome_index,
                            genome: genome.to_string(),
                        })?;
                        let unit_count = produce_genome_tasks(
                            &task_sender,
                            &cancelled,
                            genome_index,
                            genome,
                            scaffolds,
                            args,
                        )?;
                        event_sender.send(QueueEvent::GenomeFinished {
                            genome_index,
                            unit_count,
                        })?;
                    }
                    Ok(())
                });
                if let Err(error) = result {
                    if !cancelled.swap(true, Ordering::AcqRel) {
                        let _ = event_sender.send(QueueEvent::Failed(error));
                    }
                }
            });
        }

        for _ in 0..args.cores {
            let receiver = Arc::clone(&task_receiver);
            let event_sender = event_sender.clone();
            let cancelled = Arc::clone(&cancelled);
            let query = query.as_str();
            scope.spawn(move || {
                let result = crate::parallel::catch_operation(|| loop {
                    if cancelled.load(Ordering::Acquire) {
                        return Ok(());
                    }
                    let task = receiver
                        .lock()
                        .unwrap_or_else(|poisoned| poisoned.into_inner())
                        .recv_timeout(std::time::Duration::from_millis(50));
                    let task = match task {
                        Ok(task) => task,
                        Err(std::sync::mpsc::RecvTimeoutError::Timeout) => continue,
                        Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => return Ok(()),
                    };
                    let (temporary, file) = create_temp_file(".lastz")?;
                    drop(file);
                    crate::lastz_align::run_many_lastz(
                        lastz_bin,
                        &task.unit.target_spec,
                        query,
                        args.coverage,
                        args.identity,
                        &temporary.path.to_string_lossy(),
                    )?;
                    let has_results = temporary
                        .path
                        .metadata()
                        .with_context(|| {
                            format!("checking lastz output size {}", temporary.path.display())
                        })?
                        .len()
                        > 0;
                    event_sender.send(QueueEvent::UnitFinished(UnitResult {
                        genome_index: task.genome_index,
                        unit_index: task.unit_index,
                        output: has_results.then_some(temporary),
                    }))?;
                });
                if let Err(error) = result {
                    if !cancelled.swap(true, Ordering::AcqRel) {
                        let _ = event_sender.send(QueueEvent::Failed(error));
                    }
                }
            });
        }
        drop(event_sender);

        let mut pending = HashMap::new();
        let mut first_error = None;
        for event in event_receiver {
            match event {
                QueueEvent::GenomeStarted {
                    genome_index,
                    genome,
                } => {
                    pending.insert(
                        genome_index,
                        PendingGenome {
                            genome,
                            expected_units: None,
                            results: BTreeMap::new(),
                        },
                    );
                }
                QueueEvent::GenomeFinished {
                    genome_index,
                    unit_count,
                } => {
                    if let Some(state) = pending.get_mut(&genome_index) {
                        state.expected_units = Some(unit_count);
                    }
                    if first_error.is_none() {
                        if let Err(error) = finalize_completed_genome(
                            &mut pending,
                            genome_index,
                            conn,
                            probefile,
                            output_dir,
                        ) {
                            cancelled.store(true, Ordering::Release);
                            first_error = Some(error);
                        }
                    }
                }
                QueueEvent::UnitFinished(result) => {
                    if let Some(state) = pending.get_mut(&result.genome_index) {
                        state.results.insert(result.unit_index, result.output);
                    }
                    if first_error.is_none() {
                        if let Err(error) = finalize_completed_genome(
                            &mut pending,
                            result.genome_index,
                            conn,
                            probefile,
                            output_dir,
                        ) {
                            cancelled.store(true, Ordering::Release);
                            first_error = Some(error);
                        }
                    }
                }
                QueueEvent::Failed(error) => {
                    cancelled.store(true, Ordering::Release);
                    if first_error.is_none() {
                        first_error = Some(error);
                    }
                }
            }
        }
        if let Some(error) = first_error {
            return Err(error);
        }
        anyhow::ensure!(
            pending.is_empty(),
            "LASTZ queue ended with incomplete genomes"
        );
        Ok(())
    })
}

pub fn run(
    db: &Path,
    output_dir: &Path,
    probefile: &Path,
    args: &RunMultipleLastzsArgs,
) -> anyhow::Result<()> {
    anyhow::ensure!(args.cores > 0, "--cores must be greater than zero");
    anyhow::ensure!(
        probefile.is_file(),
        "{} is not a probe file",
        probefile.display()
    );
    std::fs::create_dir_all(output_dir)
        .with_context(|| format!("creating output directory {}", output_dir.display()))?;
    for g in args.chromolist.iter().chain(args.scaffoldlist.iter()) {
        let path = genome_path(&args.genome_base_path, args.no_dir, g);
        anyhow::ensure!(path.is_file(), "{} is not a file", path.display());
    }

    let cfg = phyluce_config::PhyluceConfig::load()?;
    let lastz_bin = cfg.get_user_path("binaries", "lastz")?;

    let conn = Connection::open(db)
        .with_context(|| format!("opening lastz results database {}", db.display()))?;
    if !args.append {
        conn.execute(
            "CREATE TABLE species (name TEXT PRIMARY KEY, description TEXT NULL, version TEXT NULL)",
            [],
        )?;
    }

    align_and_store_all(&conn, &lastz_bin, probefile, output_dir, args)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_twobit(records: &[(&str, &str)]) -> Vec<u8> {
        const SIGNATURE: u32 = 0x1A41_2743;
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&SIGNATURE.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&(records.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());

        let mut offset_positions = Vec::new();
        for (name, _) in records {
            bytes.push(name.len() as u8);
            bytes.extend_from_slice(name.as_bytes());
            offset_positions.push(bytes.len());
            bytes.extend_from_slice(&0u32.to_le_bytes());
        }

        for ((_, sequence), offset_position) in records.iter().zip(offset_positions) {
            let sequence_offset = bytes.len() as u32;
            bytes[offset_position..offset_position + 4]
                .copy_from_slice(&sequence_offset.to_le_bytes());
            bytes.extend_from_slice(&(sequence.len() as u32).to_le_bytes());
            bytes.extend_from_slice(&0u32.to_le_bytes());
            bytes.extend_from_slice(&0u32.to_le_bytes());
            bytes.extend_from_slice(&0u32.to_le_bytes());

            for chunk in sequence.as_bytes().chunks(4) {
                let mut packed = 0u8;
                for (index, base) in chunk.iter().enumerate() {
                    let code = match base {
                        b'T' => 0,
                        b'C' => 1,
                        b'A' => 2,
                        b'G' => 3,
                        _ => panic!("unsupported test base"),
                    };
                    packed |= code << (6 - 2 * index);
                }
                bytes.push(packed);
            }
        }
        bytes
    }

    #[test]
    fn clean_lastz_data_strips_percent_signs() {
        assert_eq!(clean_lastz_data("100.0%\t83.2%\tfoo"), "100.0\t83.2\tfoo");
    }

    #[test]
    fn genome_path_respects_no_dir() {
        assert_eq!(
            genome_path("/genomes", false, "gallus"),
            PathBuf::from("/genomes/gallus/gallus.2bit")
        );
        assert_eq!(
            genome_path("/genomes", true, "gallus"),
            PathBuf::from("/genomes/gallus.2bit")
        );
    }

    #[test]
    fn species_import_rolls_back_on_a_malformed_row() {
        let conn = Connection::open_in_memory().unwrap();
        conn.execute("CREATE TABLE species (name TEXT PRIMARY KEY)", [])
            .unwrap();
        let valid = vec!["1"; 19].join("\t");
        let input = format!("{valid}\nmalformed\n");
        assert!(store_species_rows(&conn, "taxon_a", &input).is_err());

        let table_count: i64 = conn
            .query_row(
                "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='taxon_a'",
                [],
                |row| row.get(0),
            )
            .unwrap();
        let species_count: i64 = conn
            .query_row("SELECT count(*) FROM species", [], |row| row.get(0))
            .unwrap();
        assert_eq!(table_count, 0);
        assert_eq!(species_count, 0);

        let path =
            std::env::temp_dir().join(format!("phyluce-malformed-lastz-{}", std::process::id()));
        std::fs::write(&path, input).unwrap();
        assert!(store_species_file(&conn, "taxon_b", &path).is_err());
        let table_count: i64 = conn
            .query_row(
                "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='taxon_b'",
                [],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(table_count, 0);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn parallel_unit_outputs_are_concatenated_in_target_order() {
        let units = (0..8)
            .map(|index| TargetUnit {
                target_spec: index.to_string(),
                _input: None,
            })
            .collect();
        let results = execute_units(units, 4, |target, output| {
            let index: u64 = target.parse()?;
            std::thread::sleep(std::time::Duration::from_millis((8 - index) * 2));
            std::fs::write(output, format!("{target}\n"))?;
            Ok(())
        })
        .unwrap();
        let (output, file) = create_temp_file(".combined.lastz").unwrap();
        drop(file);
        concatenate_results(&output.path, &results).unwrap();
        assert_eq!(
            std::fs::read_to_string(&output.path).unwrap(),
            "0\n1\n2\n3\n4\n5\n6\n7\n"
        );
    }

    #[test]
    fn multiple_genomes_share_one_parallel_queue() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let active = AtomicUsize::new(0);
        let peak = AtomicUsize::new(0);
        let genomes = ["genome-a", "genome-b"]
            .into_iter()
            .map(|genome| GenomeUnits {
                genome: genome.to_string(),
                units: vec![TargetUnit {
                    target_spec: genome.to_string(),
                    _input: None,
                }],
            })
            .collect();
        let grouped = execute_genome_queue(genomes, 2, |target, output| {
            let now = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
            std::thread::sleep(std::time::Duration::from_millis(30));
            std::fs::write(output, target)?;
            active.fetch_sub(1, Ordering::SeqCst);
            Ok(())
        })
        .unwrap();

        assert_eq!(peak.load(Ordering::SeqCst), 2);
        assert_eq!(grouped[0].0, "genome-a");
        assert_eq!(grouped[1].0, "genome-b");
        assert_eq!(grouped[0].1.len(), 1);
        assert_eq!(grouped[1].1.len(), 1);
    }

    #[test]
    fn twobit_targets_split_by_sequence_and_scaffold_size() {
        let (temporary, file) = create_temp_file(".2bit").unwrap();
        drop(file);
        std::fs::write(
            &temporary.path,
            build_twobit(&[("chr2", "TTTT"), ("chr1", "ACGT")]),
        )
        .unwrap();

        let chromosomes = chromosome_units(&temporary.path).unwrap();
        assert!(chromosomes[0].target_spec.ends_with("/chr1"));
        assert!(chromosomes[1].target_spec.ends_with("/chr2"));

        let scaffolds = scaffold_units(&temporary.path, 3).unwrap();
        assert_eq!(scaffolds.len(), 2);
        assert_eq!(
            std::fs::read_to_string(&scaffolds[0].target_spec).unwrap(),
            ">chr1\nACGT\n"
        );
        assert_eq!(
            std::fs::read_to_string(&scaffolds[1].target_spec).unwrap(),
            ">chr2\nTTTT\n"
        );
    }

    #[cfg(unix)]
    #[test]
    fn fake_lastz_runs_through_merge_clean_and_sqlite_import() {
        use std::os::unix::fs::PermissionsExt;

        let root =
            std::env::temp_dir().join(format!("phyluce-lastz-integration-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        std::fs::write(
            root.join("test.2bit"),
            build_twobit(&[("chr2", "TTTT"), ("chr1", "ACGT")]),
        )
        .unwrap();
        let probes = root.join("probes.fasta");
        std::fs::write(&probes, ">uce-1_p1\nACGT\n").unwrap();
        let fake_lastz = root.join("lastz");
        let active_lock = root.join("lastz-active");
        let overlap_marker = root.join("lastz-overlap");
        std::fs::write(
            &fake_lastz,
            format!(
                "#!/bin/sh\nout=''\nfor arg in \"$@\"; do\n  case \"$arg\" in --output=*) out=${{arg#--output=}};; esac\ndone\nif mkdir \"{}\" 2>/dev/null; then\n  sleep 0.3\n  rmdir \"{}\"\nelse\n  touch \"{}\"\n  sleep 0.3\nfi\nprintf '1\\ttarget\\t+\\t0\\t4\\t4\\tuce-1_p1\\t+\\t0\\t4\\t4\\t0/4\\t4M\\t4/4\\t100.0%%\\t4/4\\t100.0%%\\t4/4\\t100.0%%\\n' > \"$out\"\n",
                active_lock.display(),
                active_lock.display(),
                overlap_marker.display()
            ),
        )
        .unwrap();
        let mut permissions = std::fs::metadata(&fake_lastz).unwrap().permissions();
        permissions.set_mode(0o700);
        std::fs::set_permissions(&fake_lastz, permissions).unwrap();

        let output = root.join("output");
        std::fs::create_dir_all(&output).unwrap();
        let mut args = RunMultipleLastzsArgs {
            chromolist: vec!["test".to_string()],
            scaffoldlist: Vec::new(),
            append: false,
            no_dir: true,
            genome_base_path: root.to_string_lossy().to_string(),
            coverage: 83.0,
            identity: 92.5,
            cores: 1,
        };
        let serial_cleaned = align_genome(
            &fake_lastz.to_string_lossy(),
            "test",
            false,
            &probes,
            &output,
            &args,
        )
        .unwrap();
        assert_eq!(serial_cleaned.lines().count(), 2);
        assert!(
            !overlap_marker.exists(),
            "single-core LASTZ processes unexpectedly overlapped"
        );

        args.cores = 2;
        let cleaned = align_genome(
            &fake_lastz.to_string_lossy(),
            "test",
            false,
            &probes,
            &output,
            &args,
        )
        .unwrap();
        assert!(
            overlap_marker.is_file(),
            "fake LASTZ processes did not overlap"
        );
        assert_eq!(cleaned.lines().count(), 2);
        assert!(!cleaned.contains('%'));

        let conn = Connection::open_in_memory().unwrap();
        conn.execute("CREATE TABLE species (name TEXT PRIMARY KEY)", [])
            .unwrap();
        store_species_rows(&conn, "test", &cleaned).unwrap();
        let row_count: i64 = conn
            .query_row("SELECT count(*) FROM test", [], |row| row.get(0))
            .unwrap();
        assert_eq!(row_count, 2);
        assert!(output.join("probes.fasta_v_test.lastz.clean").is_file());
        assert!(!output.join("probes.fasta_v_test.lastz").exists());

        std::fs::write(root.join("one.2bit"), build_twobit(&[("chr1", "ACGT")])).unwrap();
        std::fs::write(root.join("two.2bit"), build_twobit(&[("chr1", "TTTT")])).unwrap();
        let _ = std::fs::remove_file(&overlap_marker);
        args.chromolist = vec!["one".to_string(), "two".to_string()];
        let conn = Connection::open_in_memory().unwrap();
        conn.execute("CREATE TABLE species (name TEXT PRIMARY KEY)", [])
            .unwrap();
        align_and_store_all(
            &conn,
            &fake_lastz.to_string_lossy(),
            &probes,
            &output,
            &args,
        )
        .unwrap();
        assert!(
            overlap_marker.is_file(),
            "different genomes did not overlap in the bounded global queue"
        );
        let species_count: i64 = conn
            .query_row("SELECT count(*) FROM species", [], |row| row.get(0))
            .unwrap();
        assert_eq!(species_count, 2);
        assert!(output.join("probes.fasta_v_one.lastz.clean").is_file());
        assert!(output.join("probes.fasta_v_two.lastz.clean").is_file());

        std::fs::remove_dir_all(root).unwrap();
    }
}
