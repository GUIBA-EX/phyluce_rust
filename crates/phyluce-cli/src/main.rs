use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use mimalloc::MiMalloc;

/// Hot paths across this CLI (contig/probe matching, alignment
/// concatenation, PAML partitioning) do a lot of small, short-lived
/// `String`/`Vec`/`HashMap`-entry allocations; `mimalloc` is consistently
/// faster than the system allocator for that pattern (same rationale as
/// `probebwa`'s `src/main.rs`, which measured this directly). Opting in
/// only affects this binary, not downstream consumers of the library
/// crates.
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use anyhow::Context;
use clap::{ArgAction, Parser, Subcommand};
use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;
use phyluce_io::fastq::fastq_lengths;
use phyluce_io::lastz::read_lastz;
use phyluce_io::{fasta_lengths, read_fasta, validate_fasta};
use tracing::level_filters::LevelFilter;

/// Preserve CLI stdout while also recording operational messages in the
/// optional tracing log. Commands use this rather than writing only to stdout.
#[macro_export]
macro_rules! cli_info {
    () => {{
        use std::io::Write as _;
        writeln!(std::io::stdout())?;
        tracing::info!(message = "");
    }};
    ($($arg:tt)*) => {{
        let message = format!($($arg)*);
        use std::io::Write as _;
        writeln!(std::io::stdout(), "{message}")?;
        tracing::info!(message = %message);
    }};
}

/// Preserve CLI stderr while also recording warnings in the optional log.
#[macro_export]
macro_rules! cli_warn {
    () => {{
        use std::io::Write as _;
        writeln!(std::io::stderr())?;
        tracing::warn!(message = "");
    }};
    ($($arg:tt)*) => {{
        let message = format!($($arg)*);
        use std::io::Write as _;
        writeln!(std::io::stderr(), "{message}")?;
        tracing::warn!(message = %message);
    }};
}

mod align_summary_cmd;
mod assemblo_abyss_cmd;
mod assemblo_spades_cmd;
mod assemblo_velvet_cmd;
mod bootrep_support_cmd;
mod bootstrap_count_cmd;
mod chunk_fasta_cmd;
mod combine_reads_cmd;
mod concatenate_cmd;
mod conf;
mod convert_cmd;
mod convert_degen_bases_cmd;
mod easy_lastz_cmd;
mod easy_stampy_cmd;
mod explode_alignments_cmd;
mod explode_cmd;
mod extract_contigs_to_barcodes_cmd;
mod extract_taxa_cmd;
mod extract_taxon_fasta_cmd;
mod filter_alignments_cmd;
mod filter_bed_cmd;
mod format_paml_cmd;
mod gblocks_cmd;
mod genome_sequences_from_bed_cmd;
mod get_fastas_cmd;
mod get_trimmed_cmd;
mod incomplete_matrix_estimates_cmd;
mod informative_sites_cmd;
mod lastz_align;
mod match_contigs;
mod match_contigs_to_barcodes_cmd;
mod match_counts_cmd;
mod merge_gzip_cmd;
mod merge_nextseq_cmd;
mod min_taxa_filter_cmd;
mod missing_data_cmd;
mod move_align_cmd;
mod multi_fasta_table_cmd;
mod multi_merge_table_cmd;
mod ncbi_prep_cmd;
mod output_path;
mod parallel;
mod probe_bed_from_lastz_cmd;
mod probebwa_align;
mod randomly_sample_concat_cmd;
mod reconstruct_uce_from_probe_cmd;
mod reduce_raxml_cmd;
mod remove_duplicate_hits_cmd;
mod remove_empty_taxa_cmd;
mod remove_locus_name_cmd;
mod remove_overlapping_probes_cmd;
mod rename_tree_leaves_cmd;
mod replace_links_cmd;
mod run_multiple_lastzs_sqlite_cmd;
mod ry_recode_cmd;
mod sample_reads_cmd;
mod screen_alignments_cmd;
mod screen_probes_dupes_cmd;
mod screened_loci_proximity_cmd;
mod seqcap_align_cmd;
mod slice_sequence_from_genomes_cmd;
mod smilogram_cmd;
mod sort_bootstraps_cmd;
mod split_concat_cmd;
mod stats;
mod strip_masked_loci_cmd;
mod subsets_tiled_probes_cmd;
mod taxon_locus_counts_cmd;
mod tiled_probe_from_multiple_inputs_cmd;
mod tiled_probes_cmd;
mod tree_counts_cmd;
mod trimal_cmd;
mod unmix_fasta_cmd;
mod workflow_cmd;

#[derive(Parser)]
#[command(
    name = "phyluce",
    version,
    about = "phyluce: software for UCE (and general) phylogenomics -- Rust CLI (early preview)"
)]
struct Cli {
    /// Logging verbosity for `--log-path` output.
    #[arg(long, global = true, default_value = "INFO", value_parser = ["INFO", "WARN", "CRITICAL"])]
    verbosity: String,
    /// Directory in which to write a `<command>.log` file.
    #[arg(long, global = true)]
    log_path: Option<PathBuf>,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Inspect and validate phyluce.conf resolution ($CONDA/$WORKFLOWS, user overrides).
    Config {
        #[command(subcommand)]
        action: ConfigAction,
    },
    /// FASTA/FASTQ I/O utilities.
    Io {
        #[command(subcommand)]
        action: IoAction,
    },
    /// Assembly-domain commands (mirrors `bin/assembly/phyluce_assembly_*`).
    Assembly {
        #[command(subcommand)]
        action: AssemblyAction,
    },
    /// Run/verify an external binary resolved from phyluce.conf.
    External {
        #[command(subcommand)]
        action: ExternalAction,
    },
    /// Utility-domain commands (mirrors `bin/utilities/phyluce_utilities_*`).
    Utilities {
        #[command(subcommand)]
        action: UtilitiesAction,
    },
    /// Alignment-domain commands (mirrors `bin/align/phyluce_align_*`).
    Align {
        #[command(subcommand)]
        action: AlignAction,
    },
    /// NCBI-domain commands (mirrors `bin/ncbi/phyluce_ncbi_*`).
    Ncbi {
        #[command(subcommand)]
        action: NcbiAction,
    },
    /// Genetrees-domain commands (mirrors `bin/genetrees/phyluce_genetrees_*`).
    Genetrees {
        #[command(subcommand)]
        action: GenetreesAction,
    },
    /// Equivalent to `phyluce_workflow`: run a Snakemake workflow. Untested
    /// (Snakemake isn't installed in this environment).
    Workflow {
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        workflow: String,
        #[arg(long, default_value_t = 1)]
        cores: u32,
        #[arg(long, default_value_t = false)]
        dryrun: bool,
    },
    /// Probe-domain commands (mirrors `bin/probes/phyluce_probe_*`).
    Probe {
        #[command(subcommand)]
        action: ProbeAction,
    },
}

#[derive(Subcommand)]
enum ProbeAction {
    /// Equivalent to `phyluce_probe_remove_overlapping_probes_given_config`.
    RemoveOverlappingProbesGivenConfig {
        #[arg(long)]
        probes: PathBuf,
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
    /// Equivalent to `phyluce_probe_get_probe_bed_from_lastz_files`.
    GetProbeBedFromLastzFiles {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
    /// Equivalent to `phyluce_probe_get_locus_bed_from_lastz_files`.
    GetLocusBedFromLastzFiles {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = r"^(uce-\d+)(?:_p\d+.*)")]
        regex: String,
    },
    /// Equivalent to `phyluce_probe_get_subsets_of_tiled_probes`.
    GetSubsetsOfTiledProbes {
        #[arg(long)]
        probes: PathBuf,
        #[arg(long, num_args = 1.., value_delimiter = ' ')]
        taxa: Vec<String>,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = r"^(uce-\d+)(?:_p\d+.*)")]
        regex: String,
    },
    /// Equivalent to `phyluce_probe_get_multi_fasta_table`.
    GetMultiFastaTable {
        #[arg(long)]
        fastas: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        base_taxon: String,
    },
    /// Equivalent to `phyluce_probe_query_multi_fasta_table`.
    QueryMultiFastaTable {
        #[arg(long)]
        db: PathBuf,
        #[arg(long)]
        base_taxon: String,
        #[arg(long)]
        specific_counts: Option<usize>,
        #[arg(long)]
        output: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_probe_get_multi_merge_table`.
    GetMultiMergeTable {
        #[arg(long)]
        conf: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        base_taxon: String,
    },
    /// Equivalent to `phyluce_probe_query_multi_merge_table`.
    QueryMultiMergeTable {
        #[arg(long)]
        db: PathBuf,
        #[arg(long)]
        base_taxon: String,
        #[arg(long)]
        specific_counts: Option<usize>,
        #[arg(long)]
        output: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_probe_get_screened_loci_by_proximity`. Ties
    /// within a cluster are broken deterministically (lowest locus id kept)
    /// rather than via Python's `random.choice` -- see
    /// `screened_loci_proximity_cmd` docs.
    GetScreenedLociByProximity {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 10000)]
        distance: i64,
        #[arg(long, default_value = r"^uce-(\d+)(?:_p\d+.*)")]
        regex: String,
    },
    /// Equivalent to `phyluce_probe_remove_duplicate_hits_from_probes_using_lastz`.
    RemoveDuplicateHitsFromProbesUsingLastz {
        #[arg(long)]
        fasta: PathBuf,
        #[arg(long)]
        lastz: PathBuf,
        #[arg(long)]
        probe_prefix: String,
        #[arg(long, default_value = r"^({}\d+)(?:_p\d+.*)")]
        probe_regex: String,
        #[arg(long)]
        probe_bed: Option<PathBuf>,
        #[arg(long)]
        locus_bed: Option<PathBuf>,
        #[arg(long, default_value_t = false)]
        long: bool,
    },
    /// Equivalent to `phyluce_probe_get_tiled_probe_from_multiple_inputs`.
    /// `--two-probes` breaks an odd-coordinate-count tie deterministically
    /// (always `+1`) instead of via Python's `random.choice([1, -1])` --
    /// see `tiled_probe_from_multiple_inputs_cmd` docs.
    GetTiledProbeFromMultipleInputs {
        #[arg(long)]
        fastas: PathBuf,
        #[arg(long)]
        multi_fasta_output: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        probe_prefix: String,
        #[arg(long)]
        designer: String,
        #[arg(long)]
        design: String,
        #[arg(long, default_value_t = 120)]
        probe_length: usize,
        #[arg(long, default_value_t = 2.0)]
        tiling_density: f64,
        #[arg(long)]
        masking: Option<f64>,
        #[arg(
            long = "do-not-remove-ambiguous",
            default_value_t = true,
            action = ArgAction::SetFalse
        )]
        remove_ambiguous: bool,
        #[arg(long, default_value_t = false)]
        remove_gc: bool,
        #[arg(long, default_value_t = 1)]
        start_index: usize,
        #[arg(long, default_value_t = false)]
        two_probes: bool,
    },
    /// Equivalent to `phyluce_probe_get_tiled_probes`. `--two-probes`
    /// breaks an odd-coordinate-count tie deterministically (always `+1`)
    /// instead of via Python's `random.choice([1, -1])` -- see
    /// `tiled_probes_cmd` docs.
    GetTiledProbes {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        probe_prefix: String,
        #[arg(long)]
        designer: String,
        #[arg(long)]
        design: String,
        #[arg(long, default_value_t = 120)]
        probe_length: usize,
        #[arg(long, default_value_t = 2.0)]
        tiling_density: f64,
        #[arg(long, default_value = "middle")]
        overlap: String,
        #[arg(long)]
        probe_bed: Option<PathBuf>,
        #[arg(long)]
        locus_bed: Option<PathBuf>,
        #[arg(long)]
        masking: Option<f64>,
        #[arg(
            long = "do-not-remove-ambiguous",
            default_value_t = true,
            action = ArgAction::SetFalse
        )]
        remove_ambiguous: bool,
        #[arg(long, default_value_t = false)]
        remove_gc: bool,
        #[arg(long, default_value_t = 0)]
        start_index: usize,
        #[arg(long, default_value_t = false)]
        two_probes: bool,
    },
    /// Equivalent to `phyluce_probe_reconstruct_uce_from_probe`. Multi-probe
    /// loci use MAFFT by default; MUSCLE 3 is available as an explicit legacy
    /// compatibility path.
    ReconstructUceFromProbe {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
        /// Explicitly use the legacy MUSCLE 3 `-clwstrict` path.
        #[arg(long)]
        muscle_binary: Option<String>,
        /// Override the default MAFFT executable.
        #[arg(long)]
        mafft_binary: Option<String>,
    },
    /// Equivalent to `phyluce_probe_get_genome_sequences_from_bed`. Reads
    /// genomes from UCSC `.2bit` files via a hand-rolled parser (see
    /// `phyluce_io::twobit`) rather than `bx.seq.twobit`.
    GetGenomeSequencesFromBed {
        #[arg(long)]
        bed: PathBuf,
        #[arg(long)]
        twobit: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 0.25)]
        filter_mask: f64,
        #[arg(long, default_value_t = 0)]
        max_n: usize,
        #[arg(long)]
        buffer_to: Option<i64>,
    },
    /// Equivalent to `phyluce_probe_strip_masked_loci_from_set`.
    StripMaskedLociFromSet {
        #[arg(long)]
        bed: PathBuf,
        #[arg(long)]
        twobit: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        filter_mask: Option<f64>,
        #[arg(long, default_value_t = 0)]
        max_n: usize,
        #[arg(long, default_value_t = 0)]
        min_length: i64,
    },
    /// Equivalent to `phyluce_probe_slice_sequence_from_genomes`. `--conf`
    /// must have a `[chromos]` and/or `[scaffolds]` section mapping short
    /// genome names to `.2bit` file paths.
    SliceSequenceFromGenomes {
        #[arg(long)]
        conf: PathBuf,
        #[arg(long)]
        lastz: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        name_pattern: Option<String>,
        #[arg(long, default_value = "uce-")]
        probe_prefix: String,
        #[arg(long, default_value = r"^({}\d+)(?:_p\d+.*)")]
        probe_regex: String,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        exclude: Vec<String>,
        #[arg(long, alias = "contig_orient", default_value_t = false)]
        contig_orient: bool,
        #[arg(long, required_unless_present = "probes", conflicts_with = "probes")]
        flank: Option<i64>,
        #[arg(long, required_unless_present = "flank", conflicts_with = "flank")]
        probes: Option<i64>,
    },
    /// Equivalent to `phyluce_probe_easy_lastz`.
    EasyLastz {
        #[arg(long)]
        target: PathBuf,
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 92.5)]
        identity: f64,
        #[arg(long, default_value_t = 83.0, conflicts_with = "min_match")]
        coverage: f64,
        #[arg(long, alias = "min_match", conflicts_with = "coverage")]
        min_match: Option<i64>,
    },
    /// Runs the stampy_ genome-prep + mapping workflow from
    /// `docs/tutorials/tutorial-4.rst` using `probebwa`
    /// (<https://github.com/GUIBA-EX/probebwa>), a stampy-compatible Rust
    /// mapper, in place of `stampy.py`. Chains `probebwa build-genome`,
    /// `build-hash`, and `map` in sequence; pass `--bam` to write BAM
    /// directly instead of piping SAM through `samtools view` by hand.
    /// A `build-genome`/`build-hash` step is skipped when its output file
    /// (`<index-prefix>.stidx`/`.sthash`) already exists, so mapping many
    /// samples against the same reference only pays the indexing cost
    /// once; pass `--force-rebuild-index` to rebuild anyway (e.g. after
    /// changing `--genome-files`).
    EasyStampy {
        #[arg(long)]
        species: String,
        #[arg(long)]
        assembly: String,
        #[arg(long, num_args = 1..)]
        genome_files: Vec<String>,
        /// Prefix for the `.stidx`/`.sthash` files (mirrors stampy's `-G`/`-H`).
        #[arg(long)]
        index_prefix: PathBuf,
        /// Read file(s) to map: one for single-ended, two (mate1 mate2) for paired-end.
        #[arg(long, num_args = 1..=2)]
        reads: Vec<String>,
        #[arg(long, default_value_t = 0.05)]
        substitution_rate: f64,
        #[arg(long, default_value_t = 1)]
        threads: usize,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = false)]
        bam: bool,
        /// Rebuild `.stidx`/`.sthash` even if they already exist.
        #[arg(long, default_value_t = false)]
        force_rebuild_index: bool,
    },
    /// Equivalent to `phyluce_probe_run_multiple_lastzs_sqlite`. Splits
    /// chromosome targets by sequence and scaffold targets into roughly
    /// 10 Mbp chunks, then runs LASTZ with at most `--cores` workers.
    RunMultipleLastzsSqlite {
        #[arg(long)]
        db: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        probefile: PathBuf,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        chromolist: Vec<String>,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        scaffoldlist: Vec<String>,
        #[arg(long, default_value_t = false)]
        append: bool,
        #[arg(long, default_value_t = false)]
        no_dir: bool,
        #[arg(long, default_value_t = 1)]
        cores: usize,
        #[arg(long)]
        genome_base_path: String,
        #[arg(long, default_value_t = 83.0)]
        coverage: f64,
        #[arg(long, default_value_t = 92.5)]
        identity: f64,
    },
}

#[derive(Subcommand)]
enum GenetreesAction {
    /// Equivalent to `phyluce_genetrees_rename_tree_leaves`.
    RenameTreeLeaves {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        section: String,
        #[arg(long, default_value = "left:right")]
        order: String,
        #[arg(long)]
        reroot: Option<String>,
    },
    /// Equivalent to `phyluce_genetrees_get_tree_counts`. Topology grouping
    /// is rooting-invariant rather than physically rerooting each tree
    /// (see `tree_counts_cmd` docs) -- counts/grouping match, but the
    /// printed representative Newick strings keep their as-parsed
    /// rooting.
    GetTreeCounts {
        #[arg(long)]
        trees: PathBuf,
        #[arg(long)]
        locus_support_output: PathBuf,
        #[arg(long)]
        root: String,
        #[arg(long)]
        extension: String,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        exclude: Vec<String>,
    },
    /// Equivalent to `phyluce_genetrees_get_mean_bootrep_support`. Always
    /// writes `outfile.csv` in the current directory (matching the legacy
    /// script -- there is no `--output` option).
    GetMeanBootrepSupport {
        #[arg(long)]
        trees: PathBuf,
        #[arg(long)]
        config: PathBuf,
    },
    /// Equivalent to `phyluce_genetrees_generate_multilocus_bootstrap_count`.
    /// Uses a plain-text replicate format instead of Python `pickle` --
    /// see `bootstrap_count_cmd` docs.
    GenerateMultilocusBootstrapCount {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long, alias = "bootstrap_replicates")]
        bootstrap_replicates: PathBuf,
        #[arg(long, default_value = "")]
        directory: String,
        #[arg(long, alias = "bootstrap_counts")]
        bootstrap_counts: PathBuf,
        #[arg(long, default_value_t = 100)]
        bootreps: usize,
    },
    /// Equivalent to `phyluce_genetrees_sort_multilocus_bootstraps`.
    SortMultilocusBootstraps {
        #[arg(long)]
        input: PathBuf,
        #[arg(long, alias = "bootstrap_replicates")]
        bootstrap_replicates: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
}

#[derive(Subcommand)]
enum NcbiAction {
    /// Equivalent to `phyluce_ncbi_chunk_fasta_for_ncbi`.
    ChunkFastaForNcbi {
        #[arg(long)]
        input: PathBuf,
        #[arg(long, default_value_t = 10000)]
        chunk_size: usize,
        #[arg(long, default_value = "split")]
        output_prefix: String,
        #[arg(long, default_value = "fsa")]
        output_suffix: String,
    },
    /// Equivalent to `phyluce_ncbi_prep_uce_align_files_for_ncbi`.
    /// Untested: the legacy Python command currently crashes on import in
    /// this environment (`Bio.Alphabet` was removed from modern
    /// Biopython) -- see `ncbi_prep_cmd` docs.
    PrepUceAlignFilesForNcbi {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        conf: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
    },
}

#[derive(Subcommand)]
enum ConfigAction {
    /// Print the resolved config file paths and all sections/keys.
    Inspect,
    /// Resolve a single `[program] binary:path` entry, expanding placeholders.
    Which {
        #[arg(long)]
        program: String,
        #[arg(long)]
        binary: String,
    },
}

#[derive(Subcommand)]
enum IoAction {
    /// Validate that a file is well-formed FASTA (headers, non-empty
    /// sequences, IUPAC alphabet). Exits non-zero if any issues are found.
    ValidateFasta {
        #[arg(long)]
        input: PathBuf,
    },
}

#[derive(Subcommand)]
#[allow(clippy::enum_variant_names)] // mirrors legacy `phyluce_assembly_get_*` command names
enum AssemblyAction {
    /// Equivalent to `phyluce_assembly_get_fasta_lengths`: summary length
    /// statistics for a FASTA (or FASTA.gz) file of contigs/reads.
    GetFastaLengths {
        #[arg(long)]
        input: PathBuf,
        /// Emit CSV instead of the human-readable report.
        #[arg(long, default_value_t = false)]
        csv: bool,
    },
    /// Equivalent to `phyluce_assembly_get_fastq_lengths`: summary length
    /// statistics across all `*.fastq*` files in a directory.
    GetFastqLengths {
        #[arg(long)]
        input: PathBuf,
        /// Emit CSV instead of the human-readable report.
        #[arg(long, default_value_t = false)]
        csv: bool,
    },
    /// Equivalent to `phyluce_assembly_get_bed_from_lastz`: write a BED file
    /// from LASTZ general-format alignment output.
    GetBedFromLastz {
        #[arg(long)]
        lastz: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 90.0)]
        identity: f64,
        #[arg(long, default_value_t = 90.0)]
        continuity: f64,
        #[arg(long, default_value_t = false)]
        long_format: bool,
        /// INI file with sections listing loci to keep (requires --sections).
        #[arg(long)]
        conf: Option<PathBuf>,
        #[arg(long, num_args = 1.., value_delimiter = ' ')]
        sections: Option<Vec<String>>,
    },
    /// Equivalent to `phyluce_assembly_match_contigs_to_probes`: LASTZ-align
    /// each taxon's contigs against the probe set and store results in
    /// `probe.matches.sqlite`.
    MatchContigsToProbes {
        #[arg(long)]
        contigs: PathBuf,
        #[arg(long)]
        probes: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 80)]
        min_coverage: u32,
        #[arg(long, default_value_t = 80)]
        min_identity: u32,
        /// Path to self-to-self LASTZ results for baits, to remove
        /// potential duplicate probes.
        #[arg(long)]
        dupefile: Option<PathBuf>,
        #[arg(long, default_value = r"^(uce-\d+)(?:_p\d+.*)")]
        regex: String,
        /// Optional output file listing loci that appear to be duplicates.
        #[arg(long)]
        keep_duplicates: Option<PathBuf>,
        /// Optional CSV summary output file.
        #[arg(long)]
        csv: Option<PathBuf>,
        /// Reuse pre-existing `<output>/<contig-stem>.lastz` files instead
        /// of invoking the `lastz` binary (Rust-only addition, useful in
        /// environments without lastz installed, e.g. CI replaying
        /// precomputed fixtures).
        #[arg(long, default_value_t = false)]
        skip_alignment: bool,
        /// Remove --output if it already exists instead of erroring
        /// (Rust-only addition; the legacy CLI prompts interactively).
        #[arg(long, default_value_t = false)]
        force: bool,
    },
    /// Equivalent to `phyluce_assembly_get_match_counts`: generate a matrix
    /// config or optimize complete-matrix taxon membership.
    GetMatchCounts {
        #[arg(long)]
        locus_db: PathBuf,
        #[arg(long)]
        taxon_list_config: PathBuf,
        #[arg(long)]
        taxon_group: String,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = false)]
        incomplete_matrix: bool,
        /// Find the taxon subset with the most shared UCE loci. Without
        /// `--random`, enumerate every group size and write a report.
        #[arg(long, default_value_t = false)]
        optimize: bool,
        /// Estimate the best `--sample-size` taxon subset by repeated
        /// sampling instead of exhaustive enumeration.
        #[arg(long, default_value_t = false, requires = "optimize")]
        random: bool,
        /// Number of random optimization iterations.
        #[arg(long, default_value_t = 10)]
        samples: usize,
        /// Number of taxa retained by each random optimization iteration.
        #[arg(long, default_value_t = 10)]
        sample_size: usize,
        /// Suppress writing taxa and locus names for normal/random modes.
        #[arg(long, default_value_t = false)]
        silent: bool,
        /// Write `sample_size,locus_count` for every random iteration instead
        /// of the best matrix config.
        #[arg(long, default_value_t = false, requires = "random")]
        keep_counts: bool,
        /// Reproducible random seed. A generated seed is reported when omitted.
        #[arg(long, requires = "random")]
        seed: Option<u64>,
        /// Maximum concurrent group-size searches during exhaustive optimization.
        #[arg(long, default_value_t = 6)]
        cores: usize,
        /// An additional SQLite database of probe matches, ATTACHed as
        /// `extended` (referenced via a trailing `*` on organism names).
        #[arg(long)]
        extend_locus_db: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_assembly_get_fastas_from_match_counts`:
    /// extract each taxon's matched UCE loci from its contigs into a
    /// single monolithic renamed/reoriented FASTA.
    GetFastasFromMatchCounts {
        #[arg(long)]
        contigs: PathBuf,
        #[arg(long)]
        locus_db: PathBuf,
        #[arg(long)]
        match_count_output: PathBuf,
        /// Path to write missing-locus info for an incomplete matrix.
        /// Passing this switches to incomplete-matrix mode (mirrors the
        /// legacy CLI's overloaded `--incomplete-matrix <path>` flag).
        #[arg(long)]
        incomplete_matrix: Option<PathBuf>,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        extend_locus_db: Option<PathBuf>,
        #[arg(long)]
        extend_locus_contigs: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_assembly_explode_get_fastas_file`: split a
    /// monolithic FASTA (from `get-fastas-from-match-counts`) into one file
    /// per locus (or, with `--by-taxon`, one file per taxon).
    ExplodeGetFastasFile {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = false)]
        by_taxon: bool,
        #[arg(long, default_value = "_")]
        split_char: String,
    },
    /// Equivalent to `phyluce_assembly_assemblo_spades`. Untested against
    /// a live SPAdes binary in this environment.
    AssembloSpades {
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 1)]
        cores: u32,
        #[arg(long, default_value_t = 8)]
        memory: u32,
        #[arg(long, default_value = "")]
        subfolder: String,
        #[arg(long, default_value_t = false)]
        no_clean: bool,
        #[arg(long)]
        config: Option<PathBuf>,
        #[arg(long)]
        dir: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_assembly_assemblo_velvet`. Untested against
    /// live velveth/velvetg binaries in this environment.
    AssembloVelvet {
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 31)]
        kmer: u32,
        #[arg(long, default_value = "")]
        subfolder: String,
        #[arg(long, default_value_t = false)]
        clean: bool,
        #[arg(long)]
        config: Option<PathBuf>,
        #[arg(long)]
        dir: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_assembly_assemblo_abyss`. Untested against
    /// live abyss-pe/abyss-se binaries in this environment.
    AssembloAbyss {
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 31)]
        kmer: u32,
        #[arg(long, default_value_t = 1)]
        cores: u32,
        #[arg(long, default_value = "")]
        subfolder: String,
        #[arg(long, default_value_t = false)]
        clean: bool,
        #[arg(long, default_value_t = false)]
        abyss_se: bool,
        #[arg(long)]
        config: Option<PathBuf>,
        #[arg(long)]
        dir: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_assembly_screen_probes_for_dupes`. The
    /// Python original is Python-2-only syntax and can't run under
    /// Python 3 at all -- see `screen_probes_dupes_cmd` docs.
    ScreenProbesForDupes {
        #[arg(long)]
        lastz: PathBuf,
    },
    /// Equivalent to `phyluce_assembly_extract_contigs_to_barcodes`.
    ExtractContigsToBarcodes {
        #[arg(long)]
        contigs: PathBuf,
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
    /// Equivalent to `phyluce_assembly_match_contigs_to_barcodes`.
    /// Untested (`lastz` not installed); pass `--no-bold` because the BOLD
    /// systems web-lookup step is not reproduced -- see
    /// `match_contigs_to_barcodes_cmd` docs.
    MatchContigsToBarcodes {
        #[arg(long)]
        contigs: PathBuf,
        #[arg(long)]
        barcodes: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = false)]
        no_bold: bool,
        #[arg(long, default_value = "COX1_SPECIES")]
        database: String,
    },
}

#[derive(Subcommand)]
enum UtilitiesAction {
    /// Equivalent to `phyluce_utilities_get_bed_from_fasta`: derive a BED
    /// file from FASTA headers of the form
    /// `id|contig:NAME|coords:BEGIN-END|locus:LOCUS`.
    GetBedFromFasta {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "")]
        locus_prefix: String,
    },
    /// Equivalent to `phyluce_utilities_filter_bed_by_fasta`.
    FilterBedByFasta {
        #[arg(long)]
        bed: PathBuf,
        #[arg(long)]
        fasta: PathBuf,
        #[arg(long)]
        output: Option<PathBuf>,
    },
    /// Equivalent to `phyluce_utilities_replace_many_links`.
    ReplaceManyLinks {
        #[arg(long)]
        indir: PathBuf,
        #[arg(long)]
        oldpath: String,
        #[arg(long)]
        newpath: String,
        #[arg(long)]
        outdir: PathBuf,
    },
    /// Equivalent to `phyluce_utilities_combine_reads`.
    CombineReads {
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "")]
        subfolder: String,
    },
    /// Equivalent to `phyluce_utilities_merge_multiple_gzip_files`.
    MergeMultipleGzipFiles {
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "samples")]
        section: String,
        #[arg(long, default_value_t = false)]
        trimmed: bool,
    },
    /// Equivalent to `phyluce_utilities_merge_next_seq_gzip_files`.
    MergeNextSeqGzipFiles {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "samples")]
        section: String,
        #[arg(long, default_value_t = false)]
        se: bool,
    },
    /// Equivalent to `phyluce_utilities_unmix_fasta_reads`.
    UnmixFastaReads {
        #[arg(long)]
        mixed_reads: PathBuf,
        #[arg(long)]
        singleton_reads: Option<PathBuf>,
        #[arg(long)]
        out_r1: PathBuf,
        #[arg(long)]
        out_r2: PathBuf,
        #[arg(long)]
        out_r_singleton: PathBuf,
        #[arg(long, default_value_t = false)]
        new_style: bool,
    },
    /// Equivalent to `phyluce_utilities_sample_reads_from_files`. Shells out
    /// to `seqkit sample` rather than the Python original's `seqtk sample`
    /// -- see `sample_reads_cmd` docs for why, and what's not
    /// byte-for-byte compatible about it.
    SampleReadsFromFiles {
        #[arg(long)]
        conf: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
}

#[derive(Subcommand)]
enum AlignAction {
    /// Equivalent to `phyluce_align_get_trimmed_alignments_from_untrimmed`:
    /// apply the native phyluce 3-stage edge-trimming algorithm to a
    /// directory of existing (fasta) alignments, writing NEXUS output.
    GetTrimmedAlignmentsFromUntrimmed {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 20)]
        window: usize,
        #[arg(long, default_value_t = 0.65)]
        proportion: f64,
        #[arg(long, default_value_t = 0.65)]
        threshold: f64,
        #[arg(long, alias = "max_divergence", default_value_t = 0.20)]
        max_divergence: f64,
        #[arg(long, default_value_t = 100)]
        min_length: usize,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_seqcap_align`: align each locus in a
    /// monolithic UCE FASTA with MAFFT, optionally trim, and write NEXUS
    /// output. (`--aligner muscle` is not yet implemented.)
    SeqcapAlign {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        taxa: usize,
        #[arg(long, default_value_t = false)]
        incomplete_matrix: bool,
        #[arg(long, default_value_t = false)]
        no_trim: bool,
        #[arg(long, default_value_t = false)]
        ambiguous: bool,
        #[arg(long, default_value_t = 20)]
        window: usize,
        #[arg(long, default_value_t = 0.65)]
        proportion: f64,
        #[arg(long, default_value_t = 0.65)]
        threshold: f64,
        #[arg(long, default_value_t = 0.20)]
        max_divergence: f64,
        #[arg(long, default_value_t = 100)]
        min_length: usize,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_get_informative_sites`.
    GetInformativeSites {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long, default_value = "nexus")]
        input_format: String,
    },
    /// Equivalent to `phyluce_align_get_align_summary_data`: print aggregate
    /// statistics and optionally write the per-alignment CSV.
    GetAlignSummaryData {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long)]
        output_stats: Option<PathBuf>,
        #[arg(long, default_value_t = false)]
        show_taxon_counts: bool,
        /// Number of alignment files to summarize concurrently.
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_concatenate_alignments`.
    ConcatenateAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = false)]
        nexus: bool,
        #[arg(long, default_value_t = false)]
        phylip: bool,
    },
    /// Equivalent to `phyluce_align_add_missing_data_designators`. Output
    /// format is always NEXUS (the only format tested against a fixture;
    /// see docs/rust-rewrite-plan.md's phased approach).
    AddMissingDataDesignators {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        match_count_output: PathBuf,
        #[arg(long)]
        incomplete_matrix: Option<PathBuf>,
        #[arg(long, default_value_t = 3)]
        min_taxa: usize,
        #[arg(long, default_value = "?")]
        missing_character: char,
        #[arg(long, default_value_t = false)]
        verbatim: bool,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(
            long = "no-check-missing",
            default_value_t = true,
            action = ArgAction::SetFalse
        )]
        check_missing: bool,
    },
    /// Equivalent to `phyluce_align_remove_empty_taxa`.
    RemoveEmptyTaxa {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value = "nexus")]
        output_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_get_ry_recoded_alignments`. Output
    /// format is always NEXUS.
    GetRyRecodedAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value_t = false)]
        binary: bool,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_extract_taxa_from_alignments`.
    ExtractTaxaFromAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value = "fasta")]
        output_format: String,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        exclude: Vec<String>,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        include: Vec<String>,
    },
    /// Equivalent to `phyluce_align_split_concat_nexus_to_loci`.
    SplitConcatNexusToLoci {
        #[arg(long)]
        nexus: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        output_format: String,
    },
    /// Equivalent to `phyluce_align_filter_alignments`.
    FilterAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        containing_data_for: Vec<String>,
        #[arg(long, default_value_t = 0)]
        min_length: usize,
        #[arg(long, default_value_t = 0)]
        min_taxa: usize,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_convert_one_align_to_another`
    /// (fasta/nexus only; see `convert_cmd` docs for what's not yet
    /// implemented).
    ConvertOneAlignToAnother {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value = "fasta")]
        output_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed`.
    GetGblocksTrimmedAlignmentsFromUntrimmed {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "fasta")]
        input_format: String,
        #[arg(long, default_value_t = 0.5)]
        b1: f64,
        #[arg(long, default_value_t = 0.85)]
        b2: f64,
        #[arg(long, default_value_t = 8)]
        b3: u32,
        #[arg(long, default_value_t = 10)]
        b4: u32,
        #[arg(long, default_value = "nexus")]
        output_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_get_trimal_trimmed_alignments_from_untrimmed`.
    GetTrimalTrimmedAlignmentsFromUntrimmed {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "fasta")]
        input_format: String,
        #[arg(long, default_value = "nexus")]
        output_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_convert_degen_bases`.
    ConvertDegenBases {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value = "nexus")]
        output_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_explode_alignments`.
    ExplodeAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "fasta")]
        input_format: String,
        #[arg(long)]
        conf: Option<PathBuf>,
        #[arg(long)]
        section: Option<String>,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        exclude: Vec<String>,
        #[arg(long, default_value_t = false)]
        by_taxon: bool,
        #[arg(long, default_value_t = false)]
        include_locus: bool,
    },
    /// Equivalent to `phyluce_align_extract_taxon_fasta_from_alignments`.
    ExtractTaxonFastaFromAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        taxon: String,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
    },
    /// Equivalent to `phyluce_align_format_concatenated_phylip_for_paml`.
    FormatConcatenatedPhylipForPaml {
        #[arg(long)]
        phylip_alignment: PathBuf,
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
    /// Equivalent to `phyluce_align_get_incomplete_matrix_estimates`.
    GetIncompleteMatrixEstimates {
        #[arg(long)]
        db: PathBuf,
        #[arg(long, default_value_t = 0.0)]
        min: f64,
        #[arg(long, default_value_t = 1.0)]
        max: f64,
        #[arg(long, default_value_t = 0.1)]
        step: f64,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        exclude: Vec<String>,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        include: Vec<String>,
    },
    /// Equivalent to `phyluce_align_get_only_loci_with_min_taxa`.
    GetOnlyLociWithMinTaxa {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        taxa: usize,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 0.75)]
        percent: f64,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_get_taxon_locus_counts_in_alignments`.
    GetTaxonLocusCountsInAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long, default_value = "fasta")]
        input_format: String,
        #[arg(long)]
        output: PathBuf,
    },
    /// Equivalent to `phyluce_align_move_align_by_conf_file`.
    MoveAlignByConfFile {
        #[arg(long)]
        conf: PathBuf,
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, num_args = 0.., value_delimiter = ' ')]
        sections: Vec<String>,
        #[arg(long, default_value_t = false)]
        opposite: bool,
        #[arg(long, default_value = "nex")]
        extension: String,
    },
    /// Equivalent to `phyluce_align_randomly_sample_and_concatenate`. Uses
    /// a seeded PRNG for sampling instead of `numpy.random.choice` -- see
    /// `randomly_sample_concat_cmd` docs.
    RandomlySampleAndConcatenate {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 100)]
        sample_size: usize,
        #[arg(long, default_value_t = 1)]
        replicates: usize,
    },
    /// Equivalent to `phyluce_align_reduce_alignments_with_raxml`.
    ReduceAlignmentsWithRaxml {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "phylip-relaxed")]
        input_format: String,
    },
    /// Equivalent to `phyluce_align_remove_locus_name_from_files`
    /// (fasta/nexus output only, matching `convert_cmd`'s convention).
    RemoveLocusNameFromFiles {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        taxa: Option<usize>,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value = "nexus")]
        output_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_screen_alignments_for_problems`.
    ScreenAlignmentsForProblems {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = false)]
        do_not_screen_n: bool,
        #[arg(long, default_value_t = false)]
        do_not_screen_x: bool,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long, default_value_t = 1)]
        cores: usize,
    },
    /// Equivalent to `phyluce_align_get_smilogram_from_alignments`. Ties
    /// in major-allele selection are broken deterministically instead of
    /// via `random.choice` -- see `smilogram_cmd` docs.
    GetSmilogramFromAlignments {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output_file: PathBuf,
        #[arg(long)]
        output_missing: PathBuf,
        #[arg(long)]
        output_database: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
    },
}

#[derive(Subcommand)]
enum ExternalAction {
    /// Resolve `[program] binary` via phyluce.conf and run it with `--version`.
    Check {
        #[arg(long)]
        program: String,
        #[arg(long)]
        binary: String,
    },
}

fn main() -> anyhow::Result<()> {
    let raw_args: Vec<OsString> = std::env::args_os().collect();
    let program_name = program_stem(raw_args.first()).unwrap_or_else(|| "phyluce".to_string());
    let cli = Cli::parse_from(expand_legacy_argv(raw_args.clone()));
    init_file_logging(cli.log_path.as_deref(), &program_name, &cli.verbosity)?;
    tracing::info!("{}", format!(" Starting {program_name} ").center(65, '='));
    tracing::info!("Version: {}", env!("CARGO_PKG_VERSION"));
    tracing::info!("Arguments: {}", display_args(&raw_args));
    let result = run_command(cli.command);
    match &result {
        Ok(()) => tracing::info!("{}", format!(" Completed {program_name} ").center(65, '=')),
        Err(err) => tracing::error!("{err:?}"),
    }
    result
}

trait Center {
    fn center(&self, width: usize, fill: char) -> String;
}

impl Center for str {
    fn center(&self, width: usize, fill: char) -> String {
        let len = self.chars().count();
        if len >= width {
            return self.to_string();
        }
        let pad = width - len;
        let left = pad / 2;
        let right = pad - left;
        format!(
            "{}{}{}",
            fill.to_string().repeat(left),
            self,
            fill.to_string().repeat(right)
        )
    }
}

fn run_command(command: Commands) -> anyhow::Result<()> {
    match command {
        Commands::Config { action } => run_config(action),
        Commands::Io { action } => run_io(action),
        Commands::Assembly { action } => run_assembly(action),
        Commands::External { action } => run_external(action),
        Commands::Utilities { action } => run_utilities(action),
        Commands::Align { action } => run_align(action),
        Commands::Ncbi { action } => run_ncbi(action),
        Commands::Genetrees { action } => run_genetrees(action),
        Commands::Workflow {
            config,
            output,
            workflow,
            cores,
            dryrun,
        } => workflow_cmd::run(&config, &output, &workflow, cores, dryrun),
        Commands::Probe { action } => run_probe(action),
    }
}

fn program_stem(arg0: Option<&OsString>) -> Option<String> {
    arg0.and_then(|p| Path::new(p).file_name())
        .and_then(|p| p.to_str())
        .map(|s| s.to_string())
}

fn display_args(args: &[OsString]) -> String {
    args.iter()
        .map(|a| a.to_string_lossy())
        .collect::<Vec<_>>()
        .join(" ")
}

#[derive(Clone)]
struct SharedLogWriter {
    file: Arc<Mutex<File>>,
}

impl<'a> tracing_subscriber::fmt::MakeWriter<'a> for SharedLogWriter {
    type Writer = SharedLogWriterGuard;

    fn make_writer(&'a self) -> Self::Writer {
        SharedLogWriterGuard {
            file: self.file.clone(),
        }
    }
}

struct SharedLogWriterGuard {
    file: Arc<Mutex<File>>,
}

impl std::io::Write for SharedLogWriterGuard {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.file
            .lock()
            .map_err(|_| std::io::Error::other("log file lock was poisoned"))?
            .write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.file
            .lock()
            .map_err(|_| std::io::Error::other("log file lock was poisoned"))?
            .flush()
    }
}

fn init_file_logging(
    log_path: Option<&Path>,
    program_name: &str,
    verbosity: &str,
) -> anyhow::Result<()> {
    let Some(log_path) = log_path else {
        return Ok(());
    };
    std::fs::create_dir_all(log_path)
        .with_context(|| format!("creating log directory {}", log_path.display()))?;
    let logfile = log_path.join(format!("{program_name}.log"));
    let writer = SharedLogWriter {
        file: Arc::new(Mutex::new(
            File::create(&logfile)
                .with_context(|| format!("creating log file {}", logfile.display()))?,
        )),
    };
    let max_level = match verbosity {
        "CRITICAL" => LevelFilter::ERROR,
        "WARN" => LevelFilter::WARN,
        _ => LevelFilter::INFO,
    };
    let subscriber = tracing_subscriber::fmt()
        .with_writer(writer)
        .with_ansi(false)
        .with_target(false)
        .with_max_level(max_level)
        .finish();
    tracing::subscriber::set_global_default(subscriber)
        .map_err(|err| anyhow::anyhow!("initializing tracing subscriber: {err}"))?;
    Ok(())
}

fn expand_legacy_argv(args: Vec<OsString>) -> Vec<OsString> {
    let Some(program) = args
        .first()
        .and_then(|p| Path::new(p).file_name())
        .and_then(|p| p.to_str())
    else {
        return args;
    };
    let Some(prefix) = legacy_command_prefix(program) else {
        return args;
    };

    let mut expanded = Vec::with_capacity(args.len() + prefix.len());
    expanded.push(OsString::from("phyluce"));
    expanded.extend(prefix.iter().map(OsString::from));
    expanded.extend(args.into_iter().skip(1));
    expanded
}

fn legacy_command_prefix(program: &str) -> Option<&'static [&'static str]> {
    LEGACY_COMMANDS
        .iter()
        .find_map(|(name, prefix)| (*name == program).then_some(*prefix))
}

const LEGACY_COMMANDS: &[(&str, &[&str])] = &[
    (
        "phyluce_align_add_missing_data_designators",
        &["align", "add-missing-data-designators"],
    ),
    (
        "phyluce_align_concatenate_alignments",
        &["align", "concatenate-alignments"],
    ),
    (
        "phyluce_align_convert_degen_bases",
        &["align", "convert-degen-bases"],
    ),
    (
        "phyluce_align_convert_one_align_to_another",
        &["align", "convert-one-align-to-another"],
    ),
    (
        "phyluce_align_explode_alignments",
        &["align", "explode-alignments"],
    ),
    (
        "phyluce_align_extract_taxa_from_alignments",
        &["align", "extract-taxa-from-alignments"],
    ),
    (
        "phyluce_align_extract_taxon_fasta_from_alignments",
        &["align", "extract-taxon-fasta-from-alignments"],
    ),
    (
        "phyluce_align_filter_alignments",
        &["align", "filter-alignments"],
    ),
    (
        "phyluce_align_format_concatenated_phylip_for_paml",
        &["align", "format-concatenated-phylip-for-paml"],
    ),
    (
        "phyluce_align_get_align_summary_data",
        &["align", "get-align-summary-data"],
    ),
    (
        "phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed",
        &["align", "get-gblocks-trimmed-alignments-from-untrimmed"],
    ),
    (
        "phyluce_align_get_incomplete_matrix_estimates",
        &["align", "get-incomplete-matrix-estimates"],
    ),
    (
        "phyluce_align_get_informative_sites",
        &["align", "get-informative-sites"],
    ),
    (
        "phyluce_align_get_only_loci_with_min_taxa",
        &["align", "get-only-loci-with-min-taxa"],
    ),
    (
        "phyluce_align_get_ry_recoded_alignments",
        &["align", "get-ry-recoded-alignments"],
    ),
    (
        "phyluce_align_get_smilogram_from_alignments",
        &["align", "get-smilogram-from-alignments"],
    ),
    (
        "phyluce_align_get_taxon_locus_counts_in_alignments",
        &["align", "get-taxon-locus-counts-in-alignments"],
    ),
    (
        "phyluce_align_get_trimal_trimmed_alignments_from_untrimmed",
        &["align", "get-trimal-trimmed-alignments-from-untrimmed"],
    ),
    (
        "phyluce_align_get_trimmed_alignments_from_untrimmed",
        &["align", "get-trimmed-alignments-from-untrimmed"],
    ),
    (
        "phyluce_align_move_align_by_conf_file",
        &["align", "move-align-by-conf-file"],
    ),
    (
        "phyluce_align_randomly_sample_and_concatenate",
        &["align", "randomly-sample-and-concatenate"],
    ),
    (
        "phyluce_align_reduce_alignments_with_raxml",
        &["align", "reduce-alignments-with-raxml"],
    ),
    (
        "phyluce_align_remove_empty_taxa",
        &["align", "remove-empty-taxa"],
    ),
    (
        "phyluce_align_remove_locus_name_from_files",
        &["align", "remove-locus-name-from-files"],
    ),
    (
        "phyluce_align_screen_alignments_for_problems",
        &["align", "screen-alignments-for-problems"],
    ),
    ("phyluce_align_seqcap_align", &["align", "seqcap-align"]),
    (
        "phyluce_align_split_concat_nexus_to_loci",
        &["align", "split-concat-nexus-to-loci"],
    ),
    (
        "phyluce_assembly_assemblo_abyss",
        &["assembly", "assemblo-abyss"],
    ),
    (
        "phyluce_assembly_assemblo_spades",
        &["assembly", "assemblo-spades"],
    ),
    (
        "phyluce_assembly_assemblo_velvet",
        &["assembly", "assemblo-velvet"],
    ),
    (
        "phyluce_assembly_explode_get_fastas_file",
        &["assembly", "explode-get-fastas-file"],
    ),
    (
        "phyluce_assembly_extract_contigs_to_barcodes",
        &["assembly", "extract-contigs-to-barcodes"],
    ),
    (
        "phyluce_assembly_get_bed_from_lastz",
        &["assembly", "get-bed-from-lastz"],
    ),
    (
        "phyluce_assembly_get_fasta_lengths",
        &["assembly", "get-fasta-lengths"],
    ),
    (
        "phyluce_assembly_get_fastas_from_match_counts",
        &["assembly", "get-fastas-from-match-counts"],
    ),
    (
        "phyluce_assembly_get_fastq_lengths",
        &["assembly", "get-fastq-lengths"],
    ),
    (
        "phyluce_assembly_get_match_counts",
        &["assembly", "get-match-counts"],
    ),
    (
        "phyluce_assembly_match_contigs_to_barcodes",
        &["assembly", "match-contigs-to-barcodes"],
    ),
    (
        "phyluce_assembly_match_contigs_to_probes",
        &["assembly", "match-contigs-to-probes"],
    ),
    (
        "phyluce_assembly_screen_probes_for_dupes",
        &["assembly", "screen-probes-for-dupes"],
    ),
    (
        "phyluce_genetrees_generate_multilocus_bootstrap_count",
        &["genetrees", "generate-multilocus-bootstrap-count"],
    ),
    (
        "phyluce_genetrees_get_mean_bootrep_support",
        &["genetrees", "get-mean-bootrep-support"],
    ),
    (
        "phyluce_genetrees_get_tree_counts",
        &["genetrees", "get-tree-counts"],
    ),
    (
        "phyluce_genetrees_rename_tree_leaves",
        &["genetrees", "rename-tree-leaves"],
    ),
    (
        "phyluce_genetrees_sort_multilocus_bootstraps",
        &["genetrees", "sort-multilocus-bootstraps"],
    ),
    (
        "phyluce_ncbi_chunk_fasta_for_ncbi",
        &["ncbi", "chunk-fasta-for-ncbi"],
    ),
    (
        "phyluce_ncbi_prep_uce_align_files_for_ncbi",
        &["ncbi", "prep-uce-align-files-for-ncbi"],
    ),
    ("phyluce_probe_easy_lastz", &["probe", "easy-lastz"]),
    (
        "phyluce_probe_get_genome_sequences_from_bed",
        &["probe", "get-genome-sequences-from-bed"],
    ),
    (
        "phyluce_probe_get_locus_bed_from_lastz_files",
        &["probe", "get-locus-bed-from-lastz-files"],
    ),
    (
        "phyluce_probe_get_multi_fasta_table",
        &["probe", "get-multi-fasta-table"],
    ),
    (
        "phyluce_probe_get_multi_merge_table",
        &["probe", "get-multi-merge-table"],
    ),
    (
        "phyluce_probe_get_probe_bed_from_lastz_files",
        &["probe", "get-probe-bed-from-lastz-files"],
    ),
    (
        "phyluce_probe_get_screened_loci_by_proximity",
        &["probe", "get-screened-loci-by-proximity"],
    ),
    (
        "phyluce_probe_get_subsets_of_tiled_probes",
        &["probe", "get-subsets-of-tiled-probes"],
    ),
    (
        "phyluce_probe_get_tiled_probe_from_multiple_inputs",
        &["probe", "get-tiled-probe-from-multiple-inputs"],
    ),
    (
        "phyluce_probe_get_tiled_probes",
        &["probe", "get-tiled-probes"],
    ),
    (
        "phyluce_probe_query_multi_fasta_table",
        &["probe", "query-multi-fasta-table"],
    ),
    (
        "phyluce_probe_query_multi_merge_table",
        &["probe", "query-multi-merge-table"],
    ),
    (
        "phyluce_probe_reconstruct_uce_from_probe",
        &["probe", "reconstruct-uce-from-probe"],
    ),
    (
        "phyluce_probe_remove_duplicate_hits_from_probes_using_lastz",
        &["probe", "remove-duplicate-hits-from-probes-using-lastz"],
    ),
    (
        "phyluce_probe_remove_overlapping_probes_given_config",
        &["probe", "remove-overlapping-probes-given-config"],
    ),
    (
        "phyluce_probe_run_multiple_lastzs_sqlite",
        &["probe", "run-multiple-lastzs-sqlite"],
    ),
    (
        "phyluce_probe_slice_sequence_from_genomes",
        &["probe", "slice-sequence-from-genomes"],
    ),
    (
        "phyluce_probe_strip_masked_loci_from_set",
        &["probe", "strip-masked-loci-from-set"],
    ),
    (
        "phyluce_utilities_combine_reads",
        &["utilities", "combine-reads"],
    ),
    (
        "phyluce_utilities_filter_bed_by_fasta",
        &["utilities", "filter-bed-by-fasta"],
    ),
    (
        "phyluce_utilities_get_bed_from_fasta",
        &["utilities", "get-bed-from-fasta"],
    ),
    (
        "phyluce_utilities_merge_multiple_gzip_files",
        &["utilities", "merge-multiple-gzip-files"],
    ),
    (
        "phyluce_utilities_merge_next_seq_gzip_files",
        &["utilities", "merge-next-seq-gzip-files"],
    ),
    (
        "phyluce_utilities_replace_many_links",
        &["utilities", "replace-many-links"],
    ),
    (
        "phyluce_utilities_sample_reads_from_files",
        &["utilities", "sample-reads-from-files"],
    ),
    (
        "phyluce_utilities_unmix_fasta_reads",
        &["utilities", "unmix-fasta-reads"],
    ),
    ("phyluce_workflow", &["workflow"]),
];

fn run_probe(action: ProbeAction) -> anyhow::Result<()> {
    match action {
        ProbeAction::RemoveOverlappingProbesGivenConfig {
            probes,
            config,
            output,
        } => remove_overlapping_probes_cmd::run(&probes, &config, &output),
        ProbeAction::GetProbeBedFromLastzFiles { alignments, output } => {
            probe_bed_from_lastz_cmd::run_probe_bed(&alignments, &output)
        }
        ProbeAction::GetLocusBedFromLastzFiles {
            alignments,
            output,
            regex,
        } => probe_bed_from_lastz_cmd::run_locus_bed(&alignments, &output, &regex),
        ProbeAction::GetSubsetsOfTiledProbes {
            probes,
            taxa,
            output,
            regex,
        } => subsets_tiled_probes_cmd::run(&probes, &taxa, &output, &regex),
        ProbeAction::GetMultiFastaTable {
            fastas,
            output,
            base_taxon,
        } => multi_fasta_table_cmd::run_get(&fastas, &output, &base_taxon),
        ProbeAction::QueryMultiFastaTable {
            db,
            base_taxon,
            specific_counts,
            output,
        } => multi_fasta_table_cmd::run_query(&db, &base_taxon, specific_counts, output.as_deref()),
        ProbeAction::GetMultiMergeTable {
            conf,
            output,
            base_taxon,
        } => multi_merge_table_cmd::run_get(&conf, &output, &base_taxon),
        ProbeAction::QueryMultiMergeTable {
            db,
            base_taxon,
            specific_counts,
            output,
        } => multi_merge_table_cmd::run_query(&db, &base_taxon, specific_counts, output.as_deref()),
        ProbeAction::GetScreenedLociByProximity {
            input,
            output,
            distance,
            regex,
        } => screened_loci_proximity_cmd::run(&input, &output, distance, &regex),
        ProbeAction::RemoveDuplicateHitsFromProbesUsingLastz {
            fasta,
            lastz,
            probe_prefix,
            probe_regex,
            probe_bed,
            locus_bed,
            long,
        } => remove_duplicate_hits_cmd::run(
            &fasta,
            &lastz,
            &probe_prefix,
            &probe_regex,
            probe_bed.as_deref(),
            locus_bed.as_deref(),
            long,
        ),
        ProbeAction::GetTiledProbeFromMultipleInputs {
            fastas,
            multi_fasta_output,
            output,
            probe_prefix,
            designer,
            design,
            probe_length,
            tiling_density,
            masking,
            remove_ambiguous,
            remove_gc,
            start_index,
            two_probes,
        } => {
            let args = tiled_probe_from_multiple_inputs_cmd::TilingArgs {
                probe_prefix,
                designer,
                design,
                length: probe_length,
                density: tiling_density,
                mask: masking,
                remove_ambiguous,
                remove_gc,
                start_index,
                two_probes,
            };
            tiled_probe_from_multiple_inputs_cmd::run(&fastas, &multi_fasta_output, &output, &args)
        }
        ProbeAction::GetTiledProbes {
            input,
            output,
            probe_prefix,
            designer,
            design,
            probe_length,
            tiling_density,
            overlap,
            probe_bed,
            locus_bed,
            masking,
            remove_ambiguous,
            remove_gc,
            start_index,
            two_probes,
        } => {
            anyhow::ensure!(
                overlap == "middle" || overlap == "flush-left",
                "--overlap must be 'middle' or 'flush-left'"
            );
            let args = tiled_probes_cmd::TiledProbesArgs {
                probe_prefix,
                designer,
                design,
                length: probe_length,
                density: tiling_density,
                overlap_flush_left: overlap == "flush-left",
                mask: masking,
                remove_ambiguous,
                remove_gc,
                start_index,
                two_probes,
            };
            tiled_probes_cmd::run(
                &input,
                &output,
                probe_bed.as_deref(),
                locus_bed.as_deref(),
                &args,
            )
        }
        ProbeAction::ReconstructUceFromProbe {
            input,
            output,
            muscle_binary,
            mafft_binary,
        } => reconstruct_uce_from_probe_cmd::run(
            &input,
            &output,
            muscle_binary.as_deref(),
            mafft_binary.as_deref(),
        ),
        ProbeAction::GetGenomeSequencesFromBed {
            bed,
            twobit,
            output,
            filter_mask,
            max_n,
            buffer_to,
        } => genome_sequences_from_bed_cmd::run(
            &bed,
            &twobit,
            &output,
            Some(filter_mask),
            max_n,
            buffer_to,
        ),
        ProbeAction::StripMaskedLociFromSet {
            bed,
            twobit,
            output,
            filter_mask,
            max_n,
            min_length,
        } => strip_masked_loci_cmd::run(&bed, &twobit, &output, filter_mask, max_n, min_length),
        ProbeAction::SliceSequenceFromGenomes {
            conf,
            lastz,
            output,
            name_pattern,
            probe_prefix,
            probe_regex,
            exclude,
            contig_orient,
            flank,
            probes,
        } => {
            anyhow::ensure!(
                flank.is_some() != probes.is_some(),
                "exactly one of --flank or --probes must be given"
            );
            let conf_text = std::fs::read_to_string(&conf)
                .with_context(|| format!("reading config {}", conf.display()))?;
            let sections = conf::parse_ini(&conf_text);
            let mut genomes = Vec::new();
            for section in ["chromos", "scaffolds"] {
                if let Some(entries) = sections.get(section) {
                    for (short_name, twobit_path) in entries {
                        let long_name = match &name_pattern {
                            Some(pattern) => pattern.replace("{}", short_name),
                            None => short_name.clone(),
                        };
                        genomes.push(slice_sequence_from_genomes_cmd::GenomeEntry {
                            short_name: short_name.clone(),
                            long_name,
                            twobit_path: PathBuf::from(twobit_path),
                        });
                    }
                }
            }
            let args = slice_sequence_from_genomes_cmd::SliceArgs {
                probe_regex,
                probe_prefix,
                exclude: exclude.into_iter().collect(),
                contig_orient,
                flank,
                probes,
            };
            slice_sequence_from_genomes_cmd::run(&genomes, &lastz, &output, &args)
        }
        ProbeAction::EasyLastz {
            target,
            query,
            output,
            identity,
            coverage,
            min_match,
        } => easy_lastz_cmd::run(&target, &query, &output, coverage, identity, min_match),
        ProbeAction::EasyStampy {
            species,
            assembly,
            genome_files,
            index_prefix,
            reads,
            substitution_rate,
            threads,
            output,
            bam,
            force_rebuild_index,
        } => easy_stampy_cmd::run(
            &species,
            &assembly,
            &genome_files,
            &index_prefix,
            &reads,
            substitution_rate,
            threads,
            &output,
            bam,
            force_rebuild_index,
        ),
        ProbeAction::RunMultipleLastzsSqlite {
            db,
            output,
            probefile,
            chromolist,
            scaffoldlist,
            append,
            no_dir,
            cores,
            genome_base_path,
            coverage,
            identity,
        } => {
            let args = run_multiple_lastzs_sqlite_cmd::RunMultipleLastzsArgs {
                chromolist,
                scaffoldlist,
                append,
                no_dir,
                genome_base_path,
                coverage,
                identity,
                cores,
            };
            run_multiple_lastzs_sqlite_cmd::run(&db, &output, &probefile, &args)
        }
    }
}

fn run_genetrees(action: GenetreesAction) -> anyhow::Result<()> {
    match action {
        GenetreesAction::RenameTreeLeaves {
            input,
            config,
            output,
            section,
            order,
            reroot,
        } => rename_tree_leaves_cmd::run(
            &input,
            &config,
            &output,
            &section,
            &order,
            reroot.as_deref(),
        ),
        GenetreesAction::GetTreeCounts {
            trees,
            locus_support_output,
            root,
            extension,
            exclude,
        } => tree_counts_cmd::run(&trees, &locus_support_output, &root, &extension, &exclude),
        GenetreesAction::GetMeanBootrepSupport { trees, config } => {
            bootrep_support_cmd::run(&trees, &config)
        }
        GenetreesAction::GenerateMultilocusBootstrapCount {
            alignments,
            bootstrap_replicates,
            directory,
            bootstrap_counts,
            bootreps,
        } => bootstrap_count_cmd::run(
            &alignments,
            &bootstrap_replicates,
            &directory,
            &bootstrap_counts,
            bootreps,
        ),
        GenetreesAction::SortMultilocusBootstraps {
            input,
            bootstrap_replicates,
            output,
        } => sort_bootstraps_cmd::run(&input, &bootstrap_replicates, &output),
    }
}

fn run_ncbi(action: NcbiAction) -> anyhow::Result<()> {
    match action {
        NcbiAction::ChunkFastaForNcbi {
            input,
            chunk_size,
            output_prefix,
            output_suffix,
        } => chunk_fasta_cmd::run(&input, chunk_size, &output_prefix, &output_suffix),
        NcbiAction::PrepUceAlignFilesForNcbi {
            alignments,
            conf,
            output,
            input_format,
        } => ncbi_prep_cmd::run(&alignments, &conf, &output, &input_format),
    }
}

fn run_align(action: AlignAction) -> anyhow::Result<()> {
    match action {
        AlignAction::GetTrimmedAlignmentsFromUntrimmed {
            alignments,
            output,
            window,
            proportion,
            threshold,
            max_divergence,
            min_length,
            cores,
        } => get_trimmed_cmd::run(
            &alignments,
            &output,
            window,
            proportion,
            threshold,
            max_divergence,
            min_length,
            cores,
        ),
        AlignAction::SeqcapAlign {
            input,
            output,
            taxa,
            incomplete_matrix,
            no_trim,
            ambiguous,
            window,
            proportion,
            threshold,
            max_divergence,
            min_length,
            cores,
        } => seqcap_align_cmd::run(
            &input,
            &output,
            taxa,
            incomplete_matrix,
            no_trim,
            ambiguous,
            window,
            proportion,
            threshold,
            max_divergence,
            min_length,
            cores,
        ),
        AlignAction::GetInformativeSites {
            alignments,
            output,
            input_format,
        } => informative_sites_cmd::run(&alignments, output, &input_format),
        AlignAction::GetAlignSummaryData {
            alignments,
            input_format,
            output_stats,
            show_taxon_counts,
            cores,
        } => align_summary_cmd::run(
            &alignments,
            &input_format,
            output_stats,
            show_taxon_counts,
            cores,
        ),
        AlignAction::ConcatenateAlignments {
            alignments,
            input_format,
            output,
            nexus,
            phylip,
        } => concatenate_cmd::run(&alignments, &input_format, &output, nexus, phylip),
        AlignAction::AddMissingDataDesignators {
            alignments,
            output,
            match_count_output,
            incomplete_matrix,
            min_taxa,
            missing_character,
            verbatim,
            input_format,
            check_missing,
        } => missing_data_cmd::run(
            &alignments,
            &output,
            &match_count_output,
            incomplete_matrix,
            min_taxa,
            missing_character,
            verbatim,
            &input_format,
            check_missing,
        ),
        AlignAction::RemoveEmptyTaxa {
            alignments,
            output,
            input_format,
            output_format,
            cores,
        } => remove_empty_taxa_cmd::run(&alignments, &output, &input_format, &output_format, cores),
        AlignAction::GetRyRecodedAlignments {
            alignments,
            output,
            input_format,
            binary,
            cores,
        } => ry_recode_cmd::run(&alignments, &output, &input_format, binary, cores),
        AlignAction::ExtractTaxaFromAlignments {
            alignments,
            output,
            input_format,
            output_format,
            exclude,
            include,
        } => extract_taxa_cmd::run(
            &alignments,
            &output,
            &input_format,
            &output_format,
            &exclude,
            &include,
        ),
        AlignAction::SplitConcatNexusToLoci {
            nexus,
            output,
            output_format,
        } => split_concat_cmd::run(&nexus, &output, &output_format),
        AlignAction::FilterAlignments {
            alignments,
            output,
            input_format,
            containing_data_for,
            min_length,
            min_taxa,
            cores,
        } => filter_alignments_cmd::run(
            &alignments,
            &output,
            &input_format,
            &containing_data_for,
            min_length,
            min_taxa,
            cores,
        ),
        AlignAction::ConvertOneAlignToAnother {
            alignments,
            output,
            input_format,
            output_format,
            cores,
        } => convert_cmd::run(&alignments, &output, &input_format, &output_format, cores),
        AlignAction::GetGblocksTrimmedAlignmentsFromUntrimmed {
            alignments,
            output,
            input_format,
            b1,
            b2,
            b3,
            b4,
            output_format,
            cores,
        } => gblocks_cmd::run(
            &alignments,
            &output,
            &input_format,
            b1,
            b2,
            b3,
            b4,
            &output_format,
            cores,
        ),
        AlignAction::GetTrimalTrimmedAlignmentsFromUntrimmed {
            alignments,
            output,
            input_format,
            output_format,
            cores,
        } => trimal_cmd::run(&alignments, &output, &input_format, &output_format, cores),
        AlignAction::ConvertDegenBases {
            alignments,
            output,
            input_format,
            output_format,
            cores,
        } => {
            convert_degen_bases_cmd::run(&alignments, &output, &input_format, &output_format, cores)
        }
        AlignAction::ExplodeAlignments {
            alignments,
            output,
            input_format,
            conf,
            section,
            exclude,
            by_taxon,
            include_locus,
        } => explode_alignments_cmd::run(
            &alignments,
            &output,
            &input_format,
            conf,
            section,
            &exclude,
            by_taxon,
            include_locus,
        ),
        AlignAction::ExtractTaxonFastaFromAlignments {
            alignments,
            taxon,
            output,
            input_format,
        } => extract_taxon_fasta_cmd::run(&alignments, &taxon, &output, &input_format),
        AlignAction::FormatConcatenatedPhylipForPaml {
            phylip_alignment,
            config,
            output,
        } => format_paml_cmd::run(&phylip_alignment, &config, &output),
        AlignAction::GetIncompleteMatrixEstimates {
            db,
            min,
            max,
            step,
            exclude,
            include,
        } => incomplete_matrix_estimates_cmd::run(&db, min, max, step, &exclude, &include),
        AlignAction::GetOnlyLociWithMinTaxa {
            alignments,
            taxa,
            output,
            percent,
            input_format,
            cores,
        } => min_taxa_filter_cmd::run(&alignments, taxa, &output, percent, &input_format, cores),
        AlignAction::GetTaxonLocusCountsInAlignments {
            alignments,
            input_format,
            output,
        } => taxon_locus_counts_cmd::run(&alignments, &input_format, &output),
        AlignAction::MoveAlignByConfFile {
            conf,
            alignments,
            output,
            sections,
            opposite,
            extension,
        } => move_align_cmd::run(&conf, &alignments, &output, &sections, opposite, &extension),
        AlignAction::RandomlySampleAndConcatenate {
            alignments,
            output,
            sample_size,
            replicates,
        } => randomly_sample_concat_cmd::run(&alignments, &output, sample_size, replicates),
        AlignAction::ReduceAlignmentsWithRaxml {
            alignments,
            output,
            input_format,
        } => reduce_raxml_cmd::run(&alignments, &output, &input_format),
        AlignAction::RemoveLocusNameFromFiles {
            alignments,
            output,
            taxa,
            input_format,
            output_format,
            cores,
        } => remove_locus_name_cmd::run(
            &alignments,
            &output,
            taxa,
            &input_format,
            &output_format,
            cores,
        ),
        AlignAction::ScreenAlignmentsForProblems {
            alignments,
            output,
            do_not_screen_n,
            do_not_screen_x,
            input_format,
            cores,
        } => screen_alignments_cmd::run(
            &alignments,
            &output,
            do_not_screen_n,
            do_not_screen_x,
            &input_format,
            cores,
        ),
        AlignAction::GetSmilogramFromAlignments {
            alignments,
            output_file,
            output_missing,
            output_database,
            input_format,
        } => smilogram_cmd::run(
            &alignments,
            &output_file,
            &output_missing,
            &output_database,
            &input_format,
        ),
    }
}

fn run_config(action: ConfigAction) -> anyhow::Result<()> {
    let cfg = PhyluceConfig::load()?;
    match action {
        ConfigAction::Inspect => {
            crate::cli_info!(
                "default config: {}",
                cfg.default_path
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|| "<embedded>".to_string())
            );
            crate::cli_info!(
                "user config:    {}",
                cfg.user_path
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|| "<none>".to_string())
            );
            for section in cfg.section_names() {
                crate::cli_info!("[{section}]");
                if let Some(values) = cfg.get_all_user_params(section) {
                    for v in values {
                        crate::cli_info!("  {v}");
                    }
                }
            }
        }
        ConfigAction::Which { program, binary } => {
            let resolved = cfg.get_user_path(&program, &binary)?;
            crate::cli_info!("{resolved}");
        }
    }
    Ok(())
}

fn run_io(action: IoAction) -> anyhow::Result<()> {
    match action {
        IoAction::ValidateFasta { input } => {
            let issues = validate_fasta(&input)?;
            if issues.is_empty() {
                crate::cli_info!("OK: {} is well-formed FASTA", input.display());
                Ok(())
            } else {
                for issue in &issues {
                    if issue.line > 0 {
                        crate::cli_warn!("{}:{}: {}", input.display(), issue.line, issue.message);
                    } else {
                        crate::cli_warn!("{}: {}", input.display(), issue.message);
                    }
                }
                anyhow::bail!("{} issue(s) found in {}", issues.len(), input.display());
            }
        }
    }
}

fn run_assembly(action: AssemblyAction) -> anyhow::Result<()> {
    match action {
        AssemblyAction::GetFastaLengths { input, csv } => {
            let lengths = fasta_lengths(&input)
                .with_context(|| format!("reading FASTA {}", input.display()))?;
            let report = stats::LengthReport::from_lengths(&lengths);
            if csv {
                crate::cli_info!(
                    "{}",
                    report.to_csv_row(&input.file_name().unwrap_or_default().to_string_lossy())
                );
            } else {
                let message = report.to_human_report();
                print!("{message}");
                tracing::info!(message = %message);
            }
            Ok(())
        }
        AssemblyAction::GetFastqLengths { input, csv } => {
            let mut fastq_files: Vec<PathBuf> = std::fs::read_dir(&input)
                .with_context(|| format!("reading input directory {}", input.display()))?
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| {
                    p.file_name()
                        .and_then(|n| n.to_str())
                        .is_some_and(|n| n.contains(".fastq"))
                })
                .collect();
            // Deterministic order: the legacy `glob.glob` order is
            // filesystem-dependent and effectively undefined; sorting here
            // is a documented, intentional improvement (see
            // docs/rust-rewrite-plan.md section 9's list of quirks worth
            // fixing).
            fastq_files.sort();

            let mut lengths = Vec::new();
            for f in &fastq_files {
                lengths.extend(
                    fastq_lengths(f).with_context(|| format!("reading FASTQ {}", f.display()))?,
                );
            }
            let last_basename = fastq_files
                .last()
                .and_then(|p| p.file_name())
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default();

            let report = stats::FastqLengthReport::from_lengths(&lengths);
            if csv {
                crate::cli_info!("{}", report.to_csv_row(&last_basename));
            } else {
                let message = report.to_human_report();
                print!("{message}");
                tracing::info!(message = %message);
            }
            Ok(())
        }
        AssemblyAction::GetBedFromLastz {
            lastz,
            output,
            identity,
            continuity,
            long_format,
            conf,
            sections,
        } => run_get_bed_from_lastz(
            &lastz,
            &output,
            identity,
            continuity,
            long_format,
            conf,
            sections,
        ),
        AssemblyAction::MatchContigsToProbes {
            contigs,
            probes,
            output,
            min_coverage,
            min_identity,
            dupefile,
            regex,
            keep_duplicates,
            csv,
            skip_alignment,
            force,
        } => match_contigs::run(
            &contigs,
            &probes,
            &output,
            min_coverage,
            min_identity,
            dupefile,
            &regex,
            keep_duplicates,
            csv,
            skip_alignment,
            force,
        ),
        AssemblyAction::GetMatchCounts {
            locus_db,
            taxon_list_config,
            taxon_group,
            output,
            incomplete_matrix,
            optimize,
            random,
            samples,
            sample_size,
            silent,
            keep_counts,
            seed,
            cores,
            extend_locus_db,
        } => match_counts_cmd::run(
            &locus_db,
            &taxon_list_config,
            &taxon_group,
            &output,
            incomplete_matrix,
            extend_locus_db,
            match_counts_cmd::OptimizationOptions {
                optimize,
                random,
                samples,
                sample_size,
                silent,
                keep_counts,
                seed,
                cores,
            },
        ),
        AssemblyAction::GetFastasFromMatchCounts {
            contigs,
            locus_db,
            match_count_output,
            incomplete_matrix,
            output,
            extend_locus_db,
            extend_locus_contigs,
        } => get_fastas_cmd::run(
            &contigs,
            &locus_db,
            &match_count_output,
            incomplete_matrix,
            &output,
            extend_locus_db,
            extend_locus_contigs,
        ),
        AssemblyAction::ExplodeGetFastasFile {
            input,
            output,
            by_taxon,
            split_char,
        } => explode_cmd::run(&input, &output, by_taxon, &split_char),
        AssemblyAction::AssembloSpades {
            output,
            cores,
            memory,
            subfolder,
            no_clean,
            config,
            dir,
        } => assemblo_spades_cmd::run(
            &output,
            cores,
            memory,
            &subfolder,
            no_clean,
            config.as_deref(),
            dir.as_deref(),
        ),
        AssemblyAction::AssembloVelvet {
            output,
            kmer,
            subfolder,
            clean,
            config,
            dir,
        } => assemblo_velvet_cmd::run(
            &output,
            kmer,
            &subfolder,
            clean,
            config.as_deref(),
            dir.as_deref(),
        ),
        AssemblyAction::AssembloAbyss {
            output,
            kmer,
            cores,
            subfolder,
            clean,
            abyss_se,
            config,
            dir,
        } => assemblo_abyss_cmd::run(
            &output,
            kmer,
            cores,
            &subfolder,
            clean,
            abyss_se,
            config.as_deref(),
            dir.as_deref(),
        ),
        AssemblyAction::ScreenProbesForDupes { lastz } => screen_probes_dupes_cmd::run(&lastz),
        AssemblyAction::ExtractContigsToBarcodes {
            contigs,
            config,
            output,
        } => extract_contigs_to_barcodes_cmd::run(&contigs, &config, &output),
        AssemblyAction::MatchContigsToBarcodes {
            contigs,
            barcodes,
            output,
            no_bold,
            database,
        } => match_contigs_to_barcodes_cmd::run(&contigs, &barcodes, &output, no_bold, &database),
    }
}

/// Parse a `[section]` name-per-line conf file into the set of item names
/// across the requested sections, tolerating both bare `name` lines
/// (`allow_no_value`-style) and `name:value`/`name=value` lines.
fn read_conf_section_items(path: &Path, sections: &[String]) -> anyhow::Result<HashSet<String>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading config {}", path.display()))?;
    let wanted: HashSet<&str> = sections.iter().map(|s| s.as_str()).collect();
    Ok(crate::conf::parse_ini(&text)
        .into_iter()
        .filter(|(section, _)| wanted.contains(section.as_str()))
        .flat_map(|(_, entries)| entries.into_iter().map(|(key, _)| key))
        .collect())
}

#[allow(clippy::too_many_arguments)]
fn run_get_bed_from_lastz(
    lastz_path: &Path,
    output: &Path,
    identity: f64,
    continuity: f64,
    long_format: bool,
    conf: Option<PathBuf>,
    sections: Option<Vec<String>>,
) -> anyhow::Result<()> {
    // Mirrors the legacy `if args.conf and args.sections:` gate exactly:
    // supplying --conf without --sections silently disables the filter
    // *and* the fallback "unconditional write" path, so nothing gets
    // written for matches that would otherwise pass identity/continuity.
    let items = match (&conf, &sections) {
        (Some(c), Some(s)) => Some(read_conf_section_items(c, s)?),
        _ => None,
    };

    let matches = read_lastz(lastz_path, long_format)
        .with_context(|| format!("reading lastz results {}", lastz_path.display()))?;
    let mut out = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;
    for m in &matches {
        let name = match m.name2.split('|').nth(1) {
            Some(n) => n.to_string(),
            None => m.name2.split(' ').next().unwrap_or(&m.name2).to_string(),
        };
        if m.percent_identity >= identity && m.percent_continuity >= continuity {
            let should_write = match (&conf, &items) {
                (Some(_), Some(items)) => items.contains(&name),
                (None, _) => true,
                (Some(_), None) => false,
            };
            if should_write {
                writeln!(out, "{}\t{}\t{}\t{}", m.name1, m.zstart1, m.end1, name)?;
            }
        } else {
            crate::cli_info!("{name}");
        }
    }
    Ok(())
}

fn run_utilities(action: UtilitiesAction) -> anyhow::Result<()> {
    match action {
        UtilitiesAction::GetBedFromFasta {
            input,
            output,
            locus_prefix,
        } => {
            let records = read_fasta(&input)
                .with_context(|| format!("reading input FASTA {}", input.display()))?;
            let mut out = std::fs::File::create(&output)
                .with_context(|| format!("creating output file {}", output.display()))?;
            for record in &records {
                let parts: Vec<&str> = record.id.split('|').collect();
                let contig = parts
                    .get(1)
                    .and_then(|s| s.split_once(':'))
                    .map(|(_, v)| v)
                    .ok_or_else(|| {
                        anyhow::anyhow!("record '{}': missing 'contig:' field", record.id)
                    })?;
                let coords = parts
                    .get(2)
                    .and_then(|s| s.split_once(':'))
                    .map(|(_, v)| v)
                    .ok_or_else(|| {
                        anyhow::anyhow!("record '{}': missing 'coords:' field", record.id)
                    })?;
                let locus = parts
                    .get(3)
                    .and_then(|s| s.split_once(':'))
                    .map(|(_, v)| v)
                    .ok_or_else(|| {
                        anyhow::anyhow!("record '{}': missing 'locus:' field", record.id)
                    })?;
                let (begin, end) = coords.split_once('-').ok_or_else(|| {
                    anyhow::anyhow!("record '{}': malformed coords '{}'", record.id, coords)
                })?;
                writeln!(out, "{contig}\t{begin}\t{end}\t{locus_prefix}{locus}")?;
            }
            Ok(())
        }
        UtilitiesAction::FilterBedByFasta { bed, fasta, output } => {
            filter_bed_cmd::run(&bed, &fasta, output)
        }
        UtilitiesAction::ReplaceManyLinks {
            indir,
            oldpath,
            newpath,
            outdir,
        } => replace_links_cmd::run(&indir, &oldpath, &newpath, &outdir),
        UtilitiesAction::CombineReads {
            config,
            output,
            subfolder,
        } => combine_reads_cmd::run(&config, &output, &subfolder),
        UtilitiesAction::MergeMultipleGzipFiles {
            config,
            output,
            section,
            trimmed,
        } => merge_gzip_cmd::run(&config, &output, &section, trimmed),
        UtilitiesAction::MergeNextSeqGzipFiles {
            input,
            config,
            output,
            section,
            se,
        } => merge_nextseq_cmd::run(&input, &config, &output, &section, se),
        UtilitiesAction::UnmixFastaReads {
            mixed_reads,
            singleton_reads,
            out_r1,
            out_r2,
            out_r_singleton,
            new_style,
        } => unmix_fasta_cmd::run(
            &mixed_reads,
            singleton_reads.as_deref(),
            &out_r1,
            &out_r2,
            &out_r_singleton,
            new_style,
        ),
        UtilitiesAction::SampleReadsFromFiles { conf, output } => {
            sample_reads_cmd::run(&conf, &output)
        }
    }
}

fn run_external(action: ExternalAction) -> anyhow::Result<()> {
    match action {
        ExternalAction::Check { program, binary } => {
            let cfg = PhyluceConfig::load()?;
            let resolved = cfg.get_user_path(&program, &binary)?;
            crate::cli_info!("resolved: {resolved}");
            let report = ExternalCommand::new(&resolved).arg("--version").run()?;
            crate::cli_info!("exit code: {:?}", report.exit_code);
            print!("{}", report.stdout);
            if !report.stdout.is_empty() {
                tracing::info!(external_stdout = %report.stdout);
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod cli_compat_tests {
    use super::*;

    #[test]
    fn expands_legacy_script_names() {
        let args = vec![
            OsString::from("/tmp/phyluce_align_convert_degen_bases"),
            OsString::from("--alignments"),
            OsString::from("in"),
        ];
        let expanded = expand_legacy_argv(args);
        assert_eq!(
            expanded,
            vec![
                OsString::from("phyluce"),
                OsString::from("align"),
                OsString::from("convert-degen-bases"),
                OsString::from("--alignments"),
                OsString::from("in"),
            ]
        );
    }

    #[test]
    fn maps_every_legacy_executable() {
        assert_eq!(LEGACY_COMMANDS.len(), 74);
        for (program, prefix) in LEGACY_COMMANDS {
            let mut args = vec!["phyluce"];
            args.extend_from_slice(prefix);
            args.push("--help");
            let error = match Cli::try_parse_from(args) {
                Ok(_) => panic!("legacy mapping target did not display help: {program}"),
                Err(error) => error,
            };
            assert_eq!(
                error.kind(),
                clap::error::ErrorKind::DisplayHelp,
                "{program}"
            );
        }
        assert_eq!(
            legacy_command_prefix("phyluce_assembly_get_fasta_lengths"),
            Some(&["assembly", "get-fasta-lengths"][..])
        );
        assert_eq!(
            legacy_command_prefix("phyluce_probe_easy_lastz"),
            Some(&["probe", "easy-lastz"][..])
        );
        assert_eq!(
            legacy_command_prefix("phyluce_workflow"),
            Some(&["workflow"][..])
        );
    }

    #[test]
    fn leaves_native_cli_args_unchanged() {
        let args = vec![
            OsString::from("phyluce"),
            OsString::from("align"),
            OsString::from("convert-degen-bases"),
        ];
        assert_eq!(expand_legacy_argv(args.clone()), args);
    }

    #[test]
    fn inverse_boolean_flags_change_their_documented_defaults() {
        let base = [
            "phyluce",
            "probe",
            "get-tiled-probes",
            "--input",
            "in.fasta",
            "--output",
            "out.fasta",
            "--probe-prefix",
            "uce-",
            "--designer",
            "tester",
            "--design",
            "test",
        ];
        let cli = Cli::try_parse_from(base).unwrap();
        assert!(matches!(
            cli.command,
            Commands::Probe {
                action: ProbeAction::GetTiledProbes {
                    remove_ambiguous: true,
                    ..
                }
            }
        ));

        let mut with_flag = base.to_vec();
        with_flag.push("--do-not-remove-ambiguous");
        let cli = Cli::try_parse_from(with_flag).unwrap();
        assert!(matches!(
            cli.command,
            Commands::Probe {
                action: ProbeAction::GetTiledProbes {
                    remove_ambiguous: false,
                    ..
                }
            }
        ));
    }

    #[test]
    fn no_check_missing_disables_missing_locus_validation() {
        let cli = Cli::try_parse_from([
            "phyluce",
            "align",
            "add-missing-data-designators",
            "--alignments",
            "in",
            "--output",
            "out",
            "--match-count-output",
            "matches.conf",
            "--no-check-missing",
        ])
        .unwrap();
        assert!(matches!(
            cli.command,
            Commands::Align {
                action: AlignAction::AddMissingDataDesignators {
                    check_missing: false,
                    ..
                }
            }
        ));
    }

    #[test]
    fn accepts_legacy_underscore_option_names() {
        Cli::try_parse_from([
            "phyluce",
            "align",
            "get-trimmed-alignments-from-untrimmed",
            "--alignments",
            "in",
            "--output",
            "out",
            "--max_divergence",
            "0.1",
        ])
        .unwrap();

        Cli::try_parse_from([
            "phyluce",
            "probe",
            "easy-lastz",
            "--target",
            "target.fasta",
            "--query",
            "query.fasta",
            "--output",
            "out.lastz",
            "--min_match",
            "100",
        ])
        .unwrap();
    }

    #[test]
    fn slice_requires_exactly_one_flank_mode() {
        let base = [
            "phyluce",
            "probe",
            "slice-sequence-from-genomes",
            "--conf",
            "conf.ini",
            "--lastz",
            "matches",
            "--output",
            "out",
        ];
        assert!(Cli::try_parse_from(base).is_err());

        let mut both = base.to_vec();
        both.extend_from_slice(&["--flank", "500", "--probes", "probes.fasta"]);
        assert!(Cli::try_parse_from(both).is_err());

        let mut flank = base.to_vec();
        flank.extend_from_slice(&["--flank", "500", "--contig_orient"]);
        Cli::try_parse_from(flank).unwrap();
    }

    #[test]
    fn parallel_core_counts_reach_command_actions() {
        let summary = Cli::try_parse_from([
            "phyluce",
            "align",
            "get-align-summary-data",
            "--alignments",
            "in",
            "--cores",
            "4",
        ])
        .unwrap();
        assert!(matches!(
            summary.command,
            Commands::Align {
                action: AlignAction::GetAlignSummaryData { cores: 4, .. }
            }
        ));

        let lastz = Cli::try_parse_from([
            "phyluce",
            "probe",
            "run-multiple-lastzs-sqlite",
            "--db",
            "matches.sqlite",
            "--output",
            "lastz",
            "--probefile",
            "probes.fasta",
            "--genome-base-path",
            "genomes",
            "--cores",
            "6",
        ])
        .unwrap();
        assert!(matches!(
            lastz.command,
            Commands::Probe {
                action: ProbeAction::RunMultipleLastzsSqlite { cores: 6, .. }
            }
        ));

        let trimming = Cli::try_parse_from([
            "phyluce",
            "align",
            "get-trimmed-alignments-from-untrimmed",
            "--alignments",
            "in",
            "--output",
            "out",
            "--cores",
            "3",
        ])
        .unwrap();
        assert!(matches!(
            trimming.command,
            Commands::Align {
                action: AlignAction::GetTrimmedAlignmentsFromUntrimmed { cores: 3, .. }
            }
        ));

        let conversion = Cli::try_parse_from([
            "phyluce",
            "align",
            "convert-one-align-to-another",
            "--alignments",
            "in",
            "--output",
            "out",
            "--cores",
            "5",
        ])
        .unwrap();
        assert!(matches!(
            conversion.command,
            Commands::Align {
                action: AlignAction::ConvertOneAlignToAnother { cores: 5, .. }
            }
        ));

        let filtering = Cli::try_parse_from([
            "phyluce",
            "align",
            "filter-alignments",
            "--alignments",
            "in",
            "--output",
            "out",
            "--cores",
            "7",
        ])
        .unwrap();
        assert!(matches!(
            filtering.command,
            Commands::Align {
                action: AlignAction::FilterAlignments { cores: 7, .. }
            }
        ));
    }
}
