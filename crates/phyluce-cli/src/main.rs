use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use clap::{Parser, Subcommand};
use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;
use phyluce_io::fastq::fastq_lengths;
use phyluce_io::lastz::read_lastz;
use phyluce_io::{fasta_lengths, read_fasta, validate_fasta};
use tracing::level_filters::LevelFilter;

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
mod probe_bed_from_lastz_cmd;
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
    /// Probe-domain commands (mirrors `bin/probes/phyluce_probe_*`). Only
    /// commands that don't invoke the `lastz` binary directly (`easy-lastz`,
    /// `run-multiple-lastzs-sqlite`) are still unported.
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
        #[arg(long, default_value_t = true)]
        do_not_remove_ambiguous: bool,
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
        #[arg(long, default_value_t = true)]
        do_not_remove_ambiguous: bool,
        #[arg(long, default_value_t = false)]
        remove_gc: bool,
        #[arg(long, default_value_t = 0)]
        start_index: usize,
        #[arg(long, default_value_t = false)]
        two_probes: bool,
    },
    /// Equivalent to `phyluce_probe_reconstruct_uce_from_probe`. Uses MAFFT
    /// instead of MUSCLE for multi-probe loci (MUSCLE isn't available in
    /// this environment) -- see `reconstruct_uce_from_probe_cmd` docs.
    /// `--mafft-binary` is only required if the probe set has multi-probe
    /// loci.
    ReconstructUceFromProbe {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        output: PathBuf,
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
        #[arg(long, default_value_t = false)]
        contig_orient: bool,
        #[arg(long)]
        flank: Option<i64>,
        #[arg(long)]
        probes: Option<i64>,
    },
    /// Equivalent to `phyluce_probe_easy_lastz`. Untested: `lastz` isn't
    /// installed in this environment -- see `lastz_align` docs.
    EasyLastz {
        #[arg(long)]
        target: PathBuf,
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 92.5)]
        identity: f64,
        #[arg(long, default_value_t = 83.0)]
        coverage: f64,
        #[arg(long)]
        min_match: Option<i64>,
    },
    /// Equivalent to `phyluce_probe_run_multiple_lastzs_sqlite`. Untested
    /// (`lastz` not installed) and runs one un-chunked/un-parallelized
    /// `lastz` invocation per genome instead of the Python original's
    /// `multiprocessing`-chunked runner -- see
    /// `run_multiple_lastzs_sqlite_cmd` docs. `--cores` is accepted but
    /// unused.
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
        cores: u32,
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
    /// Equivalent to `phyluce_genetrees_rename_tree_leaves`. `--reroot` is
    /// not yet implemented.
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
        #[arg(long)]
        bootstrap_replicates: PathBuf,
        #[arg(long, default_value = "")]
        directory: String,
        #[arg(long)]
        bootstrap_counts: PathBuf,
        #[arg(long, default_value_t = 100)]
        bootreps: usize,
    },
    /// Equivalent to `phyluce_genetrees_sort_multilocus_bootstraps`.
    SortMultilocusBootstraps {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
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
    /// Equivalent to `phyluce_assembly_get_match_counts` (non-`--optimize`
    /// path): generate a complete- or incomplete-matrix taxon/loci config
    /// from `probe.matches.sqlite`.
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
    /// Equivalent to `phyluce_utilities_merge_multiple_gzip_files`
    /// (`--trimmed` not yet implemented).
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
    /// Equivalent to `phyluce_utilities_sample_reads_from_files`. Untested:
    /// `seqtk` isn't installed in this environment -- see
    /// `sample_reads_cmd` docs.
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
        #[arg(long, default_value_t = 0.20)]
        max_divergence: f64,
        #[arg(long, default_value_t = 100)]
        min_length: usize,
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
    /// Equivalent to `phyluce_align_get_align_summary_data`'s
    /// `--output-stats` CSV (the log-only summary lines aren't reproduced).
    GetAlignSummaryData {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long, default_value = "nexus")]
        input_format: String,
        #[arg(long)]
        output_stats: Option<PathBuf>,
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
        #[arg(long, default_value_t = true)]
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
    },
    /// Equivalent to `phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed`.
    /// Untested against a live Gblocks binary in this environment.
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
    },
    /// Equivalent to `phyluce_align_get_trimal_trimmed_alignments_from_untrimmed`.
    /// Untested against a live trimAl binary in this environment.
    GetTrimalTrimmedAlignmentsFromUntrimmed {
        #[arg(long)]
        alignments: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "fasta")]
        input_format: String,
        #[arg(long, default_value = "nexus")]
        output_format: String,
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
    /// Untested: `raxmlHPC-SSE3` isn't installed in this environment.
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
        self.file.lock().unwrap().write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.file.lock().unwrap().flush()
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
    std::fs::create_dir_all(log_path)?;
    let logfile = log_path.join(format!("{program_name}.log"));
    let writer = SharedLogWriter {
        file: Arc::new(Mutex::new(File::create(logfile)?)),
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
    let _ = tracing::subscriber::set_global_default(subscriber);
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
    match program {
        "phyluce_align_convert_degen_bases" => Some(&["align", "convert-degen-bases"]),
        "phyluce_align_explode_alignments" => Some(&["align", "explode-alignments"]),
        "phyluce_align_extract_taxon_fasta_from_alignments" => {
            Some(&["align", "extract-taxon-fasta-from-alignments"])
        }
        "phyluce_align_format_concatenated_phylip_for_paml" => {
            Some(&["align", "format-concatenated-phylip-for-paml"])
        }
        "phyluce_align_get_incomplete_matrix_estimates" => {
            Some(&["align", "get-incomplete-matrix-estimates"])
        }
        "phyluce_align_get_only_loci_with_min_taxa" => {
            Some(&["align", "get-only-loci-with-min-taxa"])
        }
        "phyluce_align_get_taxon_locus_counts_in_alignments" => {
            Some(&["align", "get-taxon-locus-counts-in-alignments"])
        }
        "phyluce_align_move_align_by_conf_file" => Some(&["align", "move-align-by-conf-file"]),
        "phyluce_align_randomly_sample_and_concatenate" => {
            Some(&["align", "randomly-sample-and-concatenate"])
        }
        "phyluce_align_reduce_alignments_with_raxml" => {
            Some(&["align", "reduce-alignments-with-raxml"])
        }
        "phyluce_align_remove_locus_name_from_files" => {
            Some(&["align", "remove-locus-name-from-files"])
        }
        "phyluce_align_screen_alignments_for_problems" => {
            Some(&["align", "screen-alignments-for-problems"])
        }
        "phyluce_align_get_smilogram_from_alignments" => {
            Some(&["align", "get-smilogram-from-alignments"])
        }
        "phyluce_assembly_screen_probes_for_dupes" => {
            Some(&["assembly", "screen-probes-for-dupes"])
        }
        "phyluce_assembly_extract_contigs_to_barcodes" => {
            Some(&["assembly", "extract-contigs-to-barcodes"])
        }
        "phyluce_assembly_match_contigs_to_barcodes" => {
            Some(&["assembly", "match-contigs-to-barcodes"])
        }
        _ => None,
    }
}

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
            do_not_remove_ambiguous,
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
                remove_ambiguous: do_not_remove_ambiguous,
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
            do_not_remove_ambiguous,
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
                remove_ambiguous: do_not_remove_ambiguous,
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
            mafft_binary,
        } => reconstruct_uce_from_probe_cmd::run(&input, &output, mafft_binary.as_deref()),
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
            let conf_text = std::fs::read_to_string(&conf)?;
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
        ProbeAction::RunMultipleLastzsSqlite {
            db,
            output,
            probefile,
            chromolist,
            scaffoldlist,
            append,
            no_dir,
            cores: _cores,
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
        } => get_trimmed_cmd::run(
            &alignments,
            &output,
            window,
            proportion,
            threshold,
            max_divergence,
            min_length,
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
        } => align_summary_cmd::run(&alignments, &input_format, output_stats),
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
        } => remove_empty_taxa_cmd::run(&alignments, &output, &input_format, &output_format),
        AlignAction::GetRyRecodedAlignments {
            alignments,
            output,
            input_format,
            binary,
        } => ry_recode_cmd::run(&alignments, &output, &input_format, binary),
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
        } => filter_alignments_cmd::run(
            &alignments,
            &output,
            &input_format,
            &containing_data_for,
            min_length,
            min_taxa,
        ),
        AlignAction::ConvertOneAlignToAnother {
            alignments,
            output,
            input_format,
            output_format,
        } => convert_cmd::run(&alignments, &output, &input_format, &output_format),
        AlignAction::GetGblocksTrimmedAlignmentsFromUntrimmed {
            alignments,
            output,
            input_format,
            b1,
            b2,
            b3,
            b4,
            output_format,
        } => gblocks_cmd::run(
            &alignments,
            &output,
            &input_format,
            b1,
            b2,
            b3,
            b4,
            &output_format,
        ),
        AlignAction::GetTrimalTrimmedAlignmentsFromUntrimmed {
            alignments,
            output,
            input_format,
            output_format,
        } => trimal_cmd::run(&alignments, &output, &input_format, &output_format),
        AlignAction::ConvertDegenBases {
            alignments,
            output,
            input_format,
            output_format,
        } => convert_degen_bases_cmd::run(&alignments, &output, &input_format, &output_format),
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
        } => min_taxa_filter_cmd::run(&alignments, taxa, &output, percent, &input_format),
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
        } => remove_locus_name_cmd::run(&alignments, &output, taxa, &input_format, &output_format),
        AlignAction::ScreenAlignmentsForProblems {
            alignments,
            output,
            do_not_screen_n,
            do_not_screen_x,
            input_format,
        } => screen_alignments_cmd::run(
            &alignments,
            &output,
            do_not_screen_n,
            do_not_screen_x,
            &input_format,
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
            println!(
                "default config: {}",
                cfg.default_path
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|| "<none>".to_string())
            );
            println!(
                "user config:    {}",
                cfg.user_path
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|| "<none>".to_string())
            );
            for section in cfg.section_names() {
                println!("[{section}]");
                if let Some(values) = cfg.get_all_user_params(section) {
                    for v in values {
                        println!("  {v}");
                    }
                }
            }
        }
        ConfigAction::Which { program, binary } => {
            let resolved = cfg.get_user_path(&program, &binary)?;
            println!("{resolved}");
        }
    }
    Ok(())
}

fn run_io(action: IoAction) -> anyhow::Result<()> {
    match action {
        IoAction::ValidateFasta { input } => {
            let issues = validate_fasta(&input)?;
            if issues.is_empty() {
                println!("OK: {} is well-formed FASTA", input.display());
                Ok(())
            } else {
                for issue in &issues {
                    if issue.line > 0 {
                        eprintln!("{}:{}: {}", input.display(), issue.line, issue.message);
                    } else {
                        eprintln!("{}: {}", input.display(), issue.message);
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
            let lengths = fasta_lengths(&input)?;
            let report = stats::LengthReport::from_lengths(&lengths);
            if csv {
                println!(
                    "{}",
                    report.to_csv_row(&input.file_name().unwrap_or_default().to_string_lossy())
                );
            } else {
                print!("{}", report.to_human_report());
            }
            Ok(())
        }
        AssemblyAction::GetFastqLengths { input, csv } => {
            let mut fastq_files: Vec<PathBuf> = std::fs::read_dir(&input)?
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
                lengths.extend(fastq_lengths(f)?);
            }
            let last_basename = fastq_files
                .last()
                .and_then(|p| p.file_name())
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default();

            let report = stats::FastqLengthReport::from_lengths(&lengths);
            if csv {
                println!("{}", report.to_csv_row(&last_basename));
            } else {
                print!("{}", report.to_human_report());
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
            extend_locus_db,
        } => match_counts_cmd::run(
            &locus_db,
            &taxon_list_config,
            &taxon_group,
            &output,
            incomplete_matrix,
            extend_locus_db,
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
    let text = std::fs::read_to_string(path)?;
    let wanted: HashSet<&str> = sections.iter().map(|s| s.as_str()).collect();
    let mut items = HashSet::new();
    let mut current: Option<String> = None;
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            current = Some(line[1..line.len() - 1].trim().to_string());
            continue;
        }
        if let Some(section) = &current {
            if wanted.contains(section.as_str()) {
                let key = line
                    .split_once(':')
                    .or_else(|| line.split_once('='))
                    .map(|(k, _)| k.trim())
                    .unwrap_or(line);
                items.insert(key.to_string());
            }
        }
    }
    Ok(items)
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

    let matches = read_lastz(lastz_path, long_format)?;
    let mut out = std::fs::File::create(output)?;
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
            println!("{name}");
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
            let records = read_fasta(&input)?;
            let mut out = std::fs::File::create(&output)?;
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
            println!("resolved: {resolved}");
            let report = ExternalCommand::new(&resolved).arg("--version").run()?;
            println!("exit code: {:?}", report.exit_code);
            print!("{}", report.stdout);
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
    fn leaves_native_cli_args_unchanged() {
        let args = vec![
            OsString::from("phyluce"),
            OsString::from("align"),
            OsString::from("convert-degen-bases"),
        ];
        assert_eq!(expand_legacy_argv(args.clone()), args);
    }
}
