//! CLI wiring for `phyluce align get-informative-sites`, mirroring
//! `phyluce_align_get_informative_sites`.

use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_align::{nexus::parse_nexus, sites::compute_informative_sites, Alignment};
use phyluce_io::read_fasta;

const NEXUS_EXTENSIONS: &[&str] = &[".nexus", ".nex"];
const FASTA_EXTENSIONS: &[&str] = &[".fasta", ".fsa", ".aln", ".fa"];
const PHYLIP_EXTENSIONS: &[&str] = &[".phylip", ".phy"];
const PHYLIP_RELAXED_EXTENSIONS: &[&str] = &[".phylip", ".phy", ".phylip-relaxed"];
const PHYLIP_SEQUENTIAL_EXTENSIONS: &[&str] = &[".phylip", ".phy", ".phylip-sequential"];
const CLUSTAL_EXTENSIONS: &[&str] = &[".clustal", ".clw"];
const EMBOSS_EXTENSIONS: &[&str] = &[".emboss"];
const STOCKHOLM_EXTENSIONS: &[&str] = &[".stockholm"];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum AlignmentFormat {
    Fasta,
    Nexus,
    Phylip,
    PhylipRelaxed,
    PhylipSequential,
    Clustal,
    Emboss,
    Stockholm,
}

impl AlignmentFormat {
    fn parse(value: &str) -> std::io::Result<Self> {
        match value {
            "fasta" => Ok(Self::Fasta),
            "nexus" => Ok(Self::Nexus),
            "phylip" => Ok(Self::Phylip),
            "phylip-relaxed" => Ok(Self::PhylipRelaxed),
            "phylip-sequential" => Ok(Self::PhylipSequential),
            "clustal" => Ok(Self::Clustal),
            "emboss" => Ok(Self::Emboss),
            "stockholm" => Ok(Self::Stockholm),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("unsupported alignment format {value:?}"),
            )),
        }
    }

    fn extensions(self) -> &'static [&'static str] {
        match self {
            Self::Nexus => NEXUS_EXTENSIONS,
            Self::Fasta => FASTA_EXTENSIONS,
            Self::Phylip => PHYLIP_EXTENSIONS,
            Self::PhylipRelaxed => PHYLIP_RELAXED_EXTENSIONS,
            Self::PhylipSequential => PHYLIP_SEQUENTIAL_EXTENSIONS,
            Self::Clustal => CLUSTAL_EXTENSIONS,
            Self::Emboss => EMBOSS_EXTENSIONS,
            Self::Stockholm => STOCKHOLM_EXTENSIONS,
        }
    }
}

pub fn find_alignment_files(dir: &Path, input_format: &str) -> std::io::Result<Vec<PathBuf>> {
    let extensions = AlignmentFormat::parse(input_format)?.extensions();
    let mut files = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if extensions.iter().any(|ext| name.ends_with(ext)) {
                files.push(path);
            }
        }
    }
    files.sort();
    Ok(files)
}

pub fn load_alignment(path: &Path, input_format: &str) -> anyhow::Result<phyluce_align::Alignment> {
    let format = AlignmentFormat::parse(input_format)?;
    match format {
        AlignmentFormat::Nexus => Ok(parse_nexus(
            &std::fs::read_to_string(path)
                .with_context(|| format!("reading alignment {}", path.display()))?,
        )?),
        AlignmentFormat::Fasta => {
            let records = read_fasta(path)
                .with_context(|| format!("reading alignment {}", path.display()))?;
            let alignment =
                Alignment::from_pairs(records.into_iter().map(|r| (r.id, r.sequence)).collect());
            alignment.validate()?;
            Ok(alignment)
        }
        AlignmentFormat::Phylip
        | AlignmentFormat::PhylipRelaxed
        | AlignmentFormat::PhylipSequential => parse_phylip(
            &std::fs::read_to_string(path)
                .with_context(|| format!("reading alignment {}", path.display()))?,
            format,
        ),
        AlignmentFormat::Clustal => parse_clustal(
            &std::fs::read_to_string(path)
                .with_context(|| format!("reading alignment {}", path.display()))?,
        ),
        AlignmentFormat::Emboss => parse_emboss(
            &std::fs::read_to_string(path)
                .with_context(|| format!("reading alignment {}", path.display()))?,
        ),
        AlignmentFormat::Stockholm => parse_stockholm(
            &std::fs::read_to_string(path)
                .with_context(|| format!("reading alignment {}", path.display()))?,
        ),
    }
}

fn alignment_from_chunks(
    order: Vec<String>,
    chunks: std::collections::HashMap<String, String>,
) -> anyhow::Result<Alignment> {
    let alignment = Alignment::from_pairs(
        order
            .into_iter()
            .map(|id| {
                let sequence = chunks.get(&id).cloned().unwrap_or_default();
                (id, sequence)
            })
            .collect(),
    );
    alignment.validate()?;
    Ok(alignment)
}

fn parse_phylip(text: &str, format: AlignmentFormat) -> anyhow::Result<Alignment> {
    let mut lines = text.lines();
    let header = lines
        .find(|line| !line.trim().is_empty())
        .ok_or_else(|| anyhow::anyhow!("empty PHYLIP file"))?;
    let dimensions: Vec<&str> = header.split_whitespace().collect();
    anyhow::ensure!(dimensions.len() >= 2, "invalid PHYLIP dimensions");
    let ntax: usize = dimensions[0].parse()?;
    let nchar: usize = dimensions[1].parse()?;
    anyhow::ensure!(ntax > 0, "PHYLIP ntax must be greater than zero");

    let relaxed = format != AlignmentFormat::Phylip;
    let mut order = Vec::with_capacity(ntax);
    let mut chunks: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut continuation_index = 0usize;

    for raw in lines.filter(|line| !line.trim().is_empty()) {
        let candidate = if relaxed {
            let mut fields = raw.split_whitespace();
            fields.next().and_then(|label| {
                let sequence: String = fields.collect();
                (!sequence.is_empty()).then(|| (label.to_string(), sequence))
            })
        } else if raw.len() >= 10 {
            let (label, sequence) = raw.split_at(10);
            let label = label.trim();
            let sequence: String = sequence.split_whitespace().collect();
            (!label.is_empty() && !sequence.is_empty()).then(|| (label.to_string(), sequence))
        } else {
            None
        };
        // Interleaved PHYLIP continuation rows contain only sequence chunks.
        // Biopython indents these rows, and their first chunk must not be
        // mistaken for a new taxon label. Repeated explicit labels remain
        // valid in interleaved files produced by other writers.
        let explicit = candidate.and_then(|(label, sequence)| {
            (chunks.contains_key(&label)
                || (order.len() < ntax && !raw.starts_with(char::is_whitespace)))
            .then_some((label, sequence))
        });

        if let Some((label, sequence)) = explicit {
            if let Some(existing) = chunks.get_mut(&label) {
                existing.push_str(&sequence);
                continuation_index =
                    (order.iter().position(|id| id == &label).unwrap_or(0) + 1) % ntax;
            } else if order.len() < ntax {
                order.push(label.clone());
                chunks.insert(label, sequence);
                continuation_index = order.len() % ntax;
            } else {
                anyhow::bail!("unexpected PHYLIP taxon {label:?}");
            }
        } else {
            anyhow::ensure!(!order.is_empty(), "PHYLIP continuation before first taxon");
            let sequence: String = raw.split_whitespace().collect();
            let index = if order.len() < ntax {
                order
                    .iter()
                    .position(|id| chunks[id].len() < nchar)
                    .unwrap_or(order.len() - 1)
            } else {
                continuation_index % ntax
            };
            chunks.get_mut(&order[index]).unwrap().push_str(&sequence);
            continuation_index = (index + 1) % ntax;
        }
    }

    anyhow::ensure!(
        order.len() == ntax,
        "PHYLIP expected {ntax} taxa, found {}",
        order.len()
    );
    for id in &order {
        anyhow::ensure!(
            chunks[id].len() == nchar,
            "PHYLIP taxon {id:?} has {} characters; expected {nchar}",
            chunks[id].len()
        );
    }
    alignment_from_chunks(order, chunks)
}

fn parse_clustal(text: &str) -> anyhow::Result<Alignment> {
    let mut order = Vec::new();
    let mut chunks = std::collections::HashMap::new();
    for raw in text.lines().skip(1) {
        if raw.trim().is_empty() || raw.starts_with(char::is_whitespace) {
            continue;
        }
        let fields: Vec<&str> = raw.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let id = fields[0].to_string();
        if !chunks.contains_key(&id) {
            order.push(id.clone());
        }
        chunks
            .entry(id)
            .or_insert_with(String::new)
            .push_str(fields[1]);
    }
    anyhow::ensure!(!order.is_empty(), "CLUSTAL file contains no sequences");
    alignment_from_chunks(order, chunks)
}

fn parse_stockholm(text: &str) -> anyhow::Result<Alignment> {
    let mut order = Vec::new();
    let mut chunks = std::collections::HashMap::new();
    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line == "//" {
            break;
        }
        let mut fields = line.split_whitespace();
        let id = fields.next().unwrap_or_default().to_string();
        let sequence = fields.next().unwrap_or_default();
        anyhow::ensure!(
            !id.is_empty() && !sequence.is_empty(),
            "invalid Stockholm row"
        );
        if !chunks.contains_key(&id) {
            order.push(id.clone());
        }
        chunks
            .entry(id)
            .or_insert_with(String::new)
            .push_str(sequence);
    }
    anyhow::ensure!(!order.is_empty(), "Stockholm file contains no sequences");
    alignment_from_chunks(order, chunks)
}

fn parse_emboss(text: &str) -> anyhow::Result<Alignment> {
    let mut order = Vec::new();
    let mut chunks = std::collections::HashMap::new();
    for raw in text.lines() {
        if raw.starts_with('#') || raw.trim().is_empty() || raw.starts_with(char::is_whitespace) {
            continue;
        }
        let fields: Vec<&str> = raw.split_whitespace().collect();
        if fields.len() < 4
            || fields[1].parse::<usize>().is_err()
            || fields
                .last()
                .and_then(|v| v.parse::<usize>().ok())
                .is_none()
        {
            continue;
        }
        let id = fields[0].to_string();
        let sequence = fields[2];
        if !chunks.contains_key(&id) {
            order.push(id.clone());
        }
        chunks
            .entry(id)
            .or_insert_with(String::new)
            .push_str(sequence);
    }
    anyhow::ensure!(!order.is_empty(), "EMBOSS file contains no sequences");
    alignment_from_chunks(order, chunks)
}

pub fn run(
    alignments_dir: &Path,
    output: Option<PathBuf>,
    input_format: &str,
) -> anyhow::Result<()> {
    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;
    let mut rows: Vec<(String, usize, usize, usize, usize)> = Vec::new();
    for file in &files {
        let name = file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let alignment = load_alignment(file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        let length = alignment.nchar();
        let (informative, differences, counted) = compute_informative_sites(&alignment);
        rows.push((name, length, informative, differences, counted));
    }

    if let Some(out_path) = output {
        let mut out = String::from("locus,length,informative_sites,differences,counted-bases\n");
        for (name, length, informative, differences, counted) in &rows {
            out.push_str(&format!(
                "{name},{length},{informative},{differences},{counted}\n"
            ));
        }
        std::fs::write(&out_path, out)
            .with_context(|| format!("writing output {}", out_path.display()))?;
    } else {
        crate::cli_info!("locus\tlength\tinformative_sites\tdifferences\tcounted-bases");
        for (name, length, informative, differences, counted) in &rows {
            crate::cli_info!("{name}\t{length}\t{informative}\t{differences}\t{counted}");
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_phylip_relaxed_extensions() {
        let dir =
            std::env::temp_dir().join(format!("phyluce-find-alignments-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("a.phy"), "2 4\nA ACGT\nB ACGT\n").unwrap();
        std::fs::write(dir.join("b.phylip"), "2 4\nA ACGT\nB ACGT\n").unwrap();
        std::fs::write(dir.join("c.fasta"), ">A\nACGT\n").unwrap();

        let names: Vec<String> = find_alignment_files(&dir, "phylip-relaxed")
            .unwrap()
            .into_iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
            .collect();
        assert_eq!(names, vec!["a.phy", "b.phylip"]);

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn rejects_unknown_alignment_formats() {
        assert!(find_alignment_files(Path::new("."), "nexsu").is_err());
    }

    #[test]
    fn parses_supported_text_alignment_formats() {
        let relaxed = parse_phylip(
            "2 4\ntaxon_a ACGT\ntaxon_b ACGA\n",
            AlignmentFormat::PhylipRelaxed,
        )
        .unwrap();
        assert_eq!(relaxed.rows[1].seq, b"ACGA");

        let clustal = parse_clustal("CLUSTAL W\n\ntaxon_a ACGT\ntaxon_b ACGA\n").unwrap();
        assert_eq!(clustal.rows[0].seq, b"ACGT");

        let stockholm =
            parse_stockholm("# STOCKHOLM 1.0\ntaxon_a ACGT\ntaxon_b ACGA\n//\n").unwrap();
        assert_eq!(stockholm.rows.len(), 2);

        let emboss = parse_emboss(
            "#=======================================\ntaxon_a 1 ACGT 4\ntaxon_b 1 ACGA 4\n",
        )
        .unwrap();
        assert_eq!(emboss.rows[1].seq, b"ACGA");
    }

    #[test]
    fn parses_biopython_interleaved_phylip_relaxed() {
        let text = " 3 12\ntaxon_a  ACGT ACGT\ntaxon_b  ACGA ACGA\ntaxon_c  ACGG ACGG\n\n         TTTT\n         GGGG\n         CCCC\n";
        let alignment = parse_phylip(text, AlignmentFormat::PhylipRelaxed).unwrap();
        assert_eq!(alignment.ntax(), 3);
        assert_eq!(alignment.rows[0].seq, b"ACGTACGTTTTT");
        assert_eq!(alignment.rows[1].seq, b"ACGAACGAGGGG");
        assert_eq!(alignment.rows[2].seq, b"ACGGACGGCCCC");
    }
}
