//! CLI wiring for `phyluce align get-informative-sites`, mirroring
//! `phyluce_align_get_informative_sites`.

use std::path::{Path, PathBuf};

use phyluce_align::{nexus::parse_nexus, sites::compute_informative_sites};
use phyluce_io::read_fasta;

const NEXUS_EXTENSIONS: &[&str] = &[".nexus", ".nex"];
const FASTA_EXTENSIONS: &[&str] = &[".fasta", ".fsa", ".aln", ".fa"];
const PHYLIP_EXTENSIONS: &[&str] = &[".phylip", ".phy"];
const PHYLIP_RELAXED_EXTENSIONS: &[&str] = &[".phylip", ".phy", ".phylip-relaxed"];
const PHYLIP_SEQUENTIAL_EXTENSIONS: &[&str] = &[".phylip", ".phy", ".phylip-sequential"];
const CLUSTAL_EXTENSIONS: &[&str] = &[".clustal", ".clw"];
const EMBOSS_EXTENSIONS: &[&str] = &[".emboss"];
const STOCKHOLM_EXTENSIONS: &[&str] = &[".stockholm"];

pub fn find_alignment_files(dir: &Path, input_format: &str) -> std::io::Result<Vec<PathBuf>> {
    let extensions: &[&str] = match input_format {
        "nexus" => NEXUS_EXTENSIONS,
        "fasta" => FASTA_EXTENSIONS,
        "phylip" => PHYLIP_EXTENSIONS,
        "phylip-relaxed" => PHYLIP_RELAXED_EXTENSIONS,
        "phylip-sequential" => PHYLIP_SEQUENTIAL_EXTENSIONS,
        "clustal" => CLUSTAL_EXTENSIONS,
        "emboss" => EMBOSS_EXTENSIONS,
        "stockholm" => STOCKHOLM_EXTENSIONS,
        _ => FASTA_EXTENSIONS,
    };
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
    match input_format {
        "nexus" => Ok(parse_nexus(&std::fs::read_to_string(path)?)?),
        "fasta" => {
            let records = read_fasta(path)?;
            Ok(phyluce_align::Alignment::from_pairs(
                records.into_iter().map(|r| (r.id, r.sequence)).collect(),
            ))
        }
        _ => anyhow::bail!(
            "alignment parsing for input format '{input_format}' is not supported by this command"
        ),
    }
}

pub fn run(
    alignments_dir: &Path,
    output: Option<PathBuf>,
    input_format: &str,
) -> anyhow::Result<()> {
    let files = find_alignment_files(alignments_dir, input_format)?;
    let mut rows: Vec<(String, usize, usize, usize, usize)> = Vec::new();
    for file in &files {
        let name = file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let alignment = load_alignment(file, input_format)?;
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
        std::fs::write(out_path, out)?;
    } else {
        println!("locus\tlength\tinformative_sites\tdifferences\tcounted-bases");
        for (name, length, informative, differences, counted) in &rows {
            println!("{name}\t{length}\t{informative}\t{differences}\t{counted}");
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
}
