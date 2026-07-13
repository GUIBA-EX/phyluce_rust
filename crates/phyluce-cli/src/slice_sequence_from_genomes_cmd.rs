//! CLI wiring for `phyluce probe slice-sequence-from-genomes`, mirroring
//! `phyluce_probe_slice_sequence_from_genomes` (aka
//! `slice_sequence_from_genomes2.py`).

use std::collections::{HashMap, HashSet};
use std::io::Write as _;
use std::path::{Path, PathBuf};

use phyluce_io::lastz::read_lastz;
use phyluce_io::twobit::TwoBitFile;
use regex::Regex;

pub struct GenomeEntry {
    pub short_name: String,
    pub long_name: String,
    pub twobit_path: PathBuf,
}

pub struct SliceArgs {
    pub probe_regex: String,
    pub probe_prefix: String,
    pub exclude: HashSet<String>,
    pub contig_orient: bool,
    pub flank: Option<i64>,
    pub probes: Option<i64>,
}

fn probe_name(header: &str, regex: &Regex) -> anyhow::Result<String> {
    regex
        .captures(header)
        .and_then(|c| c.get(1))
        .map(|m| m.as_str().to_string())
        .ok_or_else(|| anyhow::anyhow!("no regex match for {:?}", header))
}

fn slice_and_return(
    tb: &TwoBitFile,
    name: &str,
    min: i64,
    max: i64,
    flank: Option<i64>,
    probes: Option<i64>,
) -> anyhow::Result<(i64, i64, Vec<u8>)> {
    let seq_len = tb.sequence_len(name)? as i64;
    let (ss, se) = if let Some(probes) = probes {
        let length = (max - min).abs();
        let mut delta = probes - length;
        if delta > 0 {
            if delta % 2 != 0 {
                delta += 1;
            }
            let mut ss = min - delta / 2;
            let mut se = max + delta / 2;
            if ss < 0 {
                ss = 0;
                se += probes - se;
            }
            (ss, se)
        } else {
            (min, max)
        }
    } else {
        let flank = flank.unwrap_or(0);
        let ss = if min - flank > 0 { min - flank } else { 0 };
        let se = if max + flank < seq_len {
            max + flank
        } else {
            seq_len
        };
        (ss, se)
    };
    let seq = tb.read_slice(name, ss, se)?;
    Ok((ss, se, seq))
}

fn strip_matching(seq: &[u8], predicate: impl Fn(u8) -> bool) -> (usize, usize) {
    let mut start = 0usize;
    while start < seq.len() && predicate(seq[start]) {
        start += 1;
    }
    let mut end = seq.len();
    while end > start && predicate(seq[end - 1]) {
        end -= 1;
    }
    (start, end)
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn build_record(
    cnt: usize,
    contig: &str,
    ss: i64,
    se: i64,
    uce: &str,
    min: i64,
    max: i64,
    orient: &str,
    n_positions: usize,
    sequence: Vec<u8>,
    probes: Option<i64>,
) -> (String, i64, i64, Vec<u8>, String) {
    let (ss, se, sequence, name_start, orient) = if probes.is_none() {
        let (n_start, n_end) = strip_matching(&sequence, |b| b == b'N' || b == b'n');
        let (r_start, r_end) = strip_matching(&sequence[n_start..n_end], |b| {
            b == b'a' || b == b'c' || b == b'g' || b == b't'
        });
        let trimmed = sequence[n_start..n_end][r_start..r_end].to_vec();
        let new_ss = ss + n_start as i64 + r_start as i64;
        let new_se = se - (sequence.len() - n_end) as i64 - (n_end - n_start - r_end) as i64;
        let name_start = format!("Node_{cnt}_length_{}_cov_1000", trimmed.len());
        (new_ss, new_se, trimmed, name_start, orient.to_string())
    } else {
        let orient = if orient != "+" {
            "revcomp".to_string()
        } else {
            orient.to_string()
        };
        (ss, se, sequence, format!("slice_{cnt}"), orient)
    };
    let name = format!(
        "{name_start}|contig:{contig}|slice:{ss}-{se}|uce:{uce}|match:{min}-{max}|orient:{orient}|probes:{n_positions}"
    );
    let final_seq = if orient == "revcomp" {
        reverse_complement(&sequence)
    } else {
        sequence
    };
    (name, ss, se, final_seq, orient)
}

pub fn run(
    genomes: &[GenomeEntry],
    lastz_dir: &Path,
    output_dir: &Path,
    args: &SliceArgs,
) -> anyhow::Result<()> {
    std::fs::create_dir_all(output_dir)?;
    let regex_str = args.probe_regex.replace("{}", &args.probe_prefix);
    let regex = Regex::new(&regex_str)?;

    for genome in genomes {
        if args.exclude.contains(&genome.short_name) {
            continue;
        }
        let out_name = crate::output_path::output_file(
            output_dir,
            &format!("{}.fasta", genome.short_name.to_lowercase()),
        )?;
        let mut outf = std::fs::File::create(&out_name)?;

        let tb = TwoBitFile::open(&genome.twobit_path)?;
        let lz_path = lastz_dir.join(&genome.long_name);
        let matches_list = read_lastz(&lz_path, true)?;

        let mut all_uce_names: HashSet<String> = HashSet::new();
        // uce_name -> contig_name -> Vec<(start, end)>
        let mut uce_matches: HashMap<String, HashMap<String, Vec<(i64, i64)>>> = HashMap::new();
        // uce_name -> contig_name -> set of strand chars
        let mut orientation: HashMap<String, HashMap<String, HashSet<String>>> = HashMap::new();

        for m in &matches_list {
            let contig_name = m.name1.clone();
            let uce_name = probe_name(&m.name2, &regex)?;
            all_uce_names.insert(uce_name.clone());
            uce_matches
                .entry(uce_name.clone())
                .or_default()
                .entry(contig_name.clone())
                .or_default()
                .push((m.zstart1, m.end1));
            let strand = if args.contig_orient {
                m.strand1.clone()
            } else {
                m.strand2.clone()
            };
            orientation
                .entry(uce_name)
                .or_default()
                .entry(contig_name)
                .or_default()
                .insert(strand);
        }

        let dupes: HashSet<String> = uce_matches
            .iter()
            .filter(|(_, contigs)| contigs.len() > 1)
            .map(|(uce, _)| uce.clone())
            .collect();
        for d in &dupes {
            uce_matches.remove(d);
        }

        let mut node_count = 0usize;
        let mut orient_drop: HashSet<String> = HashSet::new();
        let mut length_drop: HashSet<String> = HashSet::new();

        let mut uce_names: Vec<&String> = uce_matches.keys().collect();
        uce_names.sort();
        for uce_name in uce_names {
            let matches = &uce_matches[uce_name];
            anyhow::ensure!(matches.len() == 1, "There are multiple UCE matches");
            let (contig_name, positions) = matches.iter().next().unwrap();
            let mut bad = false;

            let orient_set = &orientation[uce_name][contig_name];
            if orient_set.len() > 1 {
                bad = true;
                orient_drop.insert(uce_name.clone());
            }

            let mut sorted_positions = positions.clone();
            sorted_positions.sort();
            if !bad && sorted_positions.len() > 1 {
                for i in 1..sorted_positions.len() {
                    if sorted_positions[i].0 - sorted_positions[i - 1].1 > 500 {
                        bad = true;
                        length_drop.insert(uce_name.clone());
                        break;
                    }
                }
            }

            if bad {
                continue;
            }

            let min = sorted_positions[0].0;
            let max = sorted_positions[sorted_positions.len() - 1].1;
            let (ss, se, sequence) =
                slice_and_return(&tb, contig_name, min, max, args.flank, args.probes)?;
            let orient_char = orient_set.iter().next().unwrap().clone();
            let (name, _ss, _se, final_seq, _orient) = build_record(
                node_count,
                contig_name,
                ss,
                se,
                uce_name,
                min,
                max,
                &orient_char,
                sorted_positions.len(),
                sequence,
                args.probes,
            );
            if let Some(probes) = args.probes {
                if (final_seq.len() as i64) < probes {
                    continue;
                }
            }
            writeln!(outf, ">{name}\n{}", String::from_utf8_lossy(&final_seq))?;
            node_count += 1;
        }

        crate::cli_warn!(
            "{}: {} uces, {} dupes, {} non-dupes, {} orient drop, {} length drop, {node_count} written",
            genome.short_name,
            all_uce_names.len(),
            dupes.len(),
            uce_matches.len(),
            orient_drop.len(),
            length_drop.len(),
        );
    }
    Ok(())
}
