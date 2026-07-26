//! Build the longest coordinate-supported probe reference per
//! locus/source/chromosome group.
//!
//! This is a Rust-only addition intended for tools such as GeneMiner2 that
//! consume one FASTA file per locus. Unlike `reconstruct-uce-from-probe`, it
//! never aligns or collapses references from different source taxa.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_io::FastaRecord;

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct GroupKey {
    locus: String,
    source: String,
    chromosome: String,
    design: String,
    designer: String,
}

#[derive(Clone, Debug)]
struct Fragment {
    id: String,
    local_start: i64,
    local_end: i64,
    global_start: i64,
    sequence: Vec<u8>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct MergedProbe {
    local_start: i64,
    local_end: i64,
    global_start: i64,
    global_end: i64,
    sequence: Vec<u8>,
}

fn metadata(header: &str) -> HashMap<&str, &str> {
    header
        .split_once('|')
        .map(|(_, values)| values)
        .unwrap_or("")
        .split(',')
        .filter_map(|field| field.trim().split_once(':'))
        .collect()
}

fn fallback_locus(id: &str) -> &str {
    match id.rsplit_once("_p") {
        Some((locus, probe)) if !locus.is_empty() && probe.chars().all(|c| c.is_ascii_digit()) => {
            locus
        }
        _ => id,
    }
}

fn required_i64(
    values: &HashMap<&str, &str>,
    key: &str,
    record: &FastaRecord,
) -> anyhow::Result<i64> {
    values
        .get(key)
        .ok_or_else(|| anyhow::anyhow!("record {:?} is missing {key}", record.id))?
        .parse()
        .with_context(|| format!("record {:?} has invalid {key}", record.id))
}

fn parse_fragment(record: FastaRecord) -> anyhow::Result<(GroupKey, Fragment)> {
    let values = metadata(&record.description);
    let locus = values
        .get("probes-locus")
        .copied()
        .unwrap_or_else(|| fallback_locus(&record.id))
        .to_string();
    let chromosome = values
        .get("probes-global-chromo")
        .copied()
        .ok_or_else(|| anyhow::anyhow!("record {:?} is missing probes-global-chromo", record.id))?
        .to_string();
    let local_start = required_i64(&values, "probes-local-start", &record)?;
    let local_end = required_i64(&values, "probes-local-end", &record)?;
    let global_start = required_i64(&values, "probes-global-start", &record)?;
    let global_end = required_i64(&values, "probes-global-end", &record)?;
    let local_span = local_end
        .checked_sub(local_start)
        .filter(|span| *span > 0)
        .ok_or_else(|| anyhow::anyhow!("record {:?} has a non-positive local span", record.id))?;
    let global_span = global_end
        .checked_sub(global_start)
        .filter(|span| *span > 0)
        .ok_or_else(|| anyhow::anyhow!("record {:?} has a non-positive global span", record.id))?;
    anyhow::ensure!(
        global_span == local_span,
        "record {:?} global span {}-{} does not match local span {}-{}",
        record.id,
        global_start,
        global_end,
        local_start,
        local_end
    );
    anyhow::ensure!(
        record.sequence.len() == local_span as usize,
        "record {:?} sequence length {} does not match local span {}-{}",
        record.id,
        record.sequence.len(),
        local_start,
        local_end
    );
    let sequence: Vec<u8> = record
        .sequence
        .bytes()
        .map(|base| base.to_ascii_uppercase())
        .collect();
    if let Some((offset, base)) = sequence
        .iter()
        .copied()
        .enumerate()
        .find(|(_, base)| !matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        anyhow::bail!(
            "record {:?} contains unsupported base {:?} at sequence position {}; full probes require A/C/G/T",
            record.id,
            char::from(base),
            offset
        );
    }

    let key = GroupKey {
        locus,
        source: values
            .get("probes-source")
            .copied()
            .unwrap_or("unknown")
            .to_string(),
        chromosome,
        design: values
            .get("design")
            .copied()
            .unwrap_or("unknown")
            .to_string(),
        designer: values
            .get("designer")
            .copied()
            .unwrap_or("unknown")
            .to_string(),
    };
    let fragment = Fragment {
        id: record.id,
        local_start,
        local_end,
        global_start,
        sequence,
    };
    Ok((key, fragment))
}

fn merge_component(fragments: &[Fragment], key: &GroupKey) -> anyhow::Result<MergedProbe> {
    let local_start = fragments
        .iter()
        .map(|fragment| fragment.local_start)
        .min()
        .unwrap_or(0);
    let local_end = fragments
        .iter()
        .map(|fragment| fragment.local_end)
        .max()
        .unwrap_or(0);
    let sequence_length = local_end
        .checked_sub(local_start)
        .and_then(|length| usize::try_from(length).ok())
        .ok_or_else(|| anyhow::anyhow!("invalid coordinate span for locus {:?}", key.locus))?;
    let mut sequence = vec![None; sequence_length];
    let coordinate_offset = fragments
        .first()
        .and_then(|fragment| fragment.global_start.checked_sub(fragment.local_start))
        .ok_or_else(|| anyhow::anyhow!("invalid coordinate offset for locus {:?}", key.locus))?;

    for fragment in fragments {
        let fragment_offset = fragment
            .global_start
            .checked_sub(fragment.local_start)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "invalid coordinate offset for locus {:?}, record {:?}",
                    key.locus,
                    fragment.id
                )
            })?;
        anyhow::ensure!(
            fragment_offset == coordinate_offset,
            "inconsistent local/global coordinates for locus {:?}, source {:?}, chromosome {:?} (record {:?})",
            key.locus,
            key.source,
            key.chromosome,
            fragment.id
        );
        for (offset, base) in fragment.sequence.iter().copied().enumerate() {
            let index = (fragment.local_start - local_start) as usize + offset;
            if let Some(existing) = sequence[index] {
                anyhow::ensure!(
                    existing == base,
                    "conflicting overlap for locus {:?}, source {:?}, chromosome {:?} at local position {} (record {:?})",
                    key.locus,
                    key.source,
                    key.chromosome,
                    local_start + index as i64,
                    fragment.id
                );
            } else {
                sequence[index] = Some(base);
            }
        }
    }

    let sequence = sequence
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| anyhow::anyhow!("internal gap while merging locus {:?}", key.locus))?;
    let global_start = local_start
        .checked_add(coordinate_offset)
        .ok_or_else(|| anyhow::anyhow!("global coordinate overflow for locus {:?}", key.locus))?;
    let global_end = local_end
        .checked_add(coordinate_offset)
        .ok_or_else(|| anyhow::anyhow!("global coordinate overflow for locus {:?}", key.locus))?;
    Ok(MergedProbe {
        local_start,
        local_end,
        global_start,
        global_end,
        sequence,
    })
}

fn longest_contiguous(mut fragments: Vec<Fragment>, key: &GroupKey) -> anyhow::Result<MergedProbe> {
    fragments.sort_by(|left, right| {
        (left.local_start, left.local_end, &left.id).cmp(&(
            right.local_start,
            right.local_end,
            &right.id,
        ))
    });

    let mut components: Vec<Vec<Fragment>> = Vec::new();
    for fragment in fragments {
        match components.last_mut() {
            Some(component)
                if fragment.local_start
                    <= component
                        .iter()
                        .map(|item| item.local_end)
                        .max()
                        .unwrap_or(fragment.local_start) =>
            {
                component.push(fragment);
            }
            _ => components.push(vec![fragment]),
        }
    }

    components.sort_by(|left, right| {
        let span = |component: &[Fragment]| {
            let start = component
                .iter()
                .map(|fragment| fragment.local_start)
                .min()
                .unwrap_or(0);
            let end = component
                .iter()
                .map(|fragment| fragment.local_end)
                .max()
                .unwrap_or(0);
            end.saturating_sub(start)
        };
        span(right).cmp(&span(left)).then_with(|| {
            left.first()
                .map(|fragment| (&fragment.id, fragment.local_start))
                .cmp(
                    &right
                        .first()
                        .map(|fragment| (&fragment.id, fragment.local_start)),
                )
        })
    });

    let mut first_error = None;
    for component in &components {
        match merge_component(component, key) {
            Ok(merged) => return Ok(merged),
            Err(error) if first_error.is_none() => first_error = Some(error),
            Err(_) => {}
        }
    }
    Err(first_error.unwrap_or_else(|| anyhow::anyhow!("no probes found for locus {:?}", key.locus)))
}

pub fn run(input: &Path, output: &Path) -> anyhow::Result<()> {
    let records = phyluce_io::read_fasta(input)
        .with_context(|| format!("reading probe FASTA {}", input.display()))?;
    let mut groups: BTreeMap<GroupKey, Vec<Fragment>> = BTreeMap::new();
    for record in records {
        let (key, fragment) = parse_fragment(record)?;
        groups.entry(key).or_default().push(fragment);
    }

    let mut loci: BTreeMap<String, Vec<(GroupKey, MergedProbe)>> = BTreeMap::new();
    for (key, fragments) in groups {
        let merged = longest_contiguous(fragments, &key)?;
        loci.entry(key.locus.clone())
            .or_default()
            .push((key, merged));
    }

    let output_paths: BTreeMap<String, PathBuf> = loci
        .keys()
        .map(|locus| {
            let filename = format!("{locus}.fasta");
            crate::output_path::output_file(output, &filename)
                .map(|path| (locus.clone(), path))
                .with_context(|| format!("building output filename for locus {locus:?}"))
        })
        .collect::<anyhow::Result<_>>()?;
    crate::output_path::prepare_output_dir(output)?;

    let mut reference_count = 0usize;
    let mut longest = 0usize;
    for (locus, references) in &loci {
        let mut text = String::new();
        for (index, (key, merged)) in references.iter().enumerate() {
            let probe_number = index + 1;
            let id = format!("{locus}_p{probe_number}");
            text.push_str(&format!(
                ">{id} |design:{},designer:{},probes-locus:{locus},probes-probe:{probe_number},probes-source:{},probes-global-chromo:{},probes-global-start:{},probes-global-end:{},probes-local-start:{},probes-local-end:{}\n{}\n",
                key.design,
                key.designer,
                key.source,
                key.chromosome,
                merged.global_start,
                merged.global_end,
                merged.local_start,
                merged.local_end,
                String::from_utf8_lossy(&merged.sequence)
            ));
            reference_count += 1;
            longest = longest.max(merged.sequence.len());
        }
        crate::output_path::write_atomic(&output_paths[locus], text)?;
    }

    crate::cli_warn!(
        "Wrote {reference_count} full probes across {} loci to {} (longest: {longest} bp)",
        loci.len(),
        output.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fragment(id: &str, start: i64, sequence: &[u8]) -> Fragment {
        Fragment {
            id: id.to_string(),
            local_start: start,
            local_end: start + sequence.len() as i64,
            global_start: 1000 + start,
            sequence: sequence.to_vec(),
        }
    }

    fn key() -> GroupKey {
        GroupKey {
            locus: "uce-1".to_string(),
            source: "taxon".to_string(),
            chromosome: "chr1".to_string(),
            design: "test".to_string(),
            designer: "tester".to_string(),
        }
    }

    #[test]
    fn merges_overlapping_probes_into_the_longest_contiguous_sequence() {
        let full = [vec![b'A'; 40], vec![b'C'; 80], vec![b'G'; 40]].concat();
        let fragments = vec![
            fragment("uce-1_p1", 10, &full[..120]),
            fragment("uce-1_p2", 50, &full[40..]),
            fragment("uce-1_p3", 500, &[b'T'; 120]),
        ];
        let merged = longest_contiguous(fragments, &key()).unwrap();
        assert_eq!(merged.local_start, 10);
        assert_eq!(merged.local_end, 170);
        assert_eq!(merged.sequence, full);
    }

    #[test]
    fn rejects_conflicting_overlap() {
        let fragments = vec![
            fragment("uce-1_p1", 0, &[b'A'; 120]),
            fragment("uce-1_p2", 40, &[b'C'; 120]),
        ];
        assert!(longest_contiguous(fragments, &key()).is_err());
    }

    #[test]
    fn ignores_an_invalid_shorter_component() {
        let full = [vec![b'A'; 80], vec![b'C'; 80]].concat();
        let fragments = vec![
            fragment("uce-1_p1", 0, &full[..120]),
            fragment("uce-1_p2", 40, &full[40..]),
            fragment("uce-1_p3", 500, &[b'G'; 80]),
            fragment("uce-1_p4", 540, &[b'T'; 80]),
        ];
        let merged = longest_contiguous(fragments, &key()).unwrap();
        assert_eq!(merged.sequence, full);
    }

    #[test]
    fn falls_back_to_the_longest_valid_component() {
        let fragments = vec![
            fragment("uce-1_p1", 0, &[b'A'; 120]),
            fragment("uce-1_p2", 40, &[b'C'; 120]),
            fragment("uce-1_p3", 500, &[b'G'; 120]),
        ];
        let merged = longest_contiguous(fragments, &key()).unwrap();
        assert_eq!(merged.sequence, vec![b'G'; 120]);
    }

    #[test]
    fn rejects_inconsistent_coordinate_offsets() {
        let mut second = fragment("uce-1_p2", 40, &[b'A'; 120]);
        second.global_start += 10;
        let fragments = vec![fragment("uce-1_p1", 0, &[b'A'; 120]), second];
        assert!(longest_contiguous(fragments, &key()).is_err());
    }

    #[test]
    fn rejects_mismatched_global_span_and_ambiguous_bases() {
        let base_header = "design:test,designer:tester,probes-locus:uce-1,probes-probe:1,probes-source:taxon,probes-global-chromo:chr1,probes-global-start:1000,probes-global-end:1120,probes-local-start:0,probes-local-end:120";
        let mismatched_span = FastaRecord {
            id: "uce-1_p1".to_string(),
            description: format!("uce-1_p1 |{base_header}")
                .replace("probes-global-end:1120", "probes-global-end:1121"),
            sequence: "A".repeat(120),
        };
        assert!(parse_fragment(mismatched_span).is_err());

        let ambiguous = FastaRecord {
            id: "uce-1_p1".to_string(),
            description: format!("uce-1_p1 |{base_header}"),
            sequence: "A".repeat(119) + "N",
        };
        assert!(parse_fragment(ambiguous).is_err());
    }

    #[test]
    fn writes_one_original_style_fasta_per_locus() {
        let directory =
            std::env::temp_dir().join(format!("phyluce-full-probes-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&directory);
        std::fs::create_dir_all(&directory).unwrap();
        let input = directory.join("probes.fasta");
        let output = directory.join("output");
        let full = "A".repeat(40) + &"C".repeat(80) + &"G".repeat(40);
        std::fs::write(
            &input,
            format!(
                ">uce-1_p1 |design:test,designer:tester,probes-locus:uce-1,probes-probe:1,probes-source:taxon,probes-global-chromo:chr1,probes-global-start:1010,probes-global-end:1130,probes-local-start:10,probes-local-end:130\n{}\n>uce-1_p2 |design:test,designer:tester,probes-locus:uce-1,probes-probe:2,probes-source:taxon,probes-global-chromo:chr1,probes-global-start:1050,probes-global-end:1170,probes-local-start:50,probes-local-end:170\n{}\n",
                &full[..120],
                &full[40..]
            ),
        )
        .unwrap();

        run(&input, &output).unwrap();

        let result = std::fs::read_to_string(output.join("uce-1.fasta")).unwrap();
        assert!(result.starts_with(
            ">uce-1_p1 |design:test,designer:tester,probes-locus:uce-1,probes-probe:1,"
        ));
        assert!(result.contains("probes-source:taxon,probes-global-chromo:chr1"));
        assert!(result.ends_with(&format!("{full}\n")));
        std::fs::remove_dir_all(directory).unwrap();
    }
}
