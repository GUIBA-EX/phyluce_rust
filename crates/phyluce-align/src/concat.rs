//! Multi-locus concatenation mirroring `Bio.Nexus.Nexus.combine` (used by
//! `phyluce_align_concatenate_alignments`): taxa missing from a given locus
//! are padded with `?` for that locus's length; taxa appearing for the
//! first time in a later locus are padded with `?` for every column seen
//! so far and appended to the end of the taxon order.

use crate::{Alignment, AlignmentRow};

/// `HashMap`/`HashSet` keyed on `ahash` instead of the standard library's
/// SipHash. Not the point of the rewrite below -- see `concatenate`'s doc
/// comment -- but free once the linear scans are gone anyway.
type FastMap<K, V> = std::collections::HashMap<K, V, ahash::RandomState>;
type FastSet<T> = std::collections::HashSet<T, ahash::RandomState>;

/// A locus's column range within the concatenated matrix (0-indexed,
/// half-open, matching Python's `range(start, stop)`).
pub struct Charset {
    pub name: String,
    pub start: usize,
    pub stop: usize,
}

/// Mirrors `Nexus.combine`: `files` must already be sorted the way the
/// caller wants columns ordered (phyluce sorts by basename).
///
/// Builds one `FastMap` per locus instead of the `Vec::iter().find()` (scan
/// of that locus's rows, once per taxon) and `Vec::contains()` (scan of
/// every taxon seen so far, once per row in this locus) an earlier version
/// used -- those made the whole function effectively O(taxa^2) per locus.
/// `bench_concatenate_scaling_with_taxon_count` showed the old version's
/// time roughly quadrupling (not doubling) each time the taxon count
/// doubled, confirming the shape before this rewrite.
pub fn concatenate(files: &[(String, Alignment)]) -> (Alignment, Vec<Charset>) {
    if files.is_empty() {
        return (Alignment::default(), Vec::new());
    }

    let (first_name, first_aln) = &files[0];
    let mut taxa_order: Vec<String> = first_aln.rows.iter().map(|r| r.id.clone()).collect();
    let mut taxa_seen: FastSet<String> = taxa_order.iter().cloned().collect();
    let mut seqs: FastMap<String, Vec<u8>> = first_aln
        .rows
        .iter()
        .map(|r| (r.id.clone(), r.seq.clone()))
        .collect();
    let mut nchar = first_aln.nchar();
    let mut charsets = vec![Charset {
        name: first_name.clone(),
        start: 0,
        stop: nchar,
    }];

    for (name, aln) in &files[1..] {
        // `entry().or_insert()`, not `.collect()`: a locus with duplicate
        // taxon IDs (malformed input) should keep the *first* matching
        // row, matching the old `Vec::iter().find()` this replaced --
        // `.collect()` into a HashMap keeps whichever insert happens last.
        let mut locus_index: FastMap<&str, &AlignmentRow> = FastMap::default();
        for r in &aln.rows {
            locus_index.entry(r.id.as_str()).or_insert(r);
        }

        // Taxa already in the matrix: extend with this locus's sequence
        // (`both`, from the old naming) or pad with `?` if this locus
        // doesn't have them (`combined_only`) -- an O(1) lookup per taxon
        // instead of an O(locus rows) scan.
        for t in &taxa_order {
            match locus_index.get(t.as_str()) {
                Some(row) => seqs.get_mut(t).unwrap().extend_from_slice(&row.seq),
                None => seqs
                    .get_mut(t)
                    .unwrap()
                    .extend(std::iter::repeat_n(b'?', aln.nchar())),
            }
        }
        // Taxa new to this locus (`m_only`, from the old naming), in the
        // order they appear in `aln.rows` (matching the old `Vec`-filter's
        // order) -- `taxa_seen.insert` returning `true` means "not seen
        // before", an O(1) check instead of `taxa_order.contains()`.
        for row in &aln.rows {
            if taxa_seen.insert(row.id.clone()) {
                let mut padded: Vec<u8> = std::iter::repeat_n(b'?', nchar).collect();
                padded.extend_from_slice(&row.seq);
                seqs.insert(row.id.clone(), padded);
                taxa_order.push(row.id.clone());
            }
        }

        charsets.push(Charset {
            name: name.clone(),
            start: nchar,
            stop: nchar + aln.nchar(),
        });
        nchar += aln.nchar();
    }

    let rows = taxa_order
        .into_iter()
        .map(|id| {
            let seq = seqs.remove(&id).unwrap_or_default();
            AlignmentRow { id, seq }
        })
        .collect();
    (Alignment { rows }, charsets)
}

/// Render the `begin sets; ... end;` block matching
/// `Nexus.append_sets()`'s output for a single, all-contiguous
/// charpartition named `combined` (phyluce never uses codon partitions or
/// taxon sets).
pub fn format_sets_block(charsets: &[Charset]) -> String {
    if charsets.is_empty() {
        return String::new();
    }
    let mut out = String::from("\nbegin sets;\n");
    for c in charsets {
        out.push_str(&format!(
            "charset {} = {}-{};\n",
            crate::nexus::safename(&c.name),
            c.start + 1,
            c.stop
        ));
    }
    out.push_str("charpartition combined = ");
    // Note: unlike the `charset` lines above, `append_sets` does NOT run
    // partition names through `safename()` here -- reproduced as-is.
    let parts: Vec<String> = charsets
        .iter()
        .map(|c| format!("{}: {}-{}", c.name, c.start + 1, c.stop))
        .collect();
    out.push_str(&parts.join(", "));
    out.push_str(";\nend;\n");
    out
}

/// Mirrors `Nexus.export_phylip`: relaxed PHYLIP (untruncated names).
pub fn format_phylip(alignment: &Alignment) -> String {
    let mut out = format!("{} {}\n", alignment.ntax(), alignment.nchar());
    for row in &alignment.rows {
        out.push_str(&crate::nexus::safename(&row.id));
        out.push(' ');
        out.push_str(std::str::from_utf8(&row.seq).unwrap_or(""));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ad hoc benchmark for `concatenate`'s per-locus loop, which does an
    // `aln.rows.iter().find(...)` (linear scan over that locus's rows) per
    // taxon, plus `taxa_order.contains(t)` (linear scan over all taxa seen
    // so far) per new taxon -- both look like they could be O(taxa) inside
    // an O(taxa)-per-locus loop, i.e. O(taxa^2) per locus, instead of the
    // O(taxa) a HashMap lookup would give. Run with:
    //   cargo +stable test --release -p phyluce-align --lib -- --ignored --nocapture bench_
    //
    // Synthetic workload: 2000 loci x 100 taxa, ~90% present per locus
    // (some missing/new-taxon padding to exercise all three code paths),
    // 500bp sequences -- roughly a mid-size UCE phylogenomic dataset.
    fn synthetic_loci(n_loci: usize, n_taxa: usize, seq_len: usize) -> Vec<(String, Alignment)> {
        let taxa: Vec<String> = (0..n_taxa).map(|i| format!("taxon_{i}")).collect();
        let mut out = Vec::with_capacity(n_loci);
        for locus in 0..n_loci {
            let seq: Vec<u8> = (0..seq_len).map(|i| b"ACGT"[(locus + i) % 4]).collect();
            let rows: Vec<AlignmentRow> = taxa
                .iter()
                .enumerate()
                // Drop ~10% of taxa per locus (deterministically) so both
                // `combined_only` and `m_only` padding paths get exercised,
                // not just the fully-shared `both` path.
                .filter(|(i, _)| (i + locus) % 10 != 0)
                .map(|(_, id)| AlignmentRow {
                    id: id.clone(),
                    seq: seq.clone(),
                })
                .collect();
            out.push((format!("locus_{locus:04}.nexus"), Alignment { rows }));
        }
        out
    }

    #[test]
    #[ignore]
    fn bench_concatenate_realistic_scale() {
        let files = synthetic_loci(2_000, 100, 500);

        let start = std::time::Instant::now();
        let (combined, charsets) = concatenate(&files);
        let elapsed = start.elapsed();

        eprintln!(
            "[bench] concatenate: {} loci x ~100 taxa -> {} taxa, {} columns, {} charsets in {:?}",
            files.len(),
            combined.ntax(),
            combined.nchar(),
            charsets.len(),
            elapsed
        );
    }

    #[test]
    #[ignore]
    fn bench_concatenate_scaling_with_taxon_count() {
        // The per-locus loop's `aln.rows.iter().find(...)` and
        // `taxa_order.contains(t)` are both linear scans -- if that's
        // really an O(taxa) operation happening O(taxa) times per locus,
        // doubling the taxon count should roughly *quadruple* the time
        // (holding loci count and sequence length fixed), not just double
        // it. This checks that shape directly instead of guessing from a
        // single data point.
        for n_taxa in [100usize, 200, 400, 800] {
            let files = synthetic_loci(500, n_taxa, 500);
            let start = std::time::Instant::now();
            let _ = concatenate(&files);
            let elapsed = start.elapsed();
            eprintln!(
                "[bench] concatenate: 500 loci x {n_taxa} taxa in {:?} ({:.2} ms/locus)",
                elapsed,
                elapsed.as_secs_f64() * 1000.0 / 500.0
            );
        }
    }

    #[test]
    fn pads_taxa_missing_from_a_locus() {
        let a = Alignment::from_pairs(vec![
            ("x".to_string(), "AAAA".to_string()),
            ("y".to_string(), "CCCC".to_string()),
        ]);
        let b = Alignment::from_pairs(vec![("x".to_string(), "GG".to_string())]);
        let (combined, charsets) = concatenate(&[("a".to_string(), a), ("b".to_string(), b)]);
        let x = combined.rows.iter().find(|r| r.id == "x").unwrap();
        let y = combined.rows.iter().find(|r| r.id == "y").unwrap();
        assert_eq!(x.seq, b"AAAAGG");
        assert_eq!(y.seq, b"CCCC??");
        assert_eq!(charsets.len(), 2);
        assert_eq!((charsets[0].start, charsets[0].stop), (0, 4));
        assert_eq!((charsets[1].start, charsets[1].stop), (4, 6));
    }

    #[test]
    fn pads_and_appends_taxa_new_to_a_later_locus() {
        let a = Alignment::from_pairs(vec![("x".to_string(), "AAAA".to_string())]);
        let b = Alignment::from_pairs(vec![
            ("x".to_string(), "GG".to_string()),
            ("z".to_string(), "TT".to_string()),
        ]);
        let (combined, _) = concatenate(&[("a".to_string(), a), ("b".to_string(), b)]);
        assert_eq!(combined.rows.last().unwrap().id, "z");
        assert_eq!(combined.rows.last().unwrap().seq, b"????TT");
    }

    #[test]
    fn duplicate_taxon_id_in_a_locus_keeps_the_first_row() {
        let a = Alignment::from_pairs(vec![("x".to_string(), "AA".to_string())]);
        let b = Alignment::from_pairs(vec![
            ("x".to_string(), "CC".to_string()),
            ("x".to_string(), "GG".to_string()),
        ]);
        let (combined, _) = concatenate(&[("a".to_string(), a), ("b".to_string(), b)]);
        let x = combined.rows.iter().find(|r| r.id == "x").unwrap();
        assert_eq!(x.seq, b"AACC");
    }

    #[test]
    fn formats_sets_block() {
        let charsets = vec![
            Charset {
                name: "uce-1.nexus".to_string(),
                start: 0,
                stop: 459,
            },
            Charset {
                name: "uce-2.nexus".to_string(),
                start: 459,
                stop: 1068,
            },
        ];
        let out = format_sets_block(&charsets);
        assert_eq!(
            out,
            "\nbegin sets;\ncharset 'uce-1.nexus' = 1-459;\ncharset 'uce-2.nexus' = 460-1068;\ncharpartition combined = uce-1.nexus: 1-459, uce-2.nexus: 460-1068;\nend;\n"
        );
    }
}
