//! Multi-locus concatenation mirroring `Bio.Nexus.Nexus.combine` (used by
//! `phyluce_align_concatenate_alignments`): taxa missing from a given locus
//! are padded with `?` for that locus's length; taxa appearing for the
//! first time in a later locus are padded with `?` for every column seen
//! so far and appended to the end of the taxon order.

use crate::{Alignment, AlignmentRow};

/// A locus's column range within the concatenated matrix (0-indexed,
/// half-open, matching Python's `range(start, stop)`).
pub struct Charset {
    pub name: String,
    pub start: usize,
    pub stop: usize,
}

/// Mirrors `Nexus.combine`: `files` must already be sorted the way the
/// caller wants columns ordered (phyluce sorts by basename).
pub fn concatenate(files: &[(String, Alignment)]) -> (Alignment, Vec<Charset>) {
    if files.is_empty() {
        return (Alignment::default(), Vec::new());
    }

    let (first_name, first_aln) = &files[0];
    let mut taxa_order: Vec<String> = first_aln.rows.iter().map(|r| r.id.clone()).collect();
    let mut seqs: std::collections::HashMap<String, Vec<u8>> = first_aln
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
        let m_taxa: std::collections::HashSet<&str> =
            aln.rows.iter().map(|r| r.id.as_str()).collect();
        let both: Vec<String> = taxa_order
            .iter()
            .filter(|t| m_taxa.contains(t.as_str()))
            .cloned()
            .collect();
        let combined_only: Vec<String> = taxa_order
            .iter()
            .filter(|t| !m_taxa.contains(t.as_str()))
            .cloned()
            .collect();
        let m_only: Vec<String> = aln
            .rows
            .iter()
            .map(|r| r.id.clone())
            .filter(|t| !taxa_order.contains(t))
            .collect();

        for t in &both {
            let row_seq = &aln.rows.iter().find(|r| &r.id == t).unwrap().seq;
            seqs.get_mut(t).unwrap().extend_from_slice(row_seq);
        }
        for t in &combined_only {
            seqs.get_mut(t)
                .unwrap()
                .extend(std::iter::repeat_n(b'?', aln.nchar()));
        }
        for t in &m_only {
            let row_seq = &aln.rows.iter().find(|r| &r.id == t).unwrap().seq;
            let mut padded: Vec<u8> = std::iter::repeat_n(b'?', nchar).collect();
            padded.extend_from_slice(row_seq);
            seqs.insert(t.clone(), padded);
            taxa_order.push(t.clone());
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
