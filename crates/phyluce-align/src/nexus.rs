//! NEXUS writer matching `Bio.AlignIO`'s nexus format (via `Bio.Nexus.Nexus`),
//! field-for-field: quoting rule, name-column padding, and the
//! `interleave` switch/block width.

use crate::{Alignment, AlignmentRow};

#[derive(Debug, thiserror::Error)]
pub enum NexusError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("no 'matrix' block found in NEXUS file")]
    NoMatrixBlock,
    #[error("empty taxon label on a matrix line: {0:?}")]
    EmptyLabel(String),
}

/// Characters that force a taxon label to be single-quoted (with internal
/// single quotes doubled), matching `Bio.Nexus.Nexus.safename` using
/// `WHITESPACE + PUNCTUATION`.
const SPECIAL_CHARS: &str = " \t\n()[]{}\\,;:=*'\"`+-<>";

pub fn safename(name: &str) -> String {
    let escaped = name.replace('\'', "''");
    if escaped.chars().any(|c| SPECIAL_CHARS.contains(c)) {
        format!("'{escaped}'")
    } else {
        escaped
    }
}

/// Render an alignment as NEXUS, matching `format(alignment, "nexus")`
/// (`Bio.AlignIO`'s writer): interleaved in 70-column blocks when
/// `nchar > 1000`, single block otherwise.
pub fn format_nexus(alignment: &Alignment) -> String {
    format_nexus_with_interleave(alignment, alignment.nchar() > 1000)
}

/// Render an alignment as NEXUS with an explicit interleave choice,
/// matching `Nexus.write_nexus_data(..., interleave=...)` called directly
/// (as `phyluce_align_concatenate_alignments` does, bypassing `Bio.AlignIO`'s
/// `columns > 1000` auto-detection -- its default is always
/// non-interleaved regardless of length). No character sets/codon blocks
/// appended here; callers needing them (e.g. concatenation) append
/// separately.
pub fn format_nexus_with_interleave(alignment: &Alignment, interleave: bool) -> String {
    let ntax = alignment.ntax();
    let nchar = alignment.nchar();
    let quoted: Vec<String> = alignment.rows.iter().map(|r| safename(&r.id)).collect();
    let namelength = quoted.iter().map(|s| s.chars().count()).max().unwrap_or(0);

    let mut out = String::new();
    out.push_str("#NEXUS\nbegin data;\n");
    out.push_str(&format!("dimensions ntax={ntax} nchar={nchar};\n"));
    if interleave {
        out.push_str("format datatype=dna missing=? gap=- interleave;\n");
    } else {
        out.push_str("format datatype=dna missing=? gap=-;\n");
    }
    out.push_str("matrix\n");

    let blocksize = if interleave { 70 } else { nchar.max(1) };
    if nchar > 0 {
        let mut seek = 0;
        while seek < nchar {
            let end = (seek + blocksize).min(nchar);
            for (i, row) in alignment.rows.iter().enumerate() {
                out.push_str(&quoted[i]);
                let pad = namelength + 1 - quoted[i].chars().count();
                out.push_str(&" ".repeat(pad));
                out.push_str(std::str::from_utf8(&row.seq[seek..end]).unwrap_or(""));
                out.push('\n');
            }
            if interleave {
                out.push('\n');
            }
            seek = end;
        }
    }
    out.push_str(";\nend;\n");
    out
}

/// Parse a NEXUS alignment (interleaved or not, quoted or bare taxon
/// labels), reconstructing each taxon's full sequence by concatenating its
/// chunks across blocks in the order first encountered.
pub fn parse_nexus(text: &str) -> Result<Alignment, NexusError> {
    let lower = text.to_lowercase();
    let matrix_pos = lower.find("matrix").ok_or(NexusError::NoMatrixBlock)?;
    let after_matrix = &text[matrix_pos + "matrix".len()..];

    let mut order: Vec<String> = Vec::new();
    let mut seqs: std::collections::HashMap<String, String> = std::collections::HashMap::new();

    for raw_line in after_matrix.lines() {
        let line = raw_line.trim_end();
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == ";" {
            break;
        }

        let (label, rest) = if let Some(stripped) = trimmed.strip_prefix('\'') {
            // find the closing quote, treating '' as an escaped literal '
            let bytes = stripped.as_bytes();
            let mut i = 0usize;
            let mut label = String::new();
            loop {
                if i >= bytes.len() {
                    break;
                }
                if bytes[i] == b'\'' {
                    if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                        label.push('\'');
                        i += 2;
                        continue;
                    } else {
                        i += 1;
                        break;
                    }
                }
                label.push(bytes[i] as char);
                i += 1;
            }
            (label, stripped[i..].trim_start())
        } else {
            match trimmed.split_once(char::is_whitespace) {
                Some((label, rest)) => (label.to_string(), rest.trim_start()),
                None => (trimmed.to_string(), ""),
            }
        };

        if label.is_empty() {
            return Err(NexusError::EmptyLabel(line.to_string()));
        }
        if !seqs.contains_key(&label) {
            order.push(label.clone());
        }
        seqs.entry(label).or_default().push_str(rest);
    }

    Ok(Alignment {
        rows: order
            .into_iter()
            .map(|id| {
                let seq = seqs.remove(&id).unwrap_or_default().into_bytes();
                AlignmentRow { id, seq }
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quotes_names_with_hyphens_only() {
        assert_eq!(safename("uce-1_short"), "'uce-1_short'");
        assert_eq!(safename("a"), "a");
    }

    #[test]
    fn parses_bare_non_interleaved_names() {
        let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\nformat datatype=dna missing=? gap=-;\nmatrix\ntaxon_a ACGT\ntaxon_b ACGA\n;\nend;\n";
        let aln = parse_nexus(text).unwrap();
        assert_eq!(aln.rows.len(), 2);
        assert_eq!(aln.rows[0].id, "taxon_a");
        assert_eq!(aln.rows[0].seq, b"ACGT");
        assert_eq!(aln.rows[1].seq, b"ACGA");
    }

    #[test]
    fn parses_quoted_and_interleaved_names() {
        let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=8;\nformat datatype=dna missing=? gap=- interleave;\nmatrix\n'uce-1_a' ACGT\n'uce-1_b' ACGA\n\n'uce-1_a' TTTT\n'uce-1_b' GGGG\n\n;\nend;\n";
        let aln = parse_nexus(text).unwrap();
        assert_eq!(aln.rows.len(), 2);
        assert_eq!(aln.rows[0].id, "uce-1_a");
        assert_eq!(aln.rows[0].seq, b"ACGTTTTT");
        assert_eq!(aln.rows[1].id, "uce-1_b");
        assert_eq!(aln.rows[1].seq, b"ACGAGGGG");
    }

    #[test]
    fn round_trips_through_format_and_parse() {
        let a = Alignment::from_pairs(vec![
            ("uce-1_short".to_string(), "ACGT-ACGT".to_string()),
            ("a".to_string(), "ACGTAACGT".to_string()),
        ]);
        let text = format_nexus(&a);
        let parsed = parse_nexus(&text).unwrap();
        assert_eq!(parsed, a);
    }

    #[test]
    fn non_interleaved_matches_biopython_reference() {
        let a = Alignment::from_pairs(vec![
            ("uce-1_short".to_string(), "ACGT-ACGT".to_string()),
            ("a".to_string(), "ACGTAACGT".to_string()),
        ]);
        let expected = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=9;\nformat datatype=dna missing=? gap=-;\nmatrix\n'uce-1_short' ACGT-ACGT\na             ACGTAACGT\n;\nend;\n";
        assert_eq!(format_nexus(&a), expected);
    }
}
