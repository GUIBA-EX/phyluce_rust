//! Parsing of a NEXUS `begin sets;` block's `charset` lines, mirroring
//! `phyluce_align_split_concat_nexus_to_loci`'s use of `Bio.Nexus.Nexus`'s
//! `charsets` dict (built from `charset 'name' = start-end;` lines; any
//! `charpartition` line is skipped by the legacy script before parsing).

use crate::Alignment;

/// One `charset name = start-end;` entry, 0-indexed half-open range
/// (`stop` is exclusive), matching `Alignment::rows`/`nchar` indexing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Charset {
    pub name: String,
    pub start: usize,
    pub stop: usize,
}

/// Reverse of `nexus::safename`: strip a single layer of surrounding
/// quotes (if present) and un-double any escaped `''`.
fn unquote_name(raw: &str) -> String {
    if let Some(inner) = raw.strip_prefix('\'').and_then(|s| s.strip_suffix('\'')) {
        inner.replace("''", "'")
    } else {
        raw.to_string()
    }
}

/// Parse every `charset` line, skipping `charpartition` lines. Mirrors
/// reading a `Nexus.Nexus` object's `.charsets`.
pub fn parse_charsets(text: &str) -> Vec<Charset> {
    let mut charsets = Vec::new();
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if !line.starts_with("charset ") {
            continue;
        }
        let rest = line
            .trim_start_matches("charset ")
            .trim_end_matches(';')
            .trim();
        let Some((name_part, range_part)) = rest.split_once('=') else {
            continue;
        };
        let name = unquote_name(name_part.trim());
        let Some((start_s, stop_s)) = range_part.trim().split_once('-') else {
            continue;
        };
        let (Ok(start), Ok(stop)) = (
            start_s.trim().parse::<usize>(),
            stop_s.trim().parse::<usize>(),
        ) else {
            continue;
        };
        charsets.push(Charset {
            name,
            start: start - 1,
            stop,
        });
    }
    charsets
}

/// Slice an alignment's columns to `[start, stop)` (0-indexed, half-open),
/// keeping every row (callers typically drop empty rows afterwards).
pub fn slice_alignment(alignment: &Alignment, start: usize, stop: usize) -> Alignment {
    Alignment {
        rows: alignment
            .rows
            .iter()
            .map(|r| crate::AlignmentRow {
                id: r.id.clone(),
                seq: r.seq[start.min(r.seq.len())..stop.min(r.seq.len())].to_vec(),
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_charset_lines_and_skips_charpartition() {
        let text = "\nbegin sets;\ncharset 'uce-1.nexus' = 1-459;\ncharset 'uce-2.nexus' = 460-1068;\ncharpartition combined = uce-1.nexus: 1-459, uce-2.nexus: 460-1068;\nend;\n";
        let charsets = parse_charsets(text);
        assert_eq!(charsets.len(), 2);
        assert_eq!(charsets[0].name, "uce-1.nexus");
        assert_eq!((charsets[0].start, charsets[0].stop), (0, 459));
        assert_eq!(charsets[1].name, "uce-2.nexus");
        assert_eq!((charsets[1].start, charsets[1].stop), (459, 1068));
    }

    #[test]
    fn slices_columns() {
        let a = Alignment::from_pairs(vec![("x".to_string(), "ACGTACGT".to_string())]);
        let s = slice_alignment(&a, 2, 5);
        assert_eq!(s.rows[0].seq, b"GTA");
    }
}
