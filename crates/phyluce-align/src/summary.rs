//! Per-alignment summary statistics mirroring `phyluce/summary.py`'s
//! `get_stats`/`get_characters` (the per-locus half used by
//! `phyluce_align_get_align_summary_data`'s `--output-stats` CSV).

use std::collections::HashMap;

use crate::sites::compute_site_statistics;
use crate::Alignment;

#[derive(Debug, Clone)]
pub struct AlignSummary {
    pub length: usize,
    pub sum_informative_sites: usize,
    pub sum_differences: usize,
    pub sum_counted_sites: usize,
    /// Uppercased character -> count, across every row and column.
    pub characters: HashMap<u8, usize>,
}

impl AlignSummary {
    pub fn char_count(&self, c: u8) -> usize {
        self.characters.get(&c).copied().unwrap_or(0)
    }

    /// Mirrors `round(sum([G, C]) / sum(v for k, v where k != '-') * 100, 2)`.
    pub fn gc_content_percent(&self) -> f64 {
        let gc = self.char_count(b'G') + self.char_count(b'C');
        let denom: usize = self
            .characters
            .iter()
            .filter(|(&c, _)| c != b'-')
            .map(|(_, v)| v)
            .sum();
        if denom == 0 {
            return 0.0;
        }
        let raw = gc as f64 / denom as f64 * 100.0;
        (raw * 100.0).round() / 100.0
    }
}

/// Mirrors `summary.get_stats` (minus the file/name bookkeeping, which the
/// CLI layer handles).
pub fn compute_align_summary(alignment: &Alignment) -> AlignSummary {
    let statistics = compute_site_statistics(alignment);
    let characters = statistics
        .characters
        .iter()
        .enumerate()
        .filter(|(_, count)| **count != 0)
        .map(|(character, &count)| (character as u8, count))
        .collect();
    AlignSummary {
        length: alignment.nchar(),
        sum_informative_sites: statistics.informative,
        sum_differences: statistics.differences,
        sum_counted_sites: statistics.counted,
        characters,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_content_matches_expected_rounding() {
        let a = Alignment::from_pairs(vec![
            ("a".to_string(), "GGCC".to_string()),
            ("b".to_string(), "AATT".to_string()),
        ]);
        let s = compute_align_summary(&a);
        assert_eq!(s.gc_content_percent(), 50.0);
    }
}
