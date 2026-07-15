//! Informative-site / difference counting, mirroring `phyluce/sites.py`.

use crate::Alignment;

pub(crate) struct SiteStatistics {
    pub informative: usize,
    pub differences: usize,
    pub counted: usize,
    pub characters: [usize; 256],
}

/// Count all site metrics and alignment characters in one column-major pass.
/// A byte-indexed array avoids allocating and cloning a hash table per column.
pub(crate) fn compute_site_statistics(alignment: &Alignment) -> SiteStatistics {
    let mut statistics = SiteStatistics {
        informative: 0,
        differences: 0,
        counted: 0,
        characters: [0; 256],
    };
    let mut column = [0usize; 256];
    let mut observed = [0u8; 256];

    for col in 0..alignment.nchar() {
        let mut observed_count = 0usize;
        for row in &alignment.rows {
            let byte = row.seq[col].to_ascii_uppercase();
            let character = byte as usize;
            if column[character] == 0 {
                observed[observed_count] = byte;
                observed_count += 1;
            }
            column[character] += 1;
            statistics.characters[character] += 1;
        }

        let mut informative_states = 0usize;
        let mut difference_states = 0usize;
        let mut max_difference_count = 0usize;
        for &byte in &observed[..observed_count] {
            let count = column[byte as usize];
            if !matches!(byte, b'-' | b'N' | b'?') && count >= 2 {
                informative_states += 1;
            }
            if !matches!(byte, b'-' | b'N' | b'?' | b'X') {
                difference_states += 1;
                max_difference_count = max_difference_count.max(count);
            }
            column[byte as usize] = 0;
        }

        statistics.informative += usize::from(informative_states >= 2);
        if difference_states >= 2 {
            statistics.counted += 1;
            statistics.differences += 1;
        } else if difference_states == 1 && max_difference_count > 1 {
            statistics.counted += 1;
        }
    }
    statistics
}

/// Mirrors `compute_informative_sites`: (sum_informative_sites,
/// sum_differences, sum_counted_sites) across every column.
pub fn compute_informative_sites(alignment: &Alignment) -> (usize, usize, usize) {
    let statistics = compute_site_statistics(alignment);
    (
        statistics.informative,
        statistics.differences,
        statistics.counted,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn informative_column_requires_two_bases_twice_each() {
        let a = Alignment::from_pairs(vec![
            ("a".to_string(), "AA".to_string()),
            ("b".to_string(), "AA".to_string()),
            ("c".to_string(), "CC".to_string()),
            ("d".to_string(), "CC".to_string()),
        ]);
        let (informative, _, _) = compute_informative_sites(&a);
        assert_eq!(informative, 2);
    }

    #[test]
    fn non_informative_when_only_one_variant_repeats() {
        let a = Alignment::from_pairs(vec![
            ("a".to_string(), "A".to_string()),
            ("b".to_string(), "A".to_string()),
            ("c".to_string(), "A".to_string()),
            ("d".to_string(), "C".to_string()),
        ]);
        let (informative, diffs, counted) = compute_informative_sites(&a);
        assert_eq!(informative, 0);
        assert_eq!(counted, 1);
        assert_eq!(diffs, 1);
    }

    #[test]
    fn character_totals_are_collected_in_the_same_pass() {
        let a = Alignment::from_pairs(vec![
            ("a".to_string(), "a?-".to_string()),
            ("b".to_string(), "ANN".to_string()),
        ]);
        let statistics = compute_site_statistics(&a);
        assert_eq!(statistics.characters[b'A' as usize], 2);
        assert_eq!(statistics.characters[b'N' as usize], 2);
        assert_eq!(statistics.characters[b'?' as usize], 1);
        assert_eq!(statistics.characters[b'-' as usize], 1);
    }
}
