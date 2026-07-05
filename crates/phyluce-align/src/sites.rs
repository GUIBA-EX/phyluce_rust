//! Informative-site / difference counting, mirroring `phyluce/sites.py`.

use std::collections::HashMap;

use crate::Alignment;

/// Per-column character counts (uppercased), used by both site-counting
/// functions below. Mirrors `Counter(align[:, idx].upper())`.
fn column_counts(alignment: &Alignment, col: usize) -> HashMap<u8, usize> {
    let mut counts = HashMap::new();
    for row in &alignment.rows {
        let c = row.seq[col].to_ascii_uppercase();
        *counts.entry(c).or_insert(0) += 1;
    }
    counts
}

/// Mirrors `get_informative_sites`: after removing gap/N/? counts, a
/// column is "informative" if at least 2 distinct remaining characters
/// each occur at least twice.
fn is_informative(mut counts: HashMap<u8, usize>) -> bool {
    counts.remove(&b'-');
    counts.remove(&b'N');
    counts.remove(&b'?');
    if counts.len() >= 2 {
        let sufficient = counts.values().filter(|&&v| v >= 2).count();
        if sufficient >= 2 {
            return true;
        }
    }
    false
}

/// Mirrors `get_differences`: returns (counted, differs).
fn differences(mut counts: HashMap<u8, usize>) -> (bool, bool) {
    counts.remove(&b'-');
    counts.remove(&b'N');
    counts.remove(&b'?');
    counts.remove(&b'X');
    let sufficient_sites = counts.len();
    if sufficient_sites >= 2 {
        (true, true)
    } else if sufficient_sites >= 1 && counts.values().max().copied().unwrap_or(0) > 1 {
        (true, false)
    } else {
        (false, false)
    }
}

/// Mirrors `compute_informative_sites`: (sum_informative_sites,
/// sum_differences, sum_counted_sites) across every column.
pub fn compute_informative_sites(alignment: &Alignment) -> (usize, usize, usize) {
    let ncols = alignment.nchar();
    let mut informative = 0usize;
    let mut diffs = 0usize;
    let mut counted = 0usize;
    for col in 0..ncols {
        let counts = column_counts(alignment, col);
        if is_informative(counts.clone()) {
            informative += 1;
        }
        let (is_counted, differs) = differences(counts);
        if is_counted {
            counted += 1;
            if differs {
                diffs += 1;
            }
        }
    }
    (informative, diffs, counted)
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
}
