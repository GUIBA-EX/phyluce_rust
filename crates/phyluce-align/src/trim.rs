//! The native phyluce 3-stage edge-trimming algorithm, ported field-for-field
//! from `phyluce/generic_align.py`'s `running_average` / `stage_one_trimming`
//! / `stage_two_trimming` / `trim_alignment`.
//!
//! Faithfully reproduces two upstream quirks rather than "fixing" them,
//! since fixture-based golden tests were generated against the original
//! behavior:
//! - The `set(trim) != (["?"])` check in the Python is comparing a `set` to
//!   a `list`, which is always `True` -- it never actually filters
//!   anything. We simply omit that no-op check.
//! - In the start/end search loops, the loop variable is only converted
//!   from a "reversed-array index" to a "true index from the front" when
//!   the loop `break`s early; if the loop runs to completion without
//!   breaking, the raw (untransformed) index is used instead. Both stages
//!   preserve this asymmetry.

use crate::{Alignment, AlignmentRow};

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum TrimParameterError {
    #[error("trim window must be greater than zero")]
    ZeroWindow,
    #[error("{name} must be a finite value between 0 and 1 (inclusive), got {value}")]
    InvalidProportion { name: &'static str, value: f64 },
}

pub fn validate_trim_parameters(
    window_size: usize,
    proportion: f64,
    threshold: f64,
    max_divergence: f64,
) -> Result<(), TrimParameterError> {
    if window_size == 0 {
        return Err(TrimParameterError::ZeroWindow);
    }
    for (name, value) in [
        ("proportion", proportion),
        ("threshold", threshold),
        ("max divergence", max_divergence),
    ] {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(TrimParameterError::InvalidProportion { name, value });
        }
    }
    Ok(())
}

/// `round()`-half-to-even, matching Python 3's float rounding used for
/// `int(round(proportion * taxa, 0))`.
fn python_round(x: f64) -> i64 {
    let floor = x.floor();
    let diff = x - floor;
    let f = floor as i64;
    if (diff - 0.5).abs() < 1e-9 {
        if f % 2 == 0 {
            f
        } else {
            f + 1
        }
    } else if diff < 0.5 {
        f
    } else {
        f + 1
    }
}

/// Per-column "is this column good" classification used by stage one,
/// mirroring `running_average`'s `column_count` logic (case-sensitive,
/// gaps counted separately from the majority-base check).
fn compute_good_columns(alignment: &Alignment, majority: i64) -> Vec<bool> {
    let ncols = alignment.nchar();
    let mut good = Vec::with_capacity(ncols);
    for col in 0..ncols {
        let mut counts: std::collections::HashMap<u8, i64> = std::collections::HashMap::new();
        let mut gap_count = 0i64;
        for row in &alignment.rows {
            let c = row.seq[col];
            *counts.entry(c).or_insert(0) += 1;
            if c == b'-' {
                gap_count += 1;
            }
        }
        if gap_count <= majority {
            let max_count = counts
                .iter()
                .filter(|(&c, _)| c != b'-')
                .map(|(_, v)| *v)
                .max()
                .unwrap_or(0);
            good.push(max_count >= majority);
        } else {
            good.push(false);
        }
    }
    good
}

/// Mirrors the `start_clip` search loop in `running_average`: only
/// evaluates the window/threshold test at positions where the column
/// itself is "good".
fn find_start(good: &[bool], window_size: usize, threshold: f64) -> usize {
    let mut start_clip = 0usize;
    for i in 0..good.len() {
        start_clip = i;
        if good[i] {
            let end = (i + window_size).min(good.len());
            let window = &good[i..end];
            let proportion = window.iter().filter(|&&b| b).count() as f64 / window.len() as f64;
            if proportion > threshold {
                break;
            }
        }
    }
    start_clip
}

/// Mirrors the `end_clip` search: same shape as `find_start` but over the
/// reversed array, and only remapped to a "from the front" index when the
/// loop actually breaks (see module docs).
fn find_end(good: &[bool], window_size: usize, threshold: f64) -> usize {
    let reversed: Vec<bool> = good.iter().rev().copied().collect();
    let mut end_clip = 0usize;
    for i in 0..reversed.len() {
        end_clip = i;
        if reversed[i] {
            let e = (i + window_size).min(reversed.len());
            let window = &reversed[i..e];
            let proportion = window.iter().filter(|&&b| b).count() as f64 / window.len() as f64;
            if proportion >= threshold {
                end_clip = reversed.len() - i;
                break;
            }
        }
    }
    end_clip
}

/// Mirrors `_replace_ends`: convert the leading and trailing runs of `-`
/// into equal-length runs of `?` (the two substitutions are sequential, so
/// an all-gap sequence is fully consumed by the first one).
fn replace_ends(seq: &[u8]) -> Vec<u8> {
    let mut v = seq.to_vec();
    let lead = v.iter().take_while(|&&c| c == b'-').count();
    for b in v.iter_mut().take(lead) {
        *b = b'?';
    }
    let trail = v.iter().rev().take_while(|&&c| c == b'-').count();
    let len = v.len();
    for b in v.iter_mut().skip(len - trail) {
        *b = b'?';
    }
    v
}

fn is_all_gaps(seq: &[u8]) -> bool {
    !seq.is_empty() && seq.iter().all(|&c| c == b'-')
}

/// Mirrors `stage_one_trimming`.
fn stage_one(
    alignment: &Alignment,
    window_size: usize,
    proportion: f64,
    threshold: f64,
    min_len: usize,
    replace_ends_flag: bool,
) -> Option<Alignment> {
    let taxa = alignment.ntax();
    if taxa == 0 || alignment.nchar() == 0 {
        return None;
    }
    let majority = python_round(proportion * taxa as f64);
    let good = compute_good_columns(alignment, majority);
    let start = find_start(&good, window_size, threshold);
    let end = find_end(&good, window_size, threshold);
    if end == 0 {
        return None;
    }

    let mut out_rows = Vec::with_capacity(taxa);
    for row in &alignment.rows {
        let len = row.seq.len();
        let (s, e) = (start.min(len), end.min(len));
        let trim: Vec<u8> = if s < e {
            row.seq[s..e].to_vec()
        } else {
            Vec::new()
        };
        if !is_all_gaps(&trim) && trim.len() >= min_len {
            let final_seq = if replace_ends_flag {
                replace_ends(&trim)
            } else {
                trim
            };
            out_rows.push(AlignmentRow {
                id: row.id.clone(),
                seq: final_seq,
            });
        } else {
            return None;
        }
    }
    Some(Alignment { rows: out_rows })
}

/// Mirrors `_get_ends`: (leading gap run length, index just past the
/// trailing gap run).
fn get_ends(seq: &[u8]) -> (usize, usize) {
    let start_gap = seq.iter().take_while(|&&c| c == b'-').count();
    let end_gap = seq.iter().rev().take_while(|&&c| c == b'-').count();
    (start_gap, seq.len().saturating_sub(end_gap))
}

/// Mirrors `_alignment_consensus`'s per-column majority-base pick (ties
/// broken by first-seen order among uppercased column characters, as
/// `collections.Counter.most_common` does).
fn consensus_char(column: &[u8]) -> u8 {
    let mut order: Vec<u8> = Vec::new();
    let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
    for &c in column {
        *counts.entry(c).or_insert_with(|| {
            order.push(c);
            0
        }) += 1;
    }
    let mut best_char = order[0];
    let mut best_count = counts[&best_char];
    for &c in &order[1..] {
        let cnt = counts[&c];
        if cnt > best_count {
            best_count = cnt;
            best_char = c;
        }
    }
    if best_count as f64 / column.len() as f64 >= 0.5 {
        best_char
    } else {
        b'N'
    }
}

fn compute_consensus(alignment: &Alignment) -> Vec<u8> {
    let ncols = alignment.nchar();
    let mut result = Vec::with_capacity(ncols);
    for col in 0..ncols {
        let column: Vec<u8> = alignment
            .rows
            .iter()
            .map(|r| r.seq[col].to_ascii_uppercase())
            .collect();
        result.push(consensus_char(&column));
    }
    result
}

/// Mirrors the `bad_start` search in `stage_two_trimming` (unconditional
/// per-position evaluation, unlike stage one's `find_start`).
fn find_divergence_start(compare: &[bool], window_size: usize, max_divergence: f64) -> usize {
    if compare.is_empty() {
        return 0;
    }
    let mut bad_start = 0usize;
    for i in 0..compare.len() {
        bad_start = i;
        let e = (i + window_size).min(compare.len());
        let window = &compare[i..e];
        let divergence = window.iter().filter(|&&b| b).count() as f64 / window.len() as f64;
        if divergence < max_divergence {
            break;
        }
    }
    bad_start
}

/// Mirrors the `bad_end` search: same front/back-index remap-only-on-break
/// quirk as `find_end`.
fn find_divergence_end(compare: &[bool], window_size: usize, max_divergence: f64) -> usize {
    if compare.is_empty() {
        return 0;
    }
    let reversed: Vec<bool> = compare.iter().rev().copied().collect();
    let mut bad_end = 0usize;
    for i in 0..reversed.len() {
        bad_end = i;
        let e = (i + window_size).min(reversed.len());
        let window = &reversed[i..e];
        let divergence = window.iter().filter(|&&b| b).count() as f64 / window.len() as f64;
        if divergence < max_divergence {
            bad_end = reversed.len() - i;
            break;
        }
    }
    bad_end
}

/// Mirrors `stage_two_trimming`.
fn stage_two(
    alignment: &Alignment,
    window_size: usize,
    max_divergence: f64,
    min_len: usize,
) -> Option<Alignment> {
    if alignment.ntax() == 0 || alignment.nchar() == 0 {
        return None;
    }
    let consensus = compute_consensus(alignment);
    let mut out_rows = Vec::with_capacity(alignment.ntax());
    for row in &alignment.rows {
        let upper: Vec<u8> = row.seq.iter().map(|c| c.to_ascii_uppercase()).collect();
        let (start, end) = get_ends(&upper);
        let (seq_slice, consensus_slice): (&[u8], &[u8]) = if start < end {
            (&upper[start..end], &consensus[start..end])
        } else {
            (&[], &[])
        };
        let compare: Vec<bool> = seq_slice
            .iter()
            .zip(consensus_slice.iter())
            .map(|(a, b)| a != b)
            .collect();
        let bad_start = find_divergence_start(&compare, window_size, max_divergence);
        let bad_end = find_divergence_end(&compare, window_size, max_divergence);

        let mut orig = upper;
        let lo = (start + bad_start).min(orig.len());
        for b in orig.iter_mut().take(lo) {
            *b = b'-';
        }
        let hi = (start + bad_end).min(orig.len());
        for b in orig.iter_mut().skip(hi) {
            *b = b'-';
        }

        if !is_all_gaps(&orig) && orig.len() >= min_len {
            out_rows.push(AlignmentRow {
                id: row.id.clone(),
                seq: orig,
            });
        } else {
            return None;
        }
    }
    Some(Alignment { rows: out_rows })
}

/// Mirrors `trim_alignment(method="running", ...)`: three passes -- trim
/// block edges, trim per-row divergent edges, then re-trim block edges
/// (converting any newly-all-gap edges to `?`). Returns `None` wherever the
/// Python would leave `self.trimmed` as `None` (including on any of the
/// internal "drop this alignment" conditions the original catches with a
/// bare `except: pass`).
pub fn trim_alignment_running(
    alignment: &Alignment,
    window_size: usize,
    proportion: f64,
    threshold: f64,
    max_divergence: f64,
    min_len: usize,
) -> Option<Alignment> {
    let s1 = stage_one(
        alignment,
        window_size,
        proportion,
        threshold,
        min_len,
        false,
    )?;
    let s2 = stage_two(&s1, window_size, max_divergence, min_len)?;
    stage_one(&s2, window_size, proportion, threshold, min_len, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(rows: &[(&str, &str)]) -> Alignment {
        Alignment::from_pairs(
            rows.iter()
                .map(|(id, s)| (id.to_string(), s.to_string()))
                .collect(),
        )
    }

    #[test]
    fn python_round_matches_banker_rounding() {
        assert_eq!(python_round(2.5), 2);
        assert_eq!(python_round(3.5), 4);
        assert_eq!(python_round(3.25), 3);
    }

    #[test]
    fn rejects_invalid_trim_parameters() {
        assert_eq!(
            validate_trim_parameters(0, 0.65, 0.65, 0.4),
            Err(TrimParameterError::ZeroWindow)
        );
        assert!(validate_trim_parameters(20, f64::NAN, 0.65, 0.4).is_err());
        assert!(validate_trim_parameters(20, 0.65, 1.1, 0.4).is_err());
        assert!(validate_trim_parameters(20, 0.65, 0.65, -0.1).is_err());
        assert!(validate_trim_parameters(20, 0.65, 0.65, 0.4).is_ok());
    }

    #[test]
    fn drops_alignment_shorter_than_min_len() {
        let a = aln(&[
            ("a", "ACGTACGTAC"),
            ("b", "ACGTACGTAC"),
            ("c", "ACGTACGTAC"),
        ]);
        let result = trim_alignment_running(&a, 20, 0.65, 0.65, 0.20, 100);
        assert!(result.is_none());
    }

    #[test]
    fn consensus_picks_majority_base() {
        let a = aln(&[("a", "AAAA"), ("b", "AAAA"), ("c", "CCCC")]);
        let consensus = compute_consensus(&a);
        assert_eq!(consensus, b"AAAA".to_vec());
    }
}
