//! phyluce-align: alignment representation, the native phyluce 3-stage
//! edge-trimming algorithm, and format writers (mirrors
//! `phyluce/generic_align.py` + `Bio.AlignIO`'s nexus writer).

pub mod charset;
pub mod concat;
pub mod mafft;
pub mod nexus;
pub mod sites;
pub mod summary;
pub mod trim;

/// One row of an alignment: a sequence id and its (gapped) characters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentRow {
    pub id: String,
    pub seq: Vec<u8>,
}

/// A multiple sequence alignment: every row must be the same length.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Alignment {
    pub rows: Vec<AlignmentRow>,
}

#[derive(Debug, thiserror::Error)]
pub enum AlignmentError {
    #[error("alignment contains an empty taxon identifier")]
    EmptyIdentifier,
    #[error("alignment contains duplicate taxon {0:?}")]
    DuplicateTaxon(String),
    #[error("taxon {taxon:?} has {actual} characters; expected {expected}")]
    UnequalLength {
        taxon: String,
        expected: usize,
        actual: usize,
    },
}

impl Alignment {
    pub fn from_pairs(pairs: Vec<(String, String)>) -> Self {
        Alignment {
            rows: pairs
                .into_iter()
                .map(|(id, seq)| AlignmentRow {
                    id,
                    seq: seq.into_bytes(),
                })
                .collect(),
        }
    }

    pub fn ntax(&self) -> usize {
        self.rows.len()
    }

    pub fn nchar(&self) -> usize {
        self.rows.first().map(|r| r.seq.len()).unwrap_or(0)
    }

    pub fn validate(&self) -> Result<(), AlignmentError> {
        let expected = self.nchar();
        let mut taxa = std::collections::HashSet::new();
        for row in &self.rows {
            if row.id.is_empty() {
                return Err(AlignmentError::EmptyIdentifier);
            }
            if !taxa.insert(row.id.as_str()) {
                return Err(AlignmentError::DuplicateTaxon(row.id.clone()));
            }
            if row.seq.len() != expected {
                return Err(AlignmentError::UnequalLength {
                    taxon: row.id.clone(),
                    expected,
                    actual: row.seq.len(),
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod alignment_tests {
    use super::*;

    #[test]
    fn validates_row_lengths_and_unique_taxa() {
        let unequal = Alignment::from_pairs(vec![
            ("a".to_string(), "AAAA".to_string()),
            ("b".to_string(), "AAA".to_string()),
        ]);
        assert!(matches!(
            unequal.validate(),
            Err(AlignmentError::UnequalLength { .. })
        ));
        let duplicate = Alignment::from_pairs(vec![
            ("a".to_string(), "AAAA".to_string()),
            ("a".to_string(), "AAAA".to_string()),
        ]);
        assert!(matches!(
            duplicate.validate(),
            Err(AlignmentError::DuplicateTaxon(_))
        ));
    }
}
