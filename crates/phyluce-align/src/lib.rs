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
}
