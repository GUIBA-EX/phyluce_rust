//! Read-only parser for the UCSC `.2bit` genome format, standing in for
//! `bx.seq.twobit.TwoBitFile`. See
//! <https://genome.ucsc.edu/FAQ/FAQformat.html#format7> for the format.

use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum TwoBitError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("not a valid .2bit file (bad signature)")]
    BadSignature,
    #[error("truncated .2bit file")]
    Truncated,
    #[error("sequence {0:?} not found in .2bit file")]
    NotFound(String),
}

pub struct TwoBitFile {
    data: Vec<u8>,
    big_endian: bool,
    index: HashMap<String, u32>,
}

const SIGNATURE: u32 = 0x1A412743;

impl TwoBitFile {
    pub fn open(path: &Path) -> Result<Self, TwoBitError> {
        let data = std::fs::read(path)?;
        Self::from_bytes(data)
    }

    pub fn from_bytes(data: Vec<u8>) -> Result<Self, TwoBitError> {
        if data.len() < 16 {
            return Err(TwoBitError::Truncated);
        }
        let sig_le = u32::from_le_bytes(data[0..4].try_into().unwrap());
        let big_endian = if sig_le == SIGNATURE {
            false
        } else {
            let sig_be = u32::from_be_bytes(data[0..4].try_into().unwrap());
            if sig_be == SIGNATURE {
                true
            } else {
                return Err(TwoBitError::BadSignature);
            }
        };
        let read_u32 = |d: &[u8], off: usize| -> Result<u32, TwoBitError> {
            let b = d.get(off..off + 4).ok_or(TwoBitError::Truncated)?;
            Ok(if big_endian {
                u32::from_be_bytes(b.try_into().unwrap())
            } else {
                u32::from_le_bytes(b.try_into().unwrap())
            })
        };
        let seq_count = read_u32(&data, 8)? as usize;
        let mut offset = 16usize;
        let mut index = HashMap::with_capacity(seq_count);
        for _ in 0..seq_count {
            let name_size = *data.get(offset).ok_or(TwoBitError::Truncated)? as usize;
            offset += 1;
            let name_bytes = data
                .get(offset..offset + name_size)
                .ok_or(TwoBitError::Truncated)?;
            let name = String::from_utf8_lossy(name_bytes).to_string();
            offset += name_size;
            let seq_offset = read_u32(&data, offset)?;
            offset += 4;
            index.insert(name, seq_offset);
        }
        Ok(Self {
            data,
            big_endian,
            index,
        })
    }

    fn read_u32(&self, off: usize) -> Result<u32, TwoBitError> {
        let b = self.data.get(off..off + 4).ok_or(TwoBitError::Truncated)?;
        Ok(if self.big_endian {
            u32::from_be_bytes(b.try_into().unwrap())
        } else {
            u32::from_le_bytes(b.try_into().unwrap())
        })
    }

    pub fn names(&self) -> Vec<&str> {
        self.index.keys().map(|s| s.as_str()).collect()
    }

    pub fn contains(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }

    pub fn sequence_len(&self, name: &str) -> Result<u32, TwoBitError> {
        let off = *self
            .index
            .get(name)
            .ok_or_else(|| TwoBitError::NotFound(name.to_string()))?;
        self.read_u32(off as usize)
    }

    /// Decode the full sequence for `name`: uppercase ACGT with N-blocks
    /// written as 'N' and soft-masked regions lowercased (matching
    /// `bx.seq.twobit`'s convention).
    pub fn read_full(&self, name: &str) -> Result<Vec<u8>, TwoBitError> {
        let off = *self
            .index
            .get(name)
            .ok_or_else(|| TwoBitError::NotFound(name.to_string()))? as usize;
        let dna_size = self.read_u32(off)? as usize;
        let n_block_count = self.read_u32(off + 4)? as usize;
        let mut pos = off + 8;
        let mut n_starts = Vec::with_capacity(n_block_count);
        for i in 0..n_block_count {
            n_starts.push(self.read_u32(pos + i * 4)?);
        }
        pos += n_block_count * 4;
        let mut n_sizes = Vec::with_capacity(n_block_count);
        for i in 0..n_block_count {
            n_sizes.push(self.read_u32(pos + i * 4)?);
        }
        pos += n_block_count * 4;
        let mask_block_count = self.read_u32(pos)? as usize;
        pos += 4;
        let mut mask_starts = Vec::with_capacity(mask_block_count);
        for i in 0..mask_block_count {
            mask_starts.push(self.read_u32(pos + i * 4)?);
        }
        pos += mask_block_count * 4;
        let mut mask_sizes = Vec::with_capacity(mask_block_count);
        for i in 0..mask_block_count {
            mask_sizes.push(self.read_u32(pos + i * 4)?);
        }
        pos += mask_block_count * 4;
        pos += 4; // reserved

        let packed_len = dna_size.div_ceil(4);
        let packed = self
            .data
            .get(pos..pos + packed_len)
            .ok_or(TwoBitError::Truncated)?;

        const BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];
        let mut seq = Vec::with_capacity(dna_size);
        for i in 0..dna_size {
            let byte = packed[i / 4];
            let shift = 6 - 2 * (i % 4);
            let code = (byte >> shift) & 0x3;
            seq.push(BASES[code as usize]);
        }

        for (&s, &l) in mask_starts.iter().zip(mask_sizes.iter()) {
            let (s, l) = (s as usize, l as usize);
            for base in seq.iter_mut().skip(s).take(l) {
                *base = base.to_ascii_lowercase();
            }
        }
        for (&s, &l) in n_starts.iter().zip(n_sizes.iter()) {
            let (s, l) = (s as usize, l as usize);
            for base in seq.iter_mut().skip(s).take(l) {
                *base = b'N';
            }
        }
        Ok(seq)
    }

    /// Half-open `[start, end)` slice, clamped to the sequence bounds.
    pub fn read_slice(&self, name: &str, start: i64, end: i64) -> Result<Vec<u8>, TwoBitError> {
        let full = self.read_full(name)?;
        let len = full.len() as i64;
        let start = start.clamp(0, len) as usize;
        let end = end.clamp(0, len) as usize;
        if start >= end {
            return Ok(Vec::new());
        }
        Ok(full[start..end].to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-encode a minimal single-sequence .2bit file: "ACGTacgtNNAC"
    /// (uppercase ACGT, lowercase acgt soft-masked, an NN block, then AC).
    fn build_fixture() -> Vec<u8> {
        let name = b"chr1";
        let bases = "ACGTACGTAAACAC"; // 14 bases, mask [4,8), N block [8,10)
        let base_to_code = |b: u8| -> u8 {
            match b {
                b'T' => 0,
                b'C' => 1,
                b'A' => 2,
                b'G' => 3,
                _ => 0,
            }
        };
        let mut packed = Vec::new();
        let bytes = bases.as_bytes();
        for chunk in bytes.chunks(4) {
            let mut byte = 0u8;
            for (i, &b) in chunk.iter().enumerate() {
                byte |= base_to_code(b) << (6 - 2 * i);
            }
            packed.push(byte);
        }

        let mut buf = Vec::new();
        buf.extend_from_slice(&SIGNATURE.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes()); // version
        buf.extend_from_slice(&1u32.to_le_bytes()); // seqCount
        buf.extend_from_slice(&0u32.to_le_bytes()); // reserved
        buf.push(name.len() as u8);
        buf.extend_from_slice(name);
        let seq_offset_pos = buf.len();
        buf.extend_from_slice(&0u32.to_le_bytes()); // placeholder, patched below

        let seq_offset = buf.len() as u32;
        buf[seq_offset_pos..seq_offset_pos + 4].copy_from_slice(&seq_offset.to_le_bytes());

        buf.extend_from_slice(&(bytes.len() as u32).to_le_bytes()); // dnaSize
        buf.extend_from_slice(&1u32.to_le_bytes()); // nBlockCount
        buf.extend_from_slice(&8u32.to_le_bytes()); // nBlockStarts[0]
        buf.extend_from_slice(&2u32.to_le_bytes()); // nBlockSizes[0]
        buf.extend_from_slice(&1u32.to_le_bytes()); // maskBlockCount
        buf.extend_from_slice(&4u32.to_le_bytes()); // maskBlockStarts[0]
        buf.extend_from_slice(&4u32.to_le_bytes()); // maskBlockSizes[0]
        buf.extend_from_slice(&0u32.to_le_bytes()); // reserved
        buf.extend_from_slice(&packed);
        buf
    }

    #[test]
    fn decodes_bases_mask_and_n_blocks() {
        let tb = TwoBitFile::from_bytes(build_fixture()).unwrap();
        assert_eq!(tb.names(), vec!["chr1"]);
        assert_eq!(tb.sequence_len("chr1").unwrap(), 14);
        let full = tb.read_full("chr1").unwrap();
        assert_eq!(String::from_utf8(full).unwrap(), "ACGTacgtNNACAC");
    }

    #[test]
    fn slices_are_clamped_to_bounds() {
        let tb = TwoBitFile::from_bytes(build_fixture()).unwrap();
        assert_eq!(tb.read_slice("chr1", 0, 4).unwrap(), b"ACGT");
        assert_eq!(tb.read_slice("chr1", 10, 1000).unwrap(), b"ACAC");
        assert_eq!(tb.read_slice("chr1", -5, 3).unwrap(), b"ACG");
    }
}
