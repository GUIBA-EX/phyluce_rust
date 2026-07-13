//! Read-only parser for the UCSC `.2bit` genome format, standing in for
//! `bx.seq.twobit.TwoBitFile`. See
//! <https://genome.ucsc.edu/FAQ/FAQformat.html#format7> for the format.

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Mutex;

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
    #[error("invalid .2bit file: {0}")]
    Invalid(String),
}

pub struct TwoBitFile {
    source: Source,
    big_endian: bool,
    index: HashMap<String, u64>,
}

enum Source {
    Bytes(Vec<u8>),
    File { file: Mutex<File>, len: u64 },
}

const SIGNATURE: u32 = 0x1A412743;

impl TwoBitFile {
    pub fn open(path: &Path) -> Result<Self, TwoBitError> {
        let file = File::open(path)?;
        let len = file.metadata()?.len();
        Self::from_source(Source::File {
            file: Mutex::new(file),
            len,
        })
    }

    pub fn from_bytes(data: Vec<u8>) -> Result<Self, TwoBitError> {
        Self::from_source(Source::Bytes(data))
    }

    fn from_source(source: Source) -> Result<Self, TwoBitError> {
        if source.len() < 16 {
            return Err(TwoBitError::Truncated);
        }
        let signature: [u8; 4] = source
            .read_at(0, 4)?
            .try_into()
            .map_err(|_| TwoBitError::Truncated)?;
        let sig_le = u32::from_le_bytes(signature);
        let big_endian = if sig_le == SIGNATURE {
            false
        } else {
            let sig_be = u32::from_be_bytes(signature);
            if sig_be == SIGNATURE {
                true
            } else {
                return Err(TwoBitError::BadSignature);
            }
        };
        let read_u32 = |off| read_u32_from(&source, big_endian, off);
        let seq_count = read_u32(8)? as usize;
        let max_entries = ((source.len().saturating_sub(16)) / 5) as usize;
        if seq_count > max_entries {
            return Err(TwoBitError::Invalid(
                "sequence index exceeds file size".to_string(),
            ));
        }
        let mut offset = 16u64;
        let mut index = HashMap::new();
        for _ in 0..seq_count {
            let name_size = source.read_at(offset, 1)?[0] as usize;
            offset += 1;
            let name_bytes = source.read_at(offset, name_size)?;
            let name = String::from_utf8_lossy(&name_bytes).to_string();
            offset = checked_add(offset, name_size as u64)?;
            let seq_offset = read_u32(offset)?;
            offset = checked_add(offset, 4)?;
            index.insert(name, seq_offset as u64);
        }
        Ok(Self {
            source,
            big_endian,
            index,
        })
    }

    fn read_u32(&self, off: u64) -> Result<u32, TwoBitError> {
        read_u32_from(&self.source, self.big_endian, off)
    }

    fn sequence_offset(&self, name: &str) -> Result<u64, TwoBitError> {
        self.index
            .get(name)
            .copied()
            .ok_or_else(|| TwoBitError::NotFound(name.to_string()))
    }

    fn sequence_layout(&self, name: &str) -> Result<SequenceLayout, TwoBitError> {
        let offset = self.sequence_offset(name)?;
        let dna_size = self.read_u32(offset)? as u64;
        let n_block_count = self.read_u32(checked_add(offset, 4)?)? as usize;
        let n_starts_offset = checked_add(offset, 8)?;
        let n_sizes_offset = checked_add(n_starts_offset, bytes_for_u32s(n_block_count)?)?;
        let mask_count_offset = checked_add(n_sizes_offset, bytes_for_u32s(n_block_count)?)?;
        let mask_block_count = self.read_u32(mask_count_offset)? as usize;
        let mask_starts_offset = checked_add(mask_count_offset, 4)?;
        let mask_sizes_offset = checked_add(mask_starts_offset, bytes_for_u32s(mask_block_count)?)?;
        let reserved_offset = checked_add(mask_sizes_offset, bytes_for_u32s(mask_block_count)?)?;
        let packed_offset = checked_add(reserved_offset, 4)?;
        let packed_len = dna_size.div_ceil(4);
        self.source.ensure_range(
            packed_offset,
            usize::try_from(packed_len).map_err(|_| {
                TwoBitError::Invalid("packed sequence length exceeds address space".to_string())
            })?,
        )?;
        Ok(SequenceLayout {
            dna_size,
            n_starts_offset,
            n_sizes_offset,
            n_block_count,
            mask_starts_offset,
            mask_sizes_offset,
            mask_block_count,
            packed_offset,
        })
    }

    fn read_blocks(
        &self,
        starts_offset: u64,
        sizes_offset: u64,
        count: usize,
        dna_size: u64,
    ) -> Result<Vec<(u64, u64)>, TwoBitError> {
        let mut blocks = Vec::with_capacity(count);
        for i in 0..count {
            let delta = checked_mul(i as u64, 4)?;
            let start = self.read_u32(checked_add(starts_offset, delta)?)? as u64;
            let len = self.read_u32(checked_add(sizes_offset, delta)?)? as u64;
            if start > dna_size || len > dna_size.saturating_sub(start) {
                return Err(TwoBitError::Invalid(
                    "block exceeds sequence length".to_string(),
                ));
            }
            blocks.push((start, len));
        }
        Ok(blocks)
    }

    fn apply_blocks(seq: &mut [u8], slice_start: u64, blocks: &[(u64, u64)], value: u8) {
        let slice_end = slice_start + seq.len() as u64;
        for &(start, len) in blocks {
            let end = start + len;
            let overlap_start = start.max(slice_start);
            let overlap_end = end.min(slice_end);
            if overlap_start < overlap_end {
                for base in &mut seq
                    [(overlap_start - slice_start) as usize..(overlap_end - slice_start) as usize]
                {
                    *base = value;
                }
            }
        }
    }

    fn apply_mask_blocks(seq: &mut [u8], slice_start: u64, blocks: &[(u64, u64)]) {
        let slice_end = slice_start + seq.len() as u64;
        for &(start, len) in blocks {
            let end = start + len;
            let overlap_start = start.max(slice_start);
            let overlap_end = end.min(slice_end);
            if overlap_start < overlap_end {
                for base in &mut seq
                    [(overlap_start - slice_start) as usize..(overlap_end - slice_start) as usize]
                {
                    *base = base.to_ascii_lowercase();
                }
            }
        }
    }

    pub fn names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.index.keys().map(String::as_str).collect();
        names.sort_unstable();
        names
    }

    pub fn contains(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }

    pub fn sequence_len(&self, name: &str) -> Result<u32, TwoBitError> {
        self.read_u32(self.sequence_offset(name)?)
    }

    /// Decode the full sequence for `name`: uppercase ACGT with N-blocks
    /// written as 'N' and soft-masked regions lowercased (matching
    /// `bx.seq.twobit`'s convention).
    pub fn read_full(&self, name: &str) -> Result<Vec<u8>, TwoBitError> {
        let dna_size = self.sequence_len(name)? as i64;
        self.read_slice(name, 0, dna_size)
    }

    /// Half-open `[start, end)` slice, clamped to the sequence bounds. Only
    /// the requested packed bases are read and decoded; it does not allocate
    /// the full chromosome for every interval.
    pub fn read_slice(&self, name: &str, start: i64, end: i64) -> Result<Vec<u8>, TwoBitError> {
        let layout = self.sequence_layout(name)?;
        let sequence_len = layout.dna_size as i64;
        let start = start.clamp(0, sequence_len) as u64;
        let end = end.clamp(0, sequence_len) as u64;
        if start >= end {
            return Ok(Vec::new());
        }

        let packed_start = start / 4;
        let packed_end = end.div_ceil(4);
        let packed_len = usize::try_from(packed_end - packed_start)
            .map_err(|_| TwoBitError::Invalid("slice length exceeds address space".to_string()))?;
        let packed = self
            .source
            .read_at(checked_add(layout.packed_offset, packed_start)?, packed_len)?;

        const BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];
        let mut seq = Vec::with_capacity((end - start) as usize);
        for position in start..end {
            let byte = packed[(position / 4 - packed_start) as usize];
            let shift = 6 - 2 * (position % 4);
            seq.push(BASES[((byte >> shift) & 0x3) as usize]);
        }

        let masks = self.read_blocks(
            layout.mask_starts_offset,
            layout.mask_sizes_offset,
            layout.mask_block_count,
            layout.dna_size,
        )?;
        Self::apply_mask_blocks(&mut seq, start, &masks);
        let n_blocks = self.read_blocks(
            layout.n_starts_offset,
            layout.n_sizes_offset,
            layout.n_block_count,
            layout.dna_size,
        )?;
        Self::apply_blocks(&mut seq, start, &n_blocks, b'N');
        Ok(seq)
    }
}

struct SequenceLayout {
    dna_size: u64,
    n_starts_offset: u64,
    n_sizes_offset: u64,
    n_block_count: usize,
    mask_starts_offset: u64,
    mask_sizes_offset: u64,
    mask_block_count: usize,
    packed_offset: u64,
}

impl Source {
    fn len(&self) -> u64 {
        match self {
            Self::Bytes(data) => data.len() as u64,
            Self::File { len, .. } => *len,
        }
    }

    fn read_at(&self, offset: u64, len: usize) -> Result<Vec<u8>, TwoBitError> {
        let end = self.ensure_range(offset, len)?;
        match self {
            Self::Bytes(data) => Ok(data[offset as usize..end as usize].to_vec()),
            Self::File { file, .. } => {
                let mut file = file.lock().map_err(|_| {
                    TwoBitError::Invalid("2bit file handle lock was poisoned".to_string())
                })?;
                file.seek(SeekFrom::Start(offset))?;
                let mut buffer = vec![0; len];
                file.read_exact(&mut buffer)?;
                Ok(buffer)
            }
        }
    }

    fn ensure_range(&self, offset: u64, len: usize) -> Result<u64, TwoBitError> {
        let end = checked_add(offset, len as u64)?;
        if end > self.len() {
            return Err(TwoBitError::Truncated);
        }
        Ok(end)
    }
}

fn read_u32_from(source: &Source, big_endian: bool, offset: u64) -> Result<u32, TwoBitError> {
    let bytes = source.read_at(offset, 4)?;
    let bytes: [u8; 4] = bytes.try_into().map_err(|_| TwoBitError::Truncated)?;
    Ok(if big_endian {
        u32::from_be_bytes(bytes)
    } else {
        u32::from_le_bytes(bytes)
    })
}

fn checked_add(a: u64, b: u64) -> Result<u64, TwoBitError> {
    a.checked_add(b)
        .ok_or_else(|| TwoBitError::Invalid("offset overflow".to_string()))
}

fn checked_mul(a: u64, b: u64) -> Result<u64, TwoBitError> {
    a.checked_mul(b)
        .ok_or_else(|| TwoBitError::Invalid("offset overflow".to_string()))
}

fn bytes_for_u32s(count: usize) -> Result<u64, TwoBitError> {
    checked_mul(count as u64, 4)
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

    #[test]
    fn file_backed_slices_match_in_memory_slices() {
        let path = std::env::temp_dir().join(format!("phyluce-twobit-{}.2bit", std::process::id()));
        std::fs::write(&path, build_fixture()).unwrap();
        let tb = TwoBitFile::open(&path).unwrap();
        assert_eq!(tb.read_slice("chr1", 3, 11).unwrap(), b"TacgtNNA");
        std::fs::remove_file(path).unwrap();
    }
}
