//! CLI wiring for `phyluce align get-ry-recoded-alignments`, mirroring
//! `phyluce_align_get_ry_recoded_alignments`.

use std::path::Path;

use phyluce_align::nexus::format_nexus;
use phyluce_align::{Alignment, AlignmentRow};

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

/// Mirrors the RY-recode translation table:
/// `AGRCTYSWKMBDHV -> RRRYYY????????` (and lowercase).
fn ry_translate(c: u8) -> u8 {
    match c {
        b'A' | b'G' | b'R' => b'R',
        b'C' | b'T' | b'Y' => b'Y',
        b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => b'?',
        b'a' | b'g' | b'r' => b'r',
        b'c' | b't' | b'y' => b'y',
        b's' | b'w' | b'k' | b'm' | b'b' | b'd' | b'h' | b'v' => b'?',
        other => other,
    }
}

/// Mirrors the binary-recode translation table:
/// `AGRCTYSWKMBDHVN -> 111000?????????` (and lowercase).
fn binary_translate(c: u8) -> u8 {
    match c {
        b'A' | b'G' | b'R' => b'1',
        b'C' | b'T' | b'Y' => b'0',
        b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N' => b'?',
        b'a' | b'g' | b'r' => b'1',
        b'c' | b't' | b'y' => b'0',
        b's' | b'w' | b'k' | b'm' | b'b' | b'd' | b'h' | b'v' | b'n' => b'?',
        other => other,
    }
}

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    binary: bool,
) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)?;
    let translate = if binary {
        binary_translate
    } else {
        ry_translate
    };

    for file in &files {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        let alignment = load_alignment(file, input_format)?;
        let recoded = Alignment {
            rows: alignment
                .rows
                .into_iter()
                .map(|r| AlignmentRow {
                    id: r.id,
                    seq: r.seq.iter().map(|&c| translate(c)).collect(),
                })
                .collect(),
        };
        let out_path = output_dir.join(format!("{stem}.nexus"));
        std::fs::write(out_path, format_nexus(&recoded))?;
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
