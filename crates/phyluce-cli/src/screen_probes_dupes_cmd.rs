//! CLI wiring for `phyluce assembly screen-probes-for-dupes`, mirroring
//! `phyluce_assembly_screen_probes_for_dupes`.
//!
//! Note: the Python original (`bin/assembly/phyluce_assembly_screen_probes_for_dupes`)
//! is Python-2-only syntax (`print get_dupes(...)` without parens) and
//! raises a `SyntaxError` under Python 3 -- it cannot run at all in this
//! environment. This port re-derives the intended behavior from
//! `phyluce.helpers.get_dupes`/`get_dupe_matches` (same "does a probe's
//! lastz self-matches include any other locus" dupe-detection logic
//! already ported in `remove_duplicate_hits_cmd.rs`) and prints the dupe
//! names sorted, one per line, instead of Python's arbitrary `set` repr.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use phyluce_io::lastz::read_lastz;

fn get_name(header: &str) -> String {
    header
        .split('|')
        .next()
        .unwrap_or(header)
        .trim_start_matches('>')
        .to_string()
}

pub fn run(lastz_file: &Path) -> anyhow::Result<()> {
    let matches_list = read_lastz(lastz_file, false)?;
    let mut matches: HashMap<String, Vec<String>> = HashMap::new();
    for m in &matches_list {
        matches
            .entry(get_name(&m.name1))
            .or_default()
            .push(get_name(&m.name2));
    }

    let mut dupes: HashSet<String> = HashSet::new();
    for (k, v) in &matches {
        if v.len() > 1 {
            for i in v {
                if i != k {
                    dupes.insert(k.clone());
                    dupes.insert(i.clone());
                }
            }
        } else if &v[0] != k {
            dupes.insert(k.clone());
        }
    }

    let mut sorted: Vec<&String> = dupes.iter().collect();
    sorted.sort();
    for d in sorted {
        println!("{d}");
    }
    Ok(())
}
