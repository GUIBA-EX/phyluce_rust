//! Shared helper for invoking the `lastz` binary directly, mirroring
//! `phyluce.lastz.Align` and `phyluce.many_lastz.lastz_params`.
//!
//! Neither `easy_lastz_cmd` nor `run_multiple_lastzs_sqlite_cmd` (the only
//! two callers) can be exercised end-to-end in this environment: `lastz`
//! itself isn't installed here. The argument construction below is a
//! byte-for-byte port of the Python CLI strings, but has never actually
//! been run against a real `lastz` binary.

use phyluce_external::ExternalCommand;

const COMMON_FLAGS: &[&str] = &[
    "--strand=both",
    "--seed=12of19",
    "--transition",
    "--nogfextend",
    "--nochain",
    "--gap=400,30",
    "--xdrop=910",
    "--ydrop=8370",
    "--hspthresh=3000",
    "--gappedthresh=3000",
    "--noentropy",
];

/// Mirrors `phyluce.lastz.Align`'s CLI construction (used by `easy-lastz`):
/// `target[multiple,nameparse=full] query[nameparse=full] ... --format=general-:...,continuity`
/// (no trailing `coverage` column). Pure argument construction, no
/// process spawned -- factored out so it's testable without `lastz`.
#[allow(clippy::too_many_arguments)]
pub fn easy_lastz_args(
    target: &str,
    query: &str,
    coverage: f64,
    identity: f64,
    output: &str,
    min_match: Option<i64>,
) -> Vec<String> {
    let mut args = vec![
        format!("{target}[multiple,nameparse=full]"),
        format!("{query}[nameparse=full]"),
    ];
    args.extend(COMMON_FLAGS.iter().map(|s| s.to_string()));
    if let Some(min_match) = min_match {
        args.push(format!("--matchcount={min_match}"));
    } else {
        args.push(format!("--coverage={coverage}"));
    }
    args.push(format!("--identity={identity}"));
    args.push(format!("--output={output}"));
    args.push(
        "--format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity"
            .to_string(),
    );
    args
}

/// Mirrors `phyluce.many_lastz.lastz_params` (used by
/// `run-multiple-lastzs-sqlite`): same flags, but with a trailing
/// `coverage` output column, and always `--coverage=` (never
/// `--matchcount=`).
pub fn many_lastz_args(
    target: &str,
    query: &str,
    coverage: f64,
    identity: f64,
    output: &str,
) -> Vec<String> {
    let mut args = vec![
        format!("{target}[multiple]"),
        format!("{query}[nameparse=full]"),
    ];
    args.extend(COMMON_FLAGS.iter().map(|s| s.to_string()));
    args.push(format!("--coverage={coverage}"));
    args.push(format!("--identity={identity}"));
    args.push(format!("--output={output}"));
    args.push(
        "--format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity,coverage"
            .to_string(),
    );
    args
}

#[allow(clippy::too_many_arguments)]
pub fn run_easy_lastz(
    lastz_bin: &str,
    target: &str,
    query: &str,
    coverage: f64,
    identity: f64,
    output: &str,
    min_match: Option<i64>,
) -> anyhow::Result<()> {
    let args = easy_lastz_args(target, query, coverage, identity, output, min_match);
    ExternalCommand::new(lastz_bin).args(args).run()?;
    Ok(())
}

pub fn run_many_lastz(
    lastz_bin: &str,
    target: &str,
    query: &str,
    coverage: f64,
    identity: f64,
    output: &str,
) -> anyhow::Result<()> {
    let args = many_lastz_args(target, query, coverage, identity, output);
    ExternalCommand::new(lastz_bin).args(args).run()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn easy_lastz_args_use_coverage_when_no_min_match() {
        let args = easy_lastz_args("t.2bit", "q.fasta", 83.0, 92.5, "out.lastz", None);
        assert_eq!(args[0], "t.2bit[multiple,nameparse=full]");
        assert_eq!(args[1], "q.fasta[nameparse=full]");
        assert!(args.contains(&"--coverage=83".to_string()));
        assert!(!args.iter().any(|a| a.starts_with("--matchcount")));
        assert!(args.last().unwrap().ends_with(",continuity"));
    }

    #[test]
    fn easy_lastz_args_use_matchcount_when_given() {
        let args = easy_lastz_args("t.2bit", "q.fasta", 83.0, 92.5, "out.lastz", Some(40));
        assert!(args.contains(&"--matchcount=40".to_string()));
        assert!(!args.iter().any(|a| a.starts_with("--coverage")));
    }

    #[test]
    fn many_lastz_args_add_trailing_coverage_column() {
        let args = many_lastz_args("t.2bit", "q.fasta", 83.0, 92.5, "out.lastz");
        assert_eq!(args[0], "t.2bit[multiple]");
        assert!(args.last().unwrap().ends_with(",continuity,coverage"));
        assert!(args.contains(&"--coverage=83".to_string()));
    }
}
