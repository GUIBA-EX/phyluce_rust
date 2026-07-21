//! CLI wiring for `phyluce probe easy-lastz`, mirroring
//! `phyluce_probe_easy_lastz`.

use std::path::Path;

use phyluce_config::PhyluceConfig;

pub fn run(
    target: &Path,
    query: &Path,
    output: &Path,
    coverage: f64,
    identity: f64,
    min_match: Option<i64>,
) -> anyhow::Result<()> {
    let cfg = PhyluceConfig::load()?;
    let lastz_bin = cfg.get_user_path("binaries", "lastz")?;
    crate::lastz_align::run_easy_lastz(
        &lastz_bin,
        &target.to_string_lossy(),
        &query.to_string_lossy(),
        coverage,
        identity,
        &output.to_string_lossy(),
        min_match,
    )
}
