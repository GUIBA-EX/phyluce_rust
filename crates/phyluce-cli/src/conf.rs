//! Small shared INI-style config helpers for the `utilities` commands,
//! which use plain `key:value` (or `key=value`) conf files (as opposed to
//! the `allow_no_value` bare-item lists used by the assembly/align
//! commands).

use std::collections::HashMap;

/// Parse a `[section] key:value` conf file into section -> ordered
/// (key, value) pairs. Mirrors `allow_no_value=True`: a line with no
/// `:`/`=` delimiter is kept as a bare entry with an empty value (used by
/// several commands' plain item-list sections, e.g. `[set1]\nlocusA\n`).
pub fn parse_ini(text: &str) -> HashMap<String, Vec<(String, String)>> {
    let ini = phyluce_config::Ini::parse_allow_no_value(text);
    ini.section_names()
        .map(|section| {
            (
                section.to_string(),
                ini.entries(section).unwrap_or_default().to_vec(),
            )
        })
        .collect()
}

/// Mirrors `dict((name, dirs.split(",")) for name, dirs in
/// conf.items(section))`: each value is a comma-separated list.
pub fn read_ini_values(text: &str, section: &str) -> anyhow::Result<HashMap<String, Vec<String>>> {
    let sections = parse_ini(text);
    let entries = sections
        .get(section)
        .ok_or_else(|| anyhow::anyhow!("no [{section}] section in config"))?;
    Ok(entries
        .iter()
        .map(|(k, v)| {
            (
                k.clone(),
                v.split(',').map(|s| s.trim().to_string()).collect(),
            )
        })
        .collect())
}
