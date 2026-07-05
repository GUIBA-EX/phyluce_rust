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
    let mut sections: HashMap<String, Vec<(String, String)>> = HashMap::new();
    let mut current: Option<String> = None;
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with(';') {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            let name = line[1..line.len() - 1].trim().to_string();
            sections.entry(name.clone()).or_default();
            current = Some(name);
            continue;
        }
        if let Some(section) = &current {
            let (k, v) = match line.split_once(':').or_else(|| line.split_once('=')) {
                Some((k, v)) => (k.trim().to_string(), v.trim().to_string()),
                None => (line.to_string(), String::new()),
            };
            sections.entry(section.clone()).or_default().push((k, v));
        }
    }
    sections
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
