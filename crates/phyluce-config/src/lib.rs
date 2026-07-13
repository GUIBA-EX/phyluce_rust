//! phyluce-config: parse `config/phyluce.conf` and `~/.phyluce.conf`,
//! mirroring the behavior of the legacy `phyluce/pth.py`.
//!
//! Placeholder expansion:
//! - `$CONDA`     -> `CONDA_PREFIX` env var (falls back to the running binary's prefix)
//! - `$WORKFLOWS` -> `<conda-prefix>/phyluce/workflows`, matching `__default_workflow_dir__`

use std::collections::BTreeMap;
use std::env;
use std::path::{Path, PathBuf};

const EMBEDDED_DEFAULT_CONFIG: &str = include_str!("../../../config/phyluce.conf");

#[derive(Debug, thiserror::Error)]
pub enum ConfigError {
    #[error("config file not found: searched {0:?}")]
    NotFound(Vec<PathBuf>),
    #[error("no section [{0}] in config")]
    NoSection(String),
    #[error("no key '{1}' in section [{0}]")]
    NoKey(String, String),
    #[error("CONDA_PREFIX is not set; cannot expand $CONDA in config value")]
    NoCondaPrefix,
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// An ordered `[section] key:value` store, preserving on-disk order
/// (matters for `[headers]`, whose values get joined into one alternation regex).
#[derive(Debug, Default, Clone)]
pub struct Ini {
    // section -> ordered (key, value) pairs
    sections: BTreeMap<String, Vec<(String, String)>>,
    section_order: Vec<String>,
}

impl Ini {
    fn ensure_section(&mut self, name: &str) -> &mut Vec<(String, String)> {
        if !self.sections.contains_key(name) {
            self.section_order.push(name.to_string());
        }
        self.sections.entry(name.to_string()).or_default()
    }

    /// Parse phyluce's flavor of ini: `[section]` headers, `key:value` or
    /// `key=value` pairs, `#` full-line comments, blank lines ignored.
    pub fn parse(text: &str) -> Self {
        Self::parse_impl(text, false)
    }

    /// Parse the same INI format while retaining bare values as keys with an
    /// empty value, matching ConfigParser's `allow_no_value=True` mode.
    pub fn parse_allow_no_value(text: &str) -> Self {
        Self::parse_impl(text, true)
    }

    fn parse_impl(text: &str, allow_no_value: bool) -> Self {
        let mut ini = Ini::default();
        let mut current: Option<String> = None;
        for raw_line in text.lines() {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with(';') {
                continue;
            }
            if line.starts_with('[') && line.ends_with(']') {
                let name = line[1..line.len() - 1].trim().to_string();
                ini.ensure_section(&name);
                current = Some(name);
                continue;
            }
            let Some(section) = &current else { continue };
            let sep_pos = match (line.find(':'), line.find('=')) {
                (Some(colon), Some(equals)) => Some(colon.min(equals)),
                (Some(colon), None) => Some(colon),
                (None, Some(equals)) => Some(equals),
                (None, None) => None,
            };
            if let Some(pos) = sep_pos {
                let key = line[..pos].trim().to_string();
                let value = line[pos + 1..].trim().to_string();
                ini.ensure_section(section).push((key, value));
            } else if allow_no_value {
                ini.ensure_section(section)
                    .push((line.to_string(), String::new()));
            }
        }
        ini
    }

    /// Merge `other` on top of `self`: matching (section, key) pairs are
    /// overwritten by `other`'s value; new keys/sections are appended.
    pub fn merge(&mut self, other: &Ini) {
        for section in &other.section_order {
            let entries = self.ensure_section(section);
            for (k, v) in &other.sections[section] {
                if let Some(existing) = entries.iter_mut().find(|(ek, _)| ek == k) {
                    existing.1 = v.clone();
                } else {
                    entries.push((k.clone(), v.clone()));
                }
            }
        }
    }

    pub fn get(&self, section: &str, key: &str) -> Result<&str, ConfigError> {
        let entries = self
            .sections
            .get(section)
            .ok_or_else(|| ConfigError::NoSection(section.to_string()))?;
        entries
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
            .ok_or_else(|| ConfigError::NoKey(section.to_string(), key.to_string()))
    }

    /// All values (not keys) in a section, in file order -- mirrors
    /// `helpers.get_contig_header_string`'s use of `get_all_user_params`.
    pub fn all_values(&self, section: &str) -> Option<Vec<&str>> {
        self.sections
            .get(section)
            .map(|entries| entries.iter().map(|(_, v)| v.as_str()).collect())
    }

    pub fn section_names(&self) -> impl Iterator<Item = &str> {
        self.section_order.iter().map(|s| s.as_str())
    }

    /// All (key, value) pairs in a section, in file order.
    pub fn entries(&self, section: &str) -> Option<&[(String, String)]> {
        self.sections.get(section).map(|v| v.as_slice())
    }
}

pub struct PhyluceConfig {
    pub default_path: Option<PathBuf>,
    pub user_path: Option<PathBuf>,
    ini: Ini,
}

impl PhyluceConfig {
    /// Locate and load the packaged default config plus `~/.phyluce.conf`
    /// (if present), mirroring `pth.get_user_path`'s two-file `configparser.read`.
    ///
    /// The packaged config location isn't fixed yet (no Rust install layout),
    /// so resolution order is:
    /// 1. `$PHYLUCE_CONFIG` env var, if set
    /// 2. `config/phyluce.conf` relative to the current working directory
    /// 3. `<repo>/config/phyluce.conf` walking up from CWD (dev checkout)
    pub fn load() -> Result<Self, ConfigError> {
        Self::load_from(None)
    }

    /// Same as [`load`], but load only the packaged/default config -- mirrors
    /// `pth.get_user_path(..., package_only=True)`.
    pub fn load_package_only() -> Result<Self, ConfigError> {
        let default_path = Self::find_default_config();
        let text = match &default_path {
            Some(path) => std::fs::read_to_string(path)?,
            None => EMBEDDED_DEFAULT_CONFIG.to_string(),
        };
        let ini = Ini::parse(&text);
        Ok(Self {
            default_path,
            user_path: None,
            ini,
        })
    }

    fn load_from(explicit: Option<PathBuf>) -> Result<Self, ConfigError> {
        let default_path = match explicit {
            Some(path) => Some(path),
            None => Self::find_default_config(),
        };
        let text = match &default_path {
            Some(path) => std::fs::read_to_string(path)?,
            None => EMBEDDED_DEFAULT_CONFIG.to_string(),
        };
        let mut ini = Ini::parse(&text);

        let user_path = dirs_home().map(|h| h.join(".phyluce.conf"));
        if let Some(up) = &user_path {
            if up.is_file() {
                let user_ini = Ini::parse(&std::fs::read_to_string(up)?);
                ini.merge(&user_ini);
            }
        }

        Ok(Self {
            default_path,
            user_path,
            ini,
        })
    }

    fn find_default_config() -> Option<PathBuf> {
        let mut candidates = Vec::new();

        if let Ok(p) = env::var("PHYLUCE_CONFIG") {
            candidates.push(PathBuf::from(p));
        }
        candidates.push(PathBuf::from("config/phyluce.conf"));

        // walk up from CWD looking for a dev checkout's config/phyluce.conf
        if let Ok(cwd) = env::current_dir() {
            let mut dir: Option<&Path> = Some(cwd.as_path());
            while let Some(d) = dir {
                let candidate = d.join("config/phyluce.conf");
                candidates.push(candidate);
                dir = d.parent();
            }
        }

        for c in &candidates {
            if c.is_file() {
                return Some(c.clone());
            }
        }
        None
    }

    /// Mirrors `pth.get_user_path`: fetch `[program] binary:path`, expanding
    /// `$CONDA` / `$WORKFLOWS` placeholders.
    pub fn get_user_path(&self, program: &str, binary: &str) -> Result<String, ConfigError> {
        let raw = self.ini.get(program, binary)?;
        self.expand_placeholders(raw)
    }

    /// Mirrors `pth.get_user_param`.
    pub fn get_user_param(&self, section: &str, param: &str) -> Result<&str, ConfigError> {
        self.ini.get(section, param)
    }

    /// Mirrors `pth.get_all_user_params`.
    pub fn get_all_user_params(&self, section: &str) -> Option<Vec<&str>> {
        self.ini.all_values(section)
    }

    /// Mirrors `helpers.get_contig_header_string`: join all `[headers]`
    /// regex fragments with `|`.
    pub fn get_contig_header_string(&self) -> Option<String> {
        self.get_all_user_params("headers").map(|v| v.join("|"))
    }

    pub fn section_names(&self) -> impl Iterator<Item = &str> {
        self.ini.section_names()
    }

    fn expand_placeholders(&self, raw: &str) -> Result<String, ConfigError> {
        if let Some(rest) = raw.strip_prefix("$CONDA") {
            let conda = conda_prefix().ok_or(ConfigError::NoCondaPrefix)?;
            Ok(format!("{}{}", conda.display(), rest))
        } else if let Some(rest) = raw.strip_prefix("$WORKFLOWS") {
            let conda = conda_prefix().ok_or(ConfigError::NoCondaPrefix)?;
            Ok(format!(
                "{}{}",
                conda.join("phyluce/workflows").display(),
                rest
            ))
        } else {
            Ok(shellexpand_home(raw))
        }
    }
}

fn conda_prefix() -> Option<PathBuf> {
    env::var("CONDA_PREFIX")
        .ok()
        .map(PathBuf::from)
        .or_else(|| {
            let executable = env::current_exe().ok()?;
            let bin_dir = executable.parent()?;
            bin_dir.parent().map(Path::to_path_buf)
        })
}

fn dirs_home() -> Option<PathBuf> {
    env::var("HOME").ok().map(PathBuf::from)
}

fn shellexpand_home(raw: &str) -> String {
    if let Some(rest) = raw.strip_prefix('~') {
        if let Some(home) = dirs_home() {
            return format!("{}{}", home.display(), rest);
        }
    }
    raw.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_sections_in_order() {
        let text =
            "[binaries]\nlastz:$CONDA/bin/lastz\n\n[headers]\ntrinity:comp\\d+\nvelvet:node_\\d+\n";
        let ini = Ini::parse(text);
        assert_eq!(ini.get("binaries", "lastz").unwrap(), "$CONDA/bin/lastz");
        let headers = ini.all_values("headers").unwrap();
        assert_eq!(headers, vec!["comp\\d+", "node_\\d+"]);
    }

    #[test]
    fn uses_the_first_delimiter_and_can_keep_bare_values() {
        let ini = Ini::parse("[binaries]\nmafft=/opt/tools:v2/mafft\n");
        assert_eq!(ini.get("binaries", "mafft").unwrap(), "/opt/tools:v2/mafft");

        let ini = Ini::parse_allow_no_value("[loci]\nuce-1\nuce-2\n");
        assert_eq!(
            ini.entries("loci").unwrap(),
            &[
                ("uce-1".to_string(), String::new()),
                ("uce-2".to_string(), String::new())
            ]
        );
    }

    #[test]
    fn user_config_overrides_default() {
        let mut base = Ini::parse("[binaries]\nlastz:$CONDA/bin/lastz\n");
        let user = Ini::parse("[binaries]\nlastz:/custom/lastz\n");
        base.merge(&user);
        assert_eq!(base.get("binaries", "lastz").unwrap(), "/custom/lastz");
    }

    #[test]
    fn expands_conda_placeholder() {
        std::env::set_var("CONDA_PREFIX", "/opt/conda/envs/phyluce");
        let ini = Ini::parse("[binaries]\nlastz:$CONDA/bin/lastz\n");
        let cfg = PhyluceConfig {
            default_path: None,
            user_path: None,
            ini,
        };
        assert_eq!(
            cfg.get_user_path("binaries", "lastz").unwrap(),
            "/opt/conda/envs/phyluce/bin/lastz"
        );
    }

    #[test]
    fn embedded_config_contains_required_sections() {
        let ini = Ini::parse(EMBEDDED_DEFAULT_CONFIG);
        assert_eq!(ini.get("binaries", "lastz").unwrap(), "$CONDA/bin/lastz");
        assert_eq!(
            ini.get("binaries", "raxmlHPC-SSE3").unwrap(),
            "$CONDA/bin/raxmlHPC-SSE3"
        );
        assert_eq!(
            ini.get("binaries", "raxml-ng").unwrap(),
            "$CONDA/bin/raxml-ng"
        );
        assert_eq!(ini.get("spades", "cov_cutoff").unwrap(), "5");
        assert!(ini
            .all_values("headers")
            .unwrap()
            .join("|")
            .contains("NODE_"));
    }
}
