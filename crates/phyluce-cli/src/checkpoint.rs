//! Content-checked checkpoints for expensive external operations.
use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::Context;

pub fn file_digest(path: &Path) -> anyhow::Result<blake3::Hash> {
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("opening {} for checkpoint validation", path.display()))?;
    let mut hasher = blake3::Hasher::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let count = file
            .read(&mut buffer)
            .with_context(|| format!("hashing {}", path.display()))?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }
    Ok(hasher.finalize())
}

/// Resolve PATH once, so the fingerprint and invocation use the same binary.
pub fn resolve_binary(binary: &str) -> anyhow::Result<PathBuf> {
    let path = Path::new(binary);
    if path.components().count() > 1 || path.is_absolute() {
        return std::fs::canonicalize(path)
            .with_context(|| format!("resolving external binary {binary}"));
    }
    for directory in std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default()) {
        let candidate = directory.join(binary);
        if !candidate.is_file() {
            continue;
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            if candidate.metadata()?.permissions().mode() & 0o111 == 0 {
                continue;
            }
        }
        return Ok(std::fs::canonicalize(candidate)?);
    }
    anyhow::bail!("external binary {binary:?} not found on PATH")
}

pub fn lastz_fingerprint(
    contigs: &Path,
    probes: &Path,
    binary: &Path,
    coverage: u32,
    identity: u32,
) -> anyhow::Result<String> {
    let mut hash = blake3::Hasher::new();
    // Bump this operation version when LASTZ flags/semantics change.
    hash.update(b"phyluce-match-contigs-lastz-v1\0");
    hash.update(env!("CARGO_PKG_VERSION").as_bytes());
    hash.update(&coverage.to_le_bytes());
    hash.update(&identity.to_le_bytes());
    for path in [contigs, probes, binary] {
        let canonical = std::fs::canonicalize(path)?;
        let name = canonical.as_os_str().as_encoded_bytes();
        hash.update(&(name.len() as u64).to_le_bytes());
        hash.update(name);
        hash.update(file_digest(path)?.as_bytes());
    }
    Ok(hash.finalize().to_hex().to_string())
}

fn manifest_path(output: &Path) -> PathBuf {
    let mut path = output.as_os_str().to_owned();
    path.push(".complete");
    PathBuf::from(path)
}

pub fn reusable(output: &Path, fingerprint: &str) -> anyhow::Result<bool> {
    if !output.is_file() {
        return Ok(false);
    }
    let manifest = match std::fs::read_to_string(manifest_path(output)) {
        Ok(text) => text,
        Err(error)
            if matches!(
                error.kind(),
                std::io::ErrorKind::NotFound | std::io::ErrorKind::InvalidData
            ) =>
        {
            return Ok(false)
        }
        Err(error) => return Err(error.into()),
    };
    let mut lines = manifest.lines();
    if lines.next() != Some(fingerprint) {
        return Ok(false);
    }
    let expected = file_digest(output)?.to_hex().to_string();
    Ok(lines.next() == Some(expected.as_str()) && lines.next().is_none())
}

pub fn invalidate(output: &Path) -> anyhow::Result<()> {
    crate::output_path::remove_stale_file(&manifest_path(output))?;
    Ok(())
}

pub fn complete(output: &Path, fingerprint: &str) -> anyhow::Result<()> {
    let contents = format!("{fingerprint}\n{}\n", file_digest(output)?.to_hex());
    crate::output_path::write_atomic(&manifest_path(output), contents)
}
