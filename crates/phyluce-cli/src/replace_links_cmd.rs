//! CLI wiring for `phyluce utilities replace-many-links`, mirroring
//! `phyluce_utilities_replace_many_links`.

use std::path::Path;

pub fn run(indir: &Path, oldpath: &str, newpath: &str, outdir: &Path) -> anyhow::Result<()> {
    let mut links = Vec::new();
    for entry in std::fs::read_dir(indir)? {
        let path = entry?.path();
        anyhow::ensure!(
            path.symlink_metadata()?.file_type().is_symlink(),
            "Not all paths are links: {}",
            path.display()
        );
        links.push(path);
    }

    for link in &links {
        let new_link_name = outdir.join(link.file_name().unwrap());
        let target = std::fs::read_link(link)?;
        let target_str = target.to_string_lossy();
        let new_target = target_str.replace(oldpath, newpath);
        anyhow::ensure!(
            Path::new(&new_target).is_file(),
            "The new target is not a file: {new_target}"
        );
        #[cfg(unix)]
        std::os::unix::fs::symlink(&new_target, &new_link_name)?;
        #[cfg(not(unix))]
        anyhow::bail!("symlink creation is only implemented on unix");
    }
    Ok(())
}
