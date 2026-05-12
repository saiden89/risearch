use std::path::Path;

use anyhow::{bail, Context, Result};

pub(crate) fn validate_readable_file(path: &Path) -> Result<()> {
    let md = fs_err::metadata(path)
        .with_context(|| format!("Failed to access input path: {}", path.display()))?;
    if !md.is_file() {
        bail!("Input path is not a file: {}", path.display());
    }
    Ok(())
}

pub(crate) fn validate_output_path(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            let md = fs_err::metadata(parent).with_context(|| {
                format!("Output directory does not exist: {}", parent.display())
            })?;
            if !md.is_dir() {
                bail!("Output parent is not a directory: {}", parent.display());
            }
        }
    }
    Ok(())
}
