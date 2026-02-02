use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::registry::TargetRegistry;

pub(crate) fn validate_readable_file(path: &Path) -> Result<()> {
    let md = std::fs::metadata(path)
        .with_context(|| format!("Failed to access input path: {}", path.display()))?;
    if !md.is_file() {
        bail!("Input path is not a file: {}", path.display());
    }
    Ok(())
}

pub(crate) fn validate_output_path(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            let md = std::fs::metadata(parent)
                .with_context(|| format!("Output directory does not exist: {}", parent.display()))?;
            if !md.is_dir() {
                bail!("Output parent is not a directory: {}", parent.display());
            }
        }
    }
    Ok(())
}

pub fn write_index_file(index: &TargetRegistry, output_file: impl AsRef<Path>) -> Result<()> {
    let encoded = bincode::serialize(index).context("Failed to serialize index")?;
    std::fs::write(&output_file, encoded).context("Failed to write index file")?;
    Ok(())
}

pub fn load_index_file(input_file: impl AsRef<Path>) -> Result<TargetRegistry> {
    let data = std::fs::read(&input_file).context("Failed to read index file")?;
    let index: TargetRegistry =
        bincode::deserialize(&data).context("Failed to deserialize index")?;
    Ok(index)
}
