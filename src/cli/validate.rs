use anyhow::Result;

use crate::cli::SearchCommand;

pub(crate) fn validate_search_command(cmd: &SearchCommand) -> Result<()> {
    cmd.opts.output.validate()?;
    Ok(())
}
