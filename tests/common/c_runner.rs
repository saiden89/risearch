//! C binary runner with type-state pattern.
//!
//! Provides `CRunner<S>` for running the legacy C risearch2 binary
//! with compile-time guarantees about index state.

use flate2::read::GzDecoder;
use log::debug;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

// =============================================================================
// TYPE-STATE MARKERS
// =============================================================================

/// Marker for an unindexed C runner.
pub struct NoIndex;

/// Marker for an indexed C runner.
pub struct Indexed {
    index_path: PathBuf,
}

// =============================================================================
// C RUNNER
// =============================================================================

/// Runner for the legacy C risearch2 binary.
///
/// Uses type-state pattern to ensure you can only search after indexing.
pub struct CRunner<S> {
    bin_path: PathBuf,
    state: S,
}

impl CRunner<NoIndex> {
    /// Create a new C runner pointing to the binary.
    pub fn new(root: &Path) -> Self {
        let bin_path = c_binary_path(root);
        Self {
            bin_path,
            state: NoIndex,
        }
    }

    /// Create an index from the target file.
    /// Consumes self and returns an indexed runner.
    pub fn create_index(self, target: &Path, index_out: &Path) -> CRunner<Indexed> {
        let output = Command::new(&self.bin_path)
            .arg("-c")
            .arg(target)
            .arg("-o")
            .arg(index_out)
            .output()
            .expect("create c index");

        if !output.status.success() {
            panic!(
                "C index creation failed: {:?}\nstderr: {}",
                output.status,
                String::from_utf8_lossy(&output.stderr)
            );
        }

        CRunner {
            bin_path: self.bin_path,
            state: Indexed {
                index_path: index_out.to_path_buf(),
            },
        }
    }
}

impl CRunner<Indexed> {
    /// Run a search against the indexed target.
    pub fn search(&self, query: &Path, args: &[&str]) -> String {
        let tmpdir = tempfile::tempdir().expect("tempdir");

        let mut final_args = vec![
            "-q",
            query.to_str().unwrap(),
            "-i",
            self.state.index_path.to_str().unwrap(),
        ];
        final_args.extend_from_slice(args);

        let out = Command::new(&self.bin_path)
            .current_dir(tmpdir.path())
            .args(&final_args)
            .output()
            .expect("run legacy C risearch2");

        if !out.status.success() {
            panic!(
                "Legacy C risearch2 failed: status={:?}\nstdout=\n{}\nstderr=\n{}\n",
                out.status,
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr)
            );
        }

        // Print C binary's debug output (stderr) for parity debugging
        let c_stderr = String::from_utf8_lossy(&out.stderr);
        if !c_stderr.is_empty() {
            for line in c_stderr.lines() {
                log::trace!("[C] {}", line);
            }
        }

        let mut out_files: Vec<PathBuf> = fs::read_dir(tmpdir.path())
            .expect("read legacy C output dir")
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.starts_with("risearch_") && n.ends_with(".out.gz"))
            })
            .collect();

        if out_files.is_empty() {
            panic!(
                "Legacy C risearch2 produced no output files!\nArgs: {:?}\nStderr:\n{}\nStdout:\n{}\n",
                final_args,
                String::from_utf8_lossy(&out.stderr),
                String::from_utf8_lossy(&out.stdout)
            );
        }

        out_files.sort();

        let mut combined = String::new();
        for path in out_files {
            let bytes = fs::read(&path).expect("read legacy output file");
            let mut decoder = GzDecoder::new(&bytes[..]);
            let mut s = String::new();
            decoder
                .read_to_string(&mut s)
                .expect("decode legacy .out.gz as utf8");
            combined.push_str(&s);
            if !combined.ends_with('\n') {
                combined.push('\n');
            }
        }

        if combined.trim().is_empty() {
            panic!(
                "Legacy C risearch2 produced empty output!\nArgs: {:?}\nStderr:\n{}\nStdout:\n{}\n",
                final_args,
                String::from_utf8_lossy(&out.stderr),
                String::from_utf8_lossy(&out.stdout)
            );
        }

        if combined.is_empty() {
            String::from_utf8_lossy(&out.stdout).to_string()
        } else {
            combined
        }
    }

    /// Get the index path.
    #[allow(dead_code)] // Useful API for debugging
    pub fn index_path(&self) -> &Path {
        &self.state.index_path
    }
}

// =============================================================================
// BINARY PATH HELPER
// =============================================================================

/// Get path to C risearch2 binary.
/// Prefers debug binary (risearch2.dbg.x) for seed boundary markers, falls back to release.
pub fn c_binary_path(root: &Path) -> PathBuf {
    let debug_bin = root.join("legacy_c/RIsearch2/bin/risearch2.dbg.x");
    let release_bin = root.join("legacy_c/RIsearch2/bin/risearch2.x");

    if debug_bin.exists() {
        debug!("[PARITY] Using debug C binary: {}", debug_bin.display());
        debug_bin
    } else {
        debug!("[PARITY] Using release C binary: {}", release_bin.display());
        release_bin
    }
}

// =============================================================================
// UNIT TESTS
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c_binary_path_structure() {
        let root = PathBuf::from("/fake/root");
        let debug_path = root.join("legacy_c/RIsearch2/bin/risearch2.dbg.x");
        let release_path = root.join("legacy_c/RIsearch2/bin/risearch2.x");

        assert!(debug_path.to_str().unwrap().contains("dbg"));
        assert!(!release_path.to_str().unwrap().contains("dbg"));
    }

    // Type-state compile-time test:
    // The following would NOT compile if uncommented:
    // ```
    // fn test_type_state_prevents_search_without_index() {
    //     let runner = CRunner::<NoIndex>::new(&PathBuf::from("/fake"));
    //     runner.search(...); // ERROR: no method `search` on `CRunner<NoIndex>`
    // }
    // ```
}
