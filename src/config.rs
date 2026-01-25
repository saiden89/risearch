use crate::seed::{MismatchSpec, SeedSpec};
use crate::types::SeedPairingMode;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Matrix {
    /// Turner 1999 RNA-RNA parameters
    T99,
    /// Turner 2004 RNA-RNA parameters (default)
    T04,

    // ========================================================================
    // TODO: Placeholder matrix types from C implementation - not yet implemented
    // ========================================================================
    /// TODO: SantaLucia 1995 RNA-DNA duplex parameters
    Su95,

    /// TODO: SantaLucia 1995 RNA-DNA modified for CRISPRoff2
    Su95c2,

    /// TODO: SantaLucia 1995 RNA-DNA with mismatches as loop size 2
    Su95wk11,

    /// TODO: SantaLucia 1995 RNA-DNA without G-U wobble pairs
    Su95NoGU,

    /// TODO: SantaLucia 2004 DNA-DNA without G-T wobble pairs
    Sl04NoGU,
}

impl Default for Matrix {
    fn default() -> Self {
        Self::T04
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    /// Report predictions in detailed format (C: -p or -p1)
    Detailed,
    /// Report predictions in a simple format together with CIGAR-like string for interaction structure (C: -p2)
    Cigar,
    /// Report predictions in a simple format together with binding site (3'->5'), flanking 5'end (3'->5') and flanking 3'end (5'->3') sequences of the target (required for post-processing of CRISPR off-target predictions) (C: -p3)
    BindingSite,

    // ========================================================================
    // TODO: Placeholder output format from C implementation - not yet implemented
    // ========================================================================
    /// TODO: Minimal format - target, start, strand, and energy only (C: -p4)
    Minimal,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputCompression {
    /// No compression (default)
    None,
    /// Gzip compression
    Gzip,
    /// Zstandard compression
    Zstd,
}

/// Arguments for seed generation
#[derive(Debug, Clone)]
pub struct SeedConfig {
    /// Seed specification (length or interval)
    pub seed: SeedSpec,

    /// DEPRECATED (will be removed in a future release): disable G-U wobble pairs within the seed
    pub no_guseed: bool,

    /// Seed pairing mode (allow_wobble or strict)
    pub pairing: SeedPairingMode,

    /// Mismatch specification (legacy -m or named overrides)
    pub mismatch: MismatchSpec,

    /// Max number of mismatches allowed in the seed (preferred)
    pub mismatch_max: Option<usize>,

    /// Min consecutive matches at seed start (prefix / 5')
    pub mismatch_prefix: Option<usize>,

    /// Min consecutive matches at seed end (suffix / 3')
    pub mismatch_suffix: Option<usize>,
}

impl SeedConfig {
    pub fn apply_pairing_overrides(&mut self, explicit_pairing: bool) {
        if self.no_guseed && !explicit_pairing {
            self.pairing = SeedPairingMode::Strict;
        }
    }

    pub fn apply_mismatch_overrides(&mut self) {
        if let Some(max) = self.mismatch_max {
            self.mismatch.max_mismatches = max;
        }
        if let Some(start) = self.mismatch_prefix {
            self.mismatch.min_prefix_matches = start;
        }
        if let Some(end) = self.mismatch_suffix {
            self.mismatch.min_suffix_matches = end;
        }
    }

    pub fn has_named_mismatch(&self) -> bool {
        self.mismatch_max.is_some()
            || self.mismatch_prefix.is_some()
            || self.mismatch_suffix.is_some()
    }
}

/// Arguments for seed extension and scoring
#[derive(Debug, Clone)]
pub struct ExtendConfig {
    /// Max extension length on the seed (do DP for max this length up- and downstream of seed)
    pub max_extension: u8, // TODO: should be strictly positive

    /// Set deltaG energy threshold (in kcal/mol) to filter predictions
    pub delta_g: f64,

    /// Energy matrix for RNA-RNA duplexes
    pub matrix: Matrix,

    /// Per-nucleotide extension penalty (in kcal/mol)
    pub penalty: f64,

    /// Energy per length threshold that filters seeds
    pub seed_energy: f64,

    /// Disable maximality check (allows redundant seeds)
    pub no_max_prune: bool,

    /// Disable shadow dedup filtering (keep hits contained by better hits)
    /// Default: true (filters contained hits). Set --no-dedup-shadow for C-compatible behavior.
    pub dedup_shadow: bool,

    // ========================================================================
    // TODO: Placeholder flags from C implementation - not yet implemented
    // ========================================================================
    /// TODO: Banded search - limits the search for bulged matches.
    /// In C: `-b band, --band=band` - Integer size of bands limiting bulge search.
    /// The minimum size is 1; use seed option to avoid any bulge.
    pub band: Option<u32>,

    /// TODO: Secondary energy matrix for custom energy parameters.
    /// In C: `-y mat2, --matrix2=mat2` - Only needed for custom energy matrices.
    pub matrix2: Option<String>,

    /// TODO: Path to directory holding custom energy matrices.
    /// In C: `-M PATH, --matpath=PATH` - Directory with energy matrix files.
    pub matpath: Option<String>,

    /// TODO: Temperature scaling for energy calculations.
    /// In C: `-K T1[,T2,T3], --temperature=T0[,T1,T2]` - Temperatures in Kelvin.
    /// T0 is the target temperature; T1/T2 only needed for custom energy parameters.
    pub temperature: Option<String>,

    /// TODO: CRISPR weighting for gRNA-target interactions.
    /// In C: `-w arr, --weights=arr` - Use "CRISPR_gRNApPAM" to weight by CRISPR/Cas9 impact.
    pub weights: Option<String>,
}

/// Options that apply to the `search` subcommand
#[derive(Debug, Clone)]
pub struct SearchArgs {
    pub seed: SeedConfig,
    pub extend: ExtendConfig,
    pub output: OutputConfig,

    // ========================================================================
    // TODO: Placeholder flags from C implementation - not yet implemented
    // ========================================================================
    /// TODO: One-vs-one mode - only print results where query name matches target name.
    /// In C: `-1, --one_vs_one` - Filters results to matching query/target names.
    pub one_vs_one: bool,

    /// TODO: 3' PAM filter - report only predictions matching a 3' PAM pattern.
    /// In C: `-3 <reg>, --three_prime_match=PC` - Regex for 3' PAM forward complement.
    /// Example for cas9: `^(.cc|.uc|.cu)` matching NGG/NAG/NGA 3' PAMs.
    /// Requires output format 3 or 4 (binding_site).
    pub three_prime_match: Option<String>,

    /// TODO: 5' PAM filter - report only predictions matching a 5' PAM pattern.
    /// In C: `-5 <reg>, --five_prime_match=PC` - Regex for 5' PAM reverse complement.
    /// Example for cas12a: `^([^a]aaa)` matching TTTV 5' PAMs.
    /// Requires output format 3 or 4 (binding_site).
    pub five_prime_match: Option<String>,
}

#[derive(Debug, Clone)]
pub struct OutputConfig {
    /// Output format
    pub format: Option<OutputFormat>,

    /// Output compression codec (overrides file extension inference; gzip/gz, zstd/zst accepted)
    pub compress: Option<OutputCompression>,

    /// Output compression level (codec-specific)
    pub level: Option<i32>,
}
