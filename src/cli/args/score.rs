use crate::config::{self, Matrix};

/// Arguments for global scoring model
#[derive(clap::Args, Debug, Clone)]
pub struct ScoreArgs {
    /// Energy matrix for RNA-RNA duplexes
    #[arg(short = 'z', long = "matrix", value_name = "MATRIX", default_value_t = Matrix::T04, value_enum)]
    pub matrix: Matrix,

    /// Per-nucleotide penalty used by the scoring model (in kcal/mol)
    #[arg(
        short = 'd',
        long = "penalty",
        value_name = "PENALTY",
        default_value_t = 0.0
    )]
    pub penalty: f64,

    /// TODO: Secondary energy matrix for custom energy parameters.
    /// In C: `-y mat2, --matrix2=mat2` - Only needed for custom energy matrices.
    #[arg(long = "matrix2", value_name = "MATRIX2", hide = true)]
    pub matrix2: Option<String>,

    /// TODO: Path to directory holding custom energy matrices.
    /// In C: `-M PATH, --matpath=PATH` - Directory with energy matrix files.
    #[arg(long = "matpath", value_name = "PATH", hide = true)]
    pub matpath: Option<String>,

    /// TODO: Temperature scaling for energy calculations.
    /// In C: `-K T1[,T2,T3], --temperature=T0[,T1,T2]` - Temperatures in Kelvin.
    /// T0 is the target temperature; T1/T2 only needed for custom energy parameters.
    #[arg(long = "temperature", value_name = "T1[,T2,T3]", hide = true)]
    pub temperature: Option<String>,

    /// TODO: CRISPR weighting for gRNA-target interactions.
    /// In C: `-w arr, --weights=arr` - Use "CRISPR_gRNApPAM" to weight by CRISPR/Cas9 impact.
    #[arg(long = "weights", value_name = "WEIGHTS", hide = true)]
    pub weights: Option<String>,
}

impl From<ScoreArgs> for config::ScoreConfig {
    fn from(value: ScoreArgs) -> Self {
        config::ScoreConfig {
            matrix: value.matrix,
            penalty: value.penalty,
            matrix2: value.matrix2,
            matpath: value.matpath,
            temperature: value.temperature,
            weights: value.weights,
        }
    }
}
