use std::fmt::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock};

use arrow_array::builder::{Float64Builder, LargeStringBuilder, StringBuilder, UInt32Builder};
use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::{ArrayRef, RecordBatch, RecordBatchIterator};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use risearch::dsm::DsmRegistry;
use risearch::{
    run_search, Energy, ExtendConfig, FilterConfig, HitSink, OutputCompression, OutputConfig,
    OutputFormat, QueryRegistry, ScoreConfig, SearchConfig, SearchHit, SeedConfig, TargetRegistry,
};

// =============================================================================
// Schema + columnar conversion
// =============================================================================

fn search_result_schema() -> &'static SchemaRef {
    static SCHEMA: OnceLock<SchemaRef> = OnceLock::new();
    SCHEMA.get_or_init(|| {
        Arc::new(Schema::new(vec![
            Field::new("query_idx", DataType::UInt32, false),
            Field::new("target_idx", DataType::UInt32, false),
            Field::new("q_start", DataType::UInt32, false),
            Field::new("q_end", DataType::UInt32, false),
            Field::new("t_start", DataType::UInt32, false),
            Field::new("t_end", DataType::UInt32, false),
            Field::new("strand", DataType::Utf8, false),
            Field::new("energy", DataType::Float64, false),
            Field::new("alignment", DataType::LargeUtf8, true),
        ]))
    })
}

#[derive(Default)]
struct HitColumns {
    query_idx: UInt32Builder,
    target_idx: UInt32Builder,
    q_start: UInt32Builder,
    q_end: UInt32Builder,
    t_start: UInt32Builder,
    t_end: UInt32Builder,
    strand: StringBuilder,
    energy: Float64Builder,
    alignment: LargeStringBuilder,
}

/// Appends each query's hits straight into the Arrow columns as that query
/// finishes, so the full hit set is never materialized at once.
#[derive(Default)]
struct ArrowSink(Mutex<HitColumns>);

impl ArrowSink {
    fn into_batch(self, schema: &SchemaRef) -> RecordBatch {
        let mut c = self.0.into_inner().unwrap();
        let columns: Vec<ArrayRef> = vec![
            Arc::new(c.query_idx.finish()),
            Arc::new(c.target_idx.finish()),
            Arc::new(c.q_start.finish()),
            Arc::new(c.q_end.finish()),
            Arc::new(c.t_start.finish()),
            Arc::new(c.t_end.finish()),
            Arc::new(c.strand.finish()),
            Arc::new(c.energy.finish()),
            Arc::new(c.alignment.finish()),
        ];
        RecordBatch::try_new(schema.clone(), columns).expect("column count and types match schema")
    }
}

impl HitSink for ArrowSink {
    fn consume(&self, _query_idx: usize, hits: Vec<SearchHit>) -> anyhow::Result<()> {
        let mut strand_buf = [0u8; 4];
        let mut c = self.0.lock().unwrap();
        for h in &hits {
            c.query_idx.append_value(h.query_idx as u32);
            c.target_idx.append_value(h.target_idx as u32);
            c.q_start.append_value(h.q_start as u32);
            c.q_end.append_value(h.q_end as u32);
            c.t_start.append_value(h.t_start as u32);
            c.t_end.append_value(h.t_end as u32);
            c.strand
                .append_value(char::from(h.strand).encode_utf8(&mut strand_buf));
            c.energy.append_value(f64::from(h.energy));
            match &h.alignment {
                // Write the fingerprint's symbols straight into the builder
                // instead of through `Alignment::fingerprint`'s String: this runs
                // under the shared lock. `append_value("")` closes the value.
                Some(a) => {
                    for &p in a.steps() {
                        let _ = c.alignment.write_char(p.symbol());
                    }
                    c.alignment.append_value("");
                }
                None => c.alignment.append_null(),
            }
        }
        Ok(())
    }
}

// =============================================================================
// PySearchResult — Arrow C Stream Interface producer
// =============================================================================

/// The result of a `search()` call.
///
/// Consumed by `pl.from_arrow(result)` via the Arrow PyCapsule Interface.
/// One-shot: the stream is consumed on the first call to `__arrow_c_stream__`.
#[pyclass(name = "SearchResult")]
struct PySearchResult {
    batch: Option<RecordBatch>,
    schema: SchemaRef,
}

impl PySearchResult {
    fn new(batch: RecordBatch, schema: SchemaRef) -> Self {
        Self {
            batch: Some(batch),
            schema,
        }
    }
}

#[pymethods]
impl PySearchResult {
    /// Arrow PyCapsule Interface producer (`__arrow_c_stream__` protocol).
    ///
    /// Called automatically by `pl.from_arrow(result)` — do not call directly.
    #[pyo3(signature = (requested_schema = None))]
    fn __arrow_c_stream__<'py>(
        &mut self,
        py: Python<'py>,
        requested_schema: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let _ = requested_schema;
        let batch = self
            .batch
            .take()
            .ok_or_else(|| pyo3::exceptions::PyIOError::new_err("Arrow stream already consumed"))?;

        let reader = RecordBatchIterator::new(std::iter::once(Ok(batch)), self.schema.clone());
        let stream = FFI_ArrowArrayStream::new(Box::new(reader));

        PyCapsule::new_with_value(py, stream, c"arrow_array_stream")
    }
}

// =============================================================================
// PyTargetRegistry
// =============================================================================

/// An mmap-backed target index for RNA-RNA interaction search.
///
/// Build once with `build_index()`, then reuse across many `search()` calls.
#[pyclass(name = "TargetRegistry")]
pub struct PyTargetRegistry(TargetRegistry);

#[pymethods]
impl PyTargetRegistry {
    /// Load a pre-built index from disk.
    ///
    /// Parameters
    /// ----------
    /// path : str | os.PathLike
    ///     Path to the `.idx` file produced by `build_index()`.
    #[staticmethod]
    fn open(path: PathBuf) -> PyResult<Self> {
        let store = TargetRegistry::open(&path)?;
        Ok(PyTargetRegistry(store))
    }

    fn __repr__(&self) -> String {
        format!("TargetRegistry(targets={})", self.0.len())
    }
}

// =============================================================================
// Module-level functions
// =============================================================================

/// Build a binary index from a FASTA file of target sequences.
///
/// Parameters
/// ----------
/// fasta : str | os.PathLike
///     Input FASTA file containing target sequences.
/// output : str | os.PathLike
///     Destination path for the binary index (`.idx`).
///
/// Notes
/// -----
/// The GIL is released during index construction so other Python threads
/// may run concurrently.
#[pyfunction]
fn build_index(py: Python<'_>, fasta: PathBuf, output: PathBuf) -> PyResult<()> {
    py.detach(|| TargetRegistry::build(&fasta, &output, None))?;
    Ok(())
}

/// Search for RNA-RNA interactions between query sequences and an indexed target set.
///
/// Parameters
/// ----------
/// query_fasta : list[str | os.PathLike]
///     One or more FASTA files containing query sequences.
/// store : TargetRegistry
///     Pre-loaded target index (from `TargetRegistry.open()`).
///
/// Returns
/// -------
/// SearchResult
///     Arrow C Stream producer; pass to `pl.from_arrow()` to obtain a DataFrame.
///
/// Notes
/// -----
/// The GIL is released during the search so other Python threads may run.
#[pyfunction]
#[pyo3(signature = (
    query_fasta,
    store,
    *,
    seed_length = None,
    seed_start = None,
    seed_end = None,
    energy_threshold = -20.0,
    mismatches = 0,
    mismatch_prefix = 1,
    mismatch_suffix = 0,
    seed_wobble = true,
    matrix = "t04",
    penalty = 3.5,
    temperature = 37,
    max_extension = 20,
    seed_energy = 0.0,
    no_max_prune = false,
    no_dedup = false,
))]
// Binding surface: each argument maps to a documented Python keyword parameter,
// so the flat signature is the public API and grouping would break it.
#[allow(clippy::too_many_arguments)]
fn search(
    py: Python<'_>,
    query_fasta: Vec<PathBuf>,
    store: &PyTargetRegistry,
    seed_length: Option<i64>,
    seed_start: Option<i64>,
    seed_end: Option<i64>,
    energy_threshold: f64,
    mismatches: usize,
    mismatch_prefix: usize,
    mismatch_suffix: usize,
    seed_wobble: bool,
    matrix: &str,
    penalty: f64,
    temperature: i32,
    max_extension: i32,
    seed_energy: f64,
    no_max_prune: bool,
    no_dedup: bool,
) -> PyResult<PySearchResult> {
    let dsm_id = DsmRegistry::parse_id(matrix)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let config = SearchConfig {
        seed: SeedConfig {
            seed_start,
            seed_end,
            seed_length,
            seed_wobble,
            max_mismatches: mismatches,
            min_prefix_matches: mismatch_prefix,
            min_suffix_matches: mismatch_suffix,
        },
        score: ScoreConfig {
            dsm_id,
            penalty: Energy::try_from(penalty).map_err(pyo3::exceptions::PyValueError::new_err)?,
            temperature,
        },
        extend: ExtendConfig {
            max_extension,
            build_alignment: true,
        },
        filter: FilterConfig {
            delta_g: Energy::try_from(energy_threshold)
                .map_err(pyo3::exceptions::PyValueError::new_err)?,
            seed_energy: Energy::try_from(seed_energy)
                .map_err(pyo3::exceptions::PyValueError::new_err)?,
            no_max_prune,
            no_dedup,
        },
        output: OutputConfig {
            format: OutputFormat::Detailed,
            compress: OutputCompression::None,
            multifile: false,
        },
    };

    let paths: Vec<&std::path::Path> = query_fasta.iter().map(|p| p.as_path()).collect();
    let queries = QueryRegistry::from_fastas(&paths, &config.seed)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let sink = ArrowSink::default();
    py.detach(|| run_search(&queries, &store.0, &config, &sink))?;
    let schema = search_result_schema().clone();
    Ok(PySearchResult::new(sink.into_batch(&schema), schema))
}

// =============================================================================
// Module registration
// =============================================================================

#[pymodule]
fn _risearch(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTargetRegistry>()?;
    m.add_function(wrap_pyfunction!(build_index, m)?)?;
    m.add_function(wrap_pyfunction!(search, m)?)?;
    Ok(())
}
