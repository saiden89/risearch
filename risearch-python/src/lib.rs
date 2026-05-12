use std::ffi::CString;
use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::{
    ArrayRef, Float64Array, RecordBatch, RecordBatchIterator, StringArray, UInt32Array,
};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use risearch::dsm::DsmRegistry;
use risearch::{
    run_search_in_memory, Energy, ExtendConfig,
    FilterConfig, OutputCompression, OutputConfig, OutputFormat, QueryRegistry,
    ScoreConfig, SearchConfig, SearchHit, SeedConfig, TargetStore,
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
            Field::new("alignment", DataType::Utf8, true),
        ]))
    })
}

fn hits_to_record_batch(hits: Vec<SearchHit>, schema: &SchemaRef) -> RecordBatch {
    let n = hits.len();
    let mut query_idx = Vec::with_capacity(n);
    let mut target_idx = Vec::with_capacity(n);
    let mut q_start = Vec::with_capacity(n);
    let mut q_end = Vec::with_capacity(n);
    let mut t_start = Vec::with_capacity(n);
    let mut t_end = Vec::with_capacity(n);
    let mut strand: Vec<String> = Vec::with_capacity(n);
    let mut energy = Vec::with_capacity(n);
    let mut alignment: Vec<Option<String>> = Vec::with_capacity(n);

    for h in hits {
        query_idx.push(h.query_idx as u32);
        target_idx.push(h.target_idx as u32);
        q_start.push(h.q_start as u32);
        q_end.push(h.q_end as u32);
        t_start.push(h.t_start as u32);
        t_end.push(h.t_end as u32);
        strand.push(h.strand.to_string());
        energy.push(f64::from(h.energy));
        alignment.push(h.alignment.as_ref().map(|a| a.fingerprint()));
    }

    let columns: Vec<ArrayRef> = vec![
        Arc::new(UInt32Array::from(query_idx)),
        Arc::new(UInt32Array::from(target_idx)),
        Arc::new(UInt32Array::from(q_start)),
        Arc::new(UInt32Array::from(q_end)),
        Arc::new(UInt32Array::from(t_start)),
        Arc::new(UInt32Array::from(t_end)),
        Arc::new(StringArray::from(strand)),
        Arc::new(Float64Array::from(energy)),
        Arc::new(StringArray::from(alignment)),
    ];

    RecordBatch::try_new(schema.clone(), columns).expect("column count and types match schema")
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

        let name = CString::new("arrow_array_stream").unwrap();
        PyCapsule::new(py, stream, Some(name))
    }
}

// =============================================================================
// PyTargetStore
// =============================================================================

/// An in-memory index of target sequences for RNA-RNA interaction search.
///
/// Build once with `build_index()`, then reuse across many `search()` calls.
#[pyclass(name = "TargetStore")]
pub struct PyTargetStore(TargetStore);

#[pymethods]
impl PyTargetStore {
    /// Load a pre-built index from disk.
    ///
    /// Parameters
    /// ----------
    /// path : str | os.PathLike
    ///     Path to the `.idx` file produced by `build_index()`.
    #[staticmethod]
    fn open(path: PathBuf) -> PyResult<Self> {
        let store = TargetStore::open(&path)?;
        Ok(PyTargetStore(store))
    }

    fn __repr__(&self) -> String {
        format!("TargetStore(targets={})", self.0.len())
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
    py.allow_threads(|| TargetStore::build(&fasta, &output))?;
    Ok(())
}

/// Search for RNA-RNA interactions between query sequences and an indexed target set.
///
/// Parameters
/// ----------
/// query_fasta : list[str | os.PathLike]
///     One or more FASTA files containing query sequences.
/// store : TargetStore
///     Pre-loaded target index (from `TargetStore.open()`).
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
))]
fn search(
    py: Python<'_>,
    query_fasta: Vec<PathBuf>,
    store: &PyTargetStore,
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
    max_extension: u8,
    seed_energy: f64,
    no_max_prune: bool,
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
            penalty: Energy::try_from(penalty)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?,
            temperature,
        },
        extend: ExtendConfig {
            max_extension,
        },
        filter: FilterConfig {
            delta_g: Energy::try_from(energy_threshold)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?,
            seed_energy: Energy::try_from(seed_energy)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?,
            no_max_prune,
        },
        output: OutputConfig {
            format: OutputFormat::Detailed,
            compress: OutputCompression::None,
            multifile: false,
        },
    };

    let paths: Vec<&std::path::Path> = query_fasta.iter().map(|p| p.as_path()).collect();
    let queries = QueryRegistry::from_fastas(&paths, &config.seed)?;

    let hits = py.allow_threads(|| run_search_in_memory(&queries, &store.0, &config))?;
    let schema = search_result_schema().clone();
    let batch = hits_to_record_batch(hits, &schema);
    Ok(PySearchResult::new(batch, schema))
}

// =============================================================================
// Module registration
// =============================================================================

#[pymodule]
fn _risearch(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTargetStore>()?;
    m.add_function(wrap_pyfunction!(build_index, m)?)?;
    m.add_function(wrap_pyfunction!(search, m)?)?;
    Ok(())
}
