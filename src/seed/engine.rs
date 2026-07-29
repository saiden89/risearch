use anyhow::{Context, Result};

use crate::config::SeedConfig;
use crate::index::sa::SuffixArray;
use crate::index::store::{TargetRegistry, SA_CHAR_PADDING};
use crate::index::RegistryView;
use crate::registry::{Query, QueryRegistry};
use crate::types::Base;

use super::parallel_sa::{traverse, SeedMatch};
use super::SeedHit;

pub struct SeedingEngine<'a> {
    queries: &'a QueryRegistry,
    tview: RegistryView<'a>,
}

impl<'a> SeedingEngine<'a> {
    pub fn new(queries: &'a QueryRegistry, targets: &'a TargetRegistry) -> Self {
        Self {
            queries,
            tview: targets.view(),
        }
    }

    /// Collect every seed grouped by query, materializing them all in memory.
    /// Retained for benchmarks and unit tests; the search paths stream one query
    /// at a time via [`seed_query`](Self::seed_query).
    pub fn run(&self, config: &SeedConfig) -> Result<Vec<(usize, Vec<SeedHit>)>> {
        let mut groups = Vec::new();
        for qi in 0..self.queries.len() {
            let mut seeds = Vec::new();
            self.seed_query(qi, config, |seed| seeds.push(seed))?;
            if !seeds.is_empty() {
                groups.push((qi, seeds));
            }
        }
        Ok(groups)
    }

    /// Seed a single query against the shared target SA: build the query's own
    /// (tiny) suffix array, traverse, emit each hit to `on_seed`, then drop it.
    /// Every seed for query `qi` is produced here, so a caller can own one query
    /// end-to-end (per-query parallelism + local dedup).
    pub fn seed_query<F: FnMut(SeedHit)>(
        &self,
        qi: usize,
        config: &SeedConfig,
        mut on_seed: F,
    ) -> Result<()> {
        // `traverse` is monomorphized on whether G-U wobble pairs are allowed.
        let query = self.queries.get(qi);
        if config.seed_wobble {
            self.seed_one::<true, _>(qi, query, config, &mut on_seed)
        } else {
            self.seed_one::<false, _>(qi, query, config, &mut on_seed)
        }
    }

    fn seed_one<const WOBBLE: bool, F: FnMut(SeedHit)>(
        &self,
        qi: usize,
        query: &Query,
        config: &SeedConfig,
        on_seed: &mut F,
    ) -> Result<()> {
        // Returning here also skips the per-query suffix array build.
        if query.min_seed_len > query.max_seed_len {
            return Ok(());
        }
        let tview = self.tview;

        // SA_CHAR_PADDING sentinels keep the branchless kernel lookups in bounds.
        let seed = query.seed_sequence();
        // A per-query seed is tens of bases; an OpenMP team would cost more than the sort.
        let mut sa = SuffixArray::build(seed, Some(1))
            .with_context(|| format!("building suffix array for query '{}'", query.name()))?
            .into_inner();
        // Suffixes too short for min_seed_len emit nothing; dropping them avoids
        // wasted target-SA descent.
        let max_valid_start = seed.len().saturating_sub(query.min_seed_len);
        sa.retain(|&p| (p as usize) <= max_valid_start);
        let sa_real_len = sa.len();
        sa.resize(sa.len() + SA_CHAR_PADDING, 0u64);
        let mut seq = vec![Base::Gap; seed.len() + SA_CHAR_PADDING];
        seq[..seed.len()].copy_from_slice(seed);
        let seq_lens = [seed.len()];
        let qview = RegistryView {
            combined_seq: &seq,
            combined_sa: &sa,
            sa_real_len,
            offsets: &[0],
            seq_lens: &seq_lens,
        };

        traverse::<WOBBLE, _>(
            qview,
            tview,
            query.min_seed_len,
            query.max_seed_len,
            config.max_mismatches,
            config.min_prefix_matches,
            config.min_suffix_matches,
            &mut |m| {
                emit_seed_match::<WOBBLE, _>(
                    qi,
                    query,
                    qview.combined_sa,
                    tview,
                    config.no_max_prune,
                    m,
                    &mut *on_seed,
                )
            },
        );
        Ok(())
    }
}

/// Emit every (query position × target position) seed for one `SeedMatch` of a
/// single query. Query positions come straight out of the single-query SA and
/// map back through [`Query::map_seed_pos`]; only the target side needs the
/// offset binary search.
///
/// Non-maximal seeds are dropped here unless `no_max_prune`: they are shorter
/// copies of a longer match, so extending them only rediscovers the same duplex.
fn emit_seed_match<const WOBBLE: bool, F: FnMut(SeedHit)>(
    qi: usize,
    query: &Query,
    query_sa: &[u64],
    tview: RegistryView<'_>,
    no_max_prune: bool,
    raw_match: SeedMatch,
    on_seed: &mut F,
) {
    let seed_len = raw_match.seed_len;
    let query_bases = query.sequence();
    let seed_interval = query.seed_interval();
    for &query_sa_pos in &query_sa[raw_match.query_interval.start..raw_match.query_interval.end] {
        let Some(query_start) = query.map_seed_pos(query_sa_pos as usize, seed_len) else {
            continue;
        };

        for &target_sa_pos in
            &tview.combined_sa[raw_match.target_interval.start..raw_match.target_interval.end]
        {
            let Some((target_idx, target_local_pos)) = remap(tview.offsets, target_sa_pos as usize)
            else {
                continue;
            };
            // Preserve block-local duplex coordinates; only identify the strand.
            let Some((strand, target_start)) = TargetRegistry::map_target_pos(
                target_local_pos,
                tview.seq_lens[target_idx],
                seed_len,
            ) else {
                continue;
            };

            let hit = SeedHit {
                query_idx: qi,
                query_start,
                target_idx,
                target_start,
                len: seed_len,
                strand,
            };
            if no_max_prune
                || hit.is_maximal(
                    query_bases,
                    tview.target(target_idx, strand),
                    &seed_interval,
                    WOBBLE,
                )
            {
                on_seed(hit);
            }
        }
    }
}

/// Remap a global SA position to an index via binary search on an offset table.
/// Returns `None` if the position falls before the first entry.
#[inline]
fn remap(offsets: &[usize], global_pos: usize) -> Option<(usize, usize)> {
    let idx = offsets
        .partition_point(|&o| o <= global_pos)
        .checked_sub(1)?;
    Some((idx, global_pos - offsets[idx]))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use crate::config::SeedConfig;
    use crate::index::store::TargetRegistry;
    use crate::registry::QueryRegistry;
    use crate::types::Strand;

    use super::*;

    fn build_store(fasta: &str) -> (TargetRegistry, tempfile::TempDir) {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");
        fs_err::write(&fasta_path, fasta).unwrap();
        TargetRegistry::build(&fasta_path, &index_path, None).unwrap();
        (TargetRegistry::open(&index_path).unwrap(), dir)
    }

    fn build_queries(fasta: &str, config: &SeedConfig) -> QueryRegistry {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(fasta.as_bytes()).unwrap();
        QueryRegistry::from_fasta(file.path(), config).unwrap()
    }

    #[test]
    fn collect_maps_seed_interval_back_to_full_query_coordinates() {
        let config = SeedConfig {
            seed_start: Some(3),
            seed_end: Some(4),
            seed_length: Some(2),
            seed_wobble: false,
            no_max_prune: true,
            ..Default::default()
        };
        let queries = build_queries(">q1\nGGAC\n", &config);
        let (targets, _dir) = build_store(">t1\nGU\n");

        let groups = SeedingEngine::new(&queries, &targets).run(&config).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0, 0);
        assert_eq!(groups[0].1.len(), 1);

        let seed = &groups[0].1[0];
        assert_eq!(seed.query_idx, 0);
        assert_eq!(seed.query_start, 2);
        assert_eq!(seed.target_idx, 0);
        assert_eq!(seed.target_start, 0);
        assert_eq!(seed.len, 2);
        assert_eq!(seed.strand, Strand::Forward);
    }

    #[test]
    fn collect_preserves_same_block_coordinates_on_both_strands() {
        let config = SeedConfig {
            seed_length: Some(2),
            seed_wobble: false,
            no_max_prune: true,
            ..Default::default()
        };
        let queries = build_queries(">q1\nCG\n", &config);
        let (targets, _dir) = build_store(">t1\nAACGU\n");

        let groups = SeedingEngine::new(&queries, &targets).run(&config).unwrap();
        assert_eq!(groups.len(), 1);
        let seeds = &groups[0].1;
        assert_eq!(seeds.len(), 2);
        assert!(seeds
            .iter()
            .any(|seed| seed.strand == Strand::Forward && seed.target_start == 1));
        assert!(seeds
            .iter()
            .any(|seed| seed.strand == Strand::Reverse && seed.target_start == 2));
    }
}
