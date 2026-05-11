use crate::config::SeedConfig;
use crate::index::store::TargetStore;
use crate::registry::QueryRegistry;

use super::parallel_sa::{traverse, SeedMatch};
use super::{SeedHit, SeedView};

pub struct SeedingEngine<'a> {
    queries: &'a QueryRegistry,
    qview: SeedView<'a>,
    tview: SeedView<'a>,
}

impl<'a> SeedingEngine<'a> {
    pub fn new(queries: &'a QueryRegistry, targets: &'a TargetStore) -> Self {
        Self {
            queries,
            qview: queries.view(),
            tview: targets.view(),
        }
    }

    pub fn run(&self, config: &SeedConfig) -> Vec<(usize, Vec<SeedHit>)> {
        if self.qview.sa_real_len == 0 {
            return Vec::new();
        }

        let Some((global_min, global_max)) = self.global_seed_bounds() else {
            return Vec::new();
        };

        let mut seeds_by_query: Vec<Vec<SeedHit>> =
            (0..self.queries.len()).map(|_| Vec::new()).collect();
        if config.seed_wobble {
            self.traverse::<true>(config, global_min, global_max, &mut seeds_by_query);
        } else {
            self.traverse::<false>(config, global_min, global_max, &mut seeds_by_query);
        }

        seeds_by_query
            .into_iter()
            .enumerate()
            .filter_map(|(qi, seeds)| (!seeds.is_empty()).then_some((qi, seeds)))
            .collect()
    }

    fn global_seed_bounds(&self) -> Option<(usize, usize)> {
        let mut global_min = usize::MAX;
        let mut global_max = 0usize;
        for q in self.queries.entries() {
            global_min = global_min.min(q.min_seed_len);
            global_max = global_max.max(q.max_seed_len);
        }
        (global_min <= global_max).then_some((global_min, global_max))
    }

    fn traverse<const WOBBLE: bool>(
        &self,
        config: &SeedConfig,
        global_min: usize,
        global_max: usize,
        seeds_by_query: &mut [Vec<SeedHit>],
    ) {
        traverse::<WOBBLE, _>(
            self.qview,
            self.tview,
            global_min,
            global_max,
            config.mismatch.max_mismatches,
            config.mismatch.min_prefix_matches,
            config.mismatch.min_suffix_matches,
            &mut |m| self.materialize_seed_match(seeds_by_query, m),
        );
    }

    fn materialize_seed_match(&self, seeds_by_query: &mut [Vec<SeedHit>], raw_match: SeedMatch) {
        let seed_len = raw_match.seed_len;
        for &query_sa_pos in
            &self.qview.combined_sa[raw_match.query_interval.start..raw_match.query_interval.end]
        {
            let Some((query_idx, query_local_pos)) =
                remap(self.qview.offsets, query_sa_pos as usize)
            else {
                continue;
            };
            if query_local_pos + seed_len > self.qview.seq_lens[query_idx] {
                continue;
            }

            let query = self.queries.get(query_idx);
            if seed_len < query.min_seed_len || seed_len > query.max_seed_len {
                continue;
            }

            let query_start = query.seed_interval.start + query_local_pos;
            if query.has_n_any() && has_n_in_range(query.n_prefix(), query_start, seed_len) {
                continue;
            }

            for &target_sa_pos in &self.tview.combined_sa
                [raw_match.target_interval.start..raw_match.target_interval.end]
            {
                let Some((target_idx, target_local_pos)) =
                    remap(self.tview.offsets, target_sa_pos as usize)
                else {
                    continue;
                };
                let Some((strand, target_start)) = TargetStore::map_target_pos(
                    target_local_pos,
                    self.tview.seq_lens[target_idx],
                    seed_len,
                ) else {
                    continue;
                };

                seeds_by_query[query_idx].push(SeedHit {
                    query_idx,
                    query_start,
                    target_idx,
                    target_start,
                    len: seed_len,
                    strand,
                });
            }
        }
    }
}

/// Check if any position in range [start, start+len) contains 'N'.
#[inline]
fn has_n_in_range(n_prefix: &[u32], start: usize, len: usize) -> bool {
    n_prefix[start + len] != n_prefix[start]
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

    use crate::config::{MismatchSpec, SeedConfig, SeedSpec};
    use crate::index::store::TargetStore;
    use crate::registry::QueryRegistry;
    use crate::types::Strand;

    use super::*;

    fn build_store(fasta: &str) -> (TargetStore, tempfile::TempDir) {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("targets.fa");
        let index_path = dir.path().join("targets.idx");
        std::fs::write(&fasta_path, fasta).unwrap();
        TargetStore::build(&fasta_path, &index_path).unwrap();
        (TargetStore::open(&index_path).unwrap(), dir)
    }

    fn build_queries(fasta: &str, config: &SeedConfig) -> QueryRegistry {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(fasta.as_bytes()).unwrap();
        QueryRegistry::from_fasta(file.path(), config).unwrap()
    }

    #[test]
    fn collect_maps_seed_interval_back_to_full_query_coordinates() {
        let config = SeedConfig::with_wobble(
            SeedSpec::Interval {
                start: 3,
                end: 4,
                length: Some(2),
            },
            MismatchSpec::exact(),
            false,
        );
        let queries = build_queries(">q1\nGGAC\n", &config);
        let (targets, _dir) = build_store(">t1\nGU\n");

        let groups = SeedingEngine::new(&queries, &targets).run(&config);
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
    fn collect_normalizes_reverse_strand_target_hits() {
        let config = SeedConfig::with_wobble(SeedSpec::LengthOnly(2), MismatchSpec::exact(), false);
        let queries = build_queries(">q1\nAC\n", &config);
        let (targets, _dir) = build_store(">t1\nAC\n");

        let groups = SeedingEngine::new(&queries, &targets).run(&config);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].1.len(), 1);

        let seed = &groups[0].1[0];
        assert_eq!(seed.query_start, 0);
        assert_eq!(seed.target_start, 0);
        assert_eq!(seed.len, 2);
        assert_eq!(seed.strand, Strand::Reverse);
    }
}
