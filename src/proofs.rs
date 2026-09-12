//! Kani harnesses for small arithmetic and coordinate invariants.
//!
//! This module is included only by the crate root under `cfg(kani)`.  The
//! harnesses intentionally use the production functions and representations;
//! they do not reimplement the properties being checked.

/// `flat_idx` must map every valid DSM coordinate into its six-by-six-by-six
/// by-six storage table.
#[kani::proof]
fn flat_idx_is_bounded_for_every_dsm_coordinate() {
    let q1: u8 = kani::any();
    let q2: u8 = kani::any();
    let t1: u8 = kani::any();
    let t2: u8 = kani::any();
    kani::assume(q1 < 6);
    kani::assume(q2 < 6);
    kani::assume(t1 < 6);
    kani::assume(t2 < 6);

    let index = crate::dsm::flat_idx(q1, q2, t1, t2);
    assert!(index < 6 * 6 * 6 * 6);
}

/// `BestScore` adds one bounded DP value to one bounded transition.  The
/// bounds are the production DP contract: at most two `MAX_EXT`-sized axes
/// contribute a transition at `MAX_TRANSITION_SCORE` per step.
#[kani::proof]
fn best_score_update_is_safe_and_selects_only_strict_improvements() {
    const MAX_PATH_SCORE: i32 =
        (2 * crate::dp::MAX_EXT as i64 * crate::dp::MAX_TRANSITION_SCORE) as i32;
    const MAX_TRANSITION: i32 = crate::dp::MAX_TRANSITION_SCORE as i32;

    let initial: i32 = kani::any();
    let value: i32 = kani::any();
    let term: i32 = kani::any();
    let q_idx: usize = kani::any();
    let t_idx: usize = kani::any();
    kani::assume(value >= -MAX_PATH_SCORE);
    kani::assume(value <= MAX_PATH_SCORE);
    kani::assume(term >= -MAX_TRANSITION);
    kani::assume(term <= MAX_TRANSITION);

    let mut best = crate::dp::BestScore {
        score: initial,
        q_idx: 0,
        t_idx: 0,
    };
    let candidate = value + term;
    best.update_if_better(value, term, q_idx, t_idx);

    if candidate > initial {
        assert_eq!(best.score, candidate);
        assert_eq!(best.q_idx, q_idx);
        assert_eq!(best.t_idx, t_idx);
    } else {
        assert_eq!(best.score, initial);
        assert_eq!(best.q_idx, 0);
        assert_eq!(best.t_idx, 0);
    }
}
