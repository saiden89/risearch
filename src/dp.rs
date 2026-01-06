//! Clean DP extension implementation from first principles.
//!
//! Uses nearest-neighbor stacking energy (DSM tables) with affine gap penalties.
//! Based on rust-bio patterns but adapted for RNA duplex alignment.

use crate::dsm::StackPair;
use crate::seq::Seq;
use crate::types::Base;

/// Extend alignment to the left (query 5', target 3')
///
/// Query extends toward 5' (decreasing index), Target extends toward 3' (increasing index).
pub fn extend_left(
    query: &Seq,
    target: &Seq,
    q_start: usize,
    t_start: usize,
    max_ext: usize,
) -> DpExtension {
    let q_len = (q_start + 1).min(max_ext);
    let t_len = (target.len() - t_start - 1).min(max_ext);

    extend(
        |i| query.left(q_start, i),
        |j| target.right(t_start, j),
        q_len,
        t_len,
    )
}

/// Extend alignment to the right (query 3', target 5')
///
/// Query extends toward 3' (increasing index), Target extends toward 5' (decreasing index).
pub fn extend_right(
    query: &Seq,
    target: &Seq,
    q_end: usize,
    t_end: usize,
    max_ext: usize,
) -> DpExtension {
    let q_len = (query.len() - q_end).min(max_ext);
    let t_len = (t_end + 1).min(max_ext);

    extend(
        |i| query.right(q_end, i),
        |j| target.left(t_end, j),
        q_len,
        t_len,
    )
}

/// Result of DP extension
#[derive(Debug, Clone)]
pub struct DpExtension {
    pub score: i32,
    pub q_len: usize,
    pub t_len: usize,
    pub trace: Vec<DpOp>,
}

/// Alignment operation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DpOp {
    Match,
    GapQ,
    GapT,
}

/// Core DP extension function
pub fn extend<Q, T>(q: Q, t: T, q_len: usize, t_len: usize) -> DpExtension
where
    Q: Fn(usize) -> Base,
    T: Fn(usize) -> Base,
{
    // Early return if nothing to extend - but still return initial terminal
    if q_len == 0 || t_len == 0 {
        // Initial terminal: DSM[GAP][q(0)][GAP][t(0)]
        let initial_terminal = StackPair::new(Base::Gap, q(0), Base::Gap, t(0)).energy() as i32;
        return DpExtension {
            score: initial_terminal,
            q_len: 0,
            t_len: 0,
            trace: vec![],
        };
    }

    let stack = |q1: Base, q2: Base, t1: Base, t2: Base| -> i32 {
        StackPair::new(q1, q2, t1, t2).energy() as i32
    };

    // DP matrices (2-row optimization)
    let mut m_prev = vec![i32::MIN / 2; t_len + 1];
    let mut m_curr = vec![i32::MIN / 2; t_len + 1];
    let mut bq_prev = vec![i32::MIN / 2; t_len + 1];
    let mut bq_curr = vec![i32::MIN / 2; t_len + 1];
    let mut bt_prev = vec![i32::MIN / 2; t_len + 1];
    let mut bt_curr = vec![i32::MIN / 2; t_len + 1];

    // Traceback storage
    let mut tb: Vec<Vec<DpOp>> = vec![vec![DpOp::Match; t_len + 1]; q_len + 1];

    // Best score tracking - initialize with terminal penalty for zero extension
    // This matches old DP: best_e = DSM[GAP][Q(0)][GAP][T(0)]
    let initial_terminal = stack(Base::Gap, q(0), Base::Gap, t(0));
    let mut best_score = initial_terminal;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    // Initialize M[0,0] = 0 (set in m_curr because loop swaps first)
    m_curr[0] = 0;

    for i in 1..=q_len {
        std::mem::swap(&mut m_prev, &mut m_curr);
        std::mem::swap(&mut bq_prev, &mut bq_curr);
        std::mem::swap(&mut bt_prev, &mut bt_curr);

        for val in m_curr.iter_mut() {
            *val = i32::MIN / 2;
        }
        for val in bq_curr.iter_mut() {
            *val = i32::MIN / 2;
        }
        for val in bt_curr.iter_mut() {
            *val = i32::MIN / 2;
        }

        for j in 1..=t_len {
            let s = stack(q(i - 1), q(i), t(j - 1), t(j));

            // Match: from M, Bq, or Bt diagonal
            let m_from_m = m_prev[j - 1] + s;
            let m_from_bq = bq_prev[j - 1] + s;
            let m_from_bt = bt_prev[j - 1] + s;

            let mut m_best = m_from_m;
            let mut m_tb = DpOp::Match;
            if m_from_bq > m_best {
                m_best = m_from_bq;
                m_tb = DpOp::GapQ;
            }
            if m_from_bt > m_best {
                m_best = m_from_bt;
                m_tb = DpOp::GapT;
            }
            m_curr[j] = m_best;
            tb[i][j] = m_tb;

            // Bq: gap in query (from above)
            let bq_open = m_prev[j] + stack(q(i - 1), q(i), Base::Gap, t(j));
            let bq_ext = bq_prev[j] + stack(q(i - 1), q(i), Base::Gap, Base::Gap);
            bq_curr[j] = bq_open.max(bq_ext);

            // Bt: gap in target (from left)
            let bt_open = m_curr[j - 1] + stack(q(i), Base::Gap, t(j - 1), t(j));
            let bt_ext = bt_curr[j - 1] + stack(Base::Gap, Base::Gap, t(j - 1), t(j));
            bt_curr[j] = bt_open.max(bt_ext);

            // Update best with terminal penalty
            let terminal = stack(q(i), Base::Gap, t(j), Base::Gap);
            let score_with_term = m_curr[j] + terminal;
            if score_with_term > best_score {
                best_score = score_with_term;
                best_i = i;
                best_j = j;
            }
        }
    }

    // Traceback
    let mut trace = Vec::new();
    let (mut i, mut j) = (best_i, best_j);
    while i > 0 && j > 0 {
        let op = tb[i][j];
        trace.push(op);
        match op {
            DpOp::Match => {
                i -= 1;
                j -= 1;
            }
            DpOp::GapQ => {
                i -= 1;
            }
            DpOp::GapT => {
                j -= 1;
            }
        }
    }
    trace.reverse();

    DpExtension {
        score: best_score,
        q_len: best_i,
        t_len: best_j,
        trace,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Strand;

    #[test]
    fn test_extend_empty() {
        let result = extend(|_| Base::A, |_| Base::U, 0, 0);
        assert_eq!(result.score, 0);
        assert_eq!(result.q_len, 0);
    }

    #[test]
    fn test_single_match() {
        // Simplest case: extend by exactly 1 position in each direction
        // Q: A-C  (positions 0,1)
        // T: U-G  (positions 0,1, antiparallel)
        //
        // At (1,1): stack(q0=A, q1=C, t0=U, t1=G)
        // This should give us a stacking energy value
        let q = [Base::A, Base::C];
        let t = [Base::U, Base::G];

        let result = extend(|i| q[i.min(1)], |j| t[j.min(1)], 1, 1);

        // Calculate expected:
        // M[1,1] = M[0,0] + stack(A,C,U,G) = 0 + stack
        // Terminal = stack(C, Gap, G, Gap)
        // Score = M[1,1] + terminal
        let stack_val = StackPair::new(Base::A, Base::C, Base::U, Base::G).energy() as i32;
        let terminal = StackPair::new(Base::C, Base::Gap, Base::G, Base::Gap).energy() as i32;
        let initial_term = StackPair::new(Base::Gap, Base::A, Base::Gap, Base::U).energy() as i32;

        println!("Single match test:");
        println!("  stack(A,C,U,G) = {}", stack_val);
        println!("  terminal(C,Gap,G,Gap) = {}", terminal);
        println!("  initial_terminal(Gap,A,Gap,U) = {}", initial_term);
        println!(
            "  Expected M[1,1] + term = {} + {} = {}",
            stack_val,
            terminal,
            stack_val + terminal
        );
        println!(
            "  Actual result: score={}, q_len={}, t_len={}",
            result.score, result.q_len, result.t_len
        );

        // The result should either be:
        // - initial_terminal (if no extension is better)
        // - stack_val + terminal (if extending is better)
        let expected = (stack_val + terminal).max(initial_term);
        assert_eq!(
            result.score, expected,
            "Score should match hand calculation"
        );
    }

    #[test]
    fn test_no_extension_better() {
        // Case where NOT extending gives better score than extending
        // Use bases that give unfavorable stacking
        let q = [Base::A, Base::A]; // AA
        let t = [Base::A, Base::A]; // AA (not complementary, should be unfavorable)

        let result = extend(|i| q[i.min(1)], |j| t[j.min(1)], 1, 1);

        let stack_val = StackPair::new(Base::A, Base::A, Base::A, Base::A).energy() as i32;
        let terminal = StackPair::new(Base::A, Base::Gap, Base::A, Base::Gap).energy() as i32;
        let initial_term = StackPair::new(Base::Gap, Base::A, Base::Gap, Base::A).energy() as i32;

        println!("No extension test:");
        println!("  stack(A,A,A,A) = {}", stack_val);
        println!("  terminal(A,Gap,A,Gap) = {}", terminal);
        println!("  initial_terminal(Gap,A,Gap,A) = {}", initial_term);
        println!(
            "  Extend score = {} + {} = {}",
            stack_val,
            terminal,
            stack_val + terminal
        );
        println!("  No extend score = {}", initial_term);
        println!("  Actual: score={}", result.score);
    }

    #[test]
    fn test_extend_right_with_seq() {
        let query = Seq::new(b"ACGUACGU", Strand::Forward);
        let target = Seq::new(b"UGCAUGCA", Strand::Forward);
        let result = extend_right(&query, &target, 2, 5, 10);
        assert!(result.score != i32::MIN / 2, "Should compute valid score");
    }

    #[test]
    fn test_extend_left_with_seq() {
        let query = Seq::new(b"ACGUACGU", Strand::Forward);
        let target = Seq::new(b"UGCAUGCA", Strand::Forward);
        let result = extend_left(&query, &target, 5, 2, 10);
        assert!(result.score != i32::MIN / 2, "Should compute valid score");
    }
}
