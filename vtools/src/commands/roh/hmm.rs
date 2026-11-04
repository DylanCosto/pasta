//! A hidden Markov model with distance-dependent transitions: a port of bcftools' `HMM.c`.
//!
//! The number of states `N` is a type parameter. A transition matrix `m` is `[[f64; N]; N]`
//! with `m[i][j]` being `HMM.c`'s `MAT(m, N, i, j)`: the probability of state `i` given state
//! `j` (the transition `j -> i`). The arithmetic follows `HMM.c` step for step, so results
//! match bcftools to the last bit.

/// A transition matrix (see the module docs).
pub type Matrix<const N: usize> = [[f64; N]; N];

/// A model state saved for a sliding HMM (`HMM.c`'s `snapshot_t`): the runs that are given it
/// fill it in at its position, and a later run can start from there.
#[derive(Clone, Copy, Debug)]
pub struct Snapshot<const N: usize> {
    /// The position the snapshot is taken at, 0 when inactive.
    snap_at_pos: u32,
    vit_prob: [f64; N],
    fwd_prob: [f64; N],
}

impl<const N: usize> Snapshot<N> {
    /// An inactive snapshot.
    pub fn new() -> Self {
        Self {
            snap_at_pos: 0,
            vit_prob: [0.0; N],
            fwd_prob: [0.0; N],
        }
    }

    /// Takes the snapshot at position `pos` of the next runs that fill it (`hmm_snapshot`).
    pub fn take_at(&mut self, pos: u32) {
        self.snap_at_pos = pos;
    }

    /// Deactivates the snapshot, so that runs start from the initial state (`hmm_reset`).
    pub fn deactivate(&mut self) {
        self.snap_at_pos = 0;
    }
}

impl<const N: usize> Default for Snapshot<N> {
    fn default() -> Self {
        Self::new()
    }
}

/// Initial or current state probabilities (`HMM.c`'s `init` and `state`).
#[derive(Clone, Copy, Debug)]
struct StateProbs<const N: usize> {
    snap_at_pos: u32,
    vit_prob: [f64; N],
    fwd_prob: [f64; N],
    bwd_prob: [f64; N],
}

/// Changes the transition matrix of one step (`HMM.c`'s `set_tprob_f`), given the positions of
/// the previous and the current site.
pub trait TransitionHook<const N: usize> {
    /// Adjusts `tprob` for the step from `prev_pos` to `pos`.
    fn set_tprob(&mut self, prev_pos: u32, pos: u32, tprob: &mut Matrix<N>);
}

/// Transitions that depend on the distance only through the precalculated matrix powers.
#[cfg(test)]
#[derive(Clone, Copy, Debug, Default)]
pub struct NoHook;

#[cfg(test)]
impl<const N: usize> TransitionHook<N> for NoHook {
    #[inline]
    fn set_tprob(&mut self, _prev_pos: u32, _pos: u32, _tprob: &mut Matrix<N>) {}
}

/// `HMM.c`'s `hmm_t`.
#[derive(Clone, Debug)]
pub struct Hmm<const N: usize> {
    vprob: [f64; N],
    vprob_tmp: [f64; N],
    /// The Viterbi path: per site, the best previous state of each state; after the traceback,
    /// entry 0 holds the state of the site.
    vpath: Vec<[u8; N]>,
    bwd: [f64; N],
    bwd_tmp: [f64; N],
    /// Forward probabilities, after an initial entry.
    fwd: Vec<[f64; N]>,
    /// The number of precalculated matrices, 0 for constant transitions (`ntprob_arr`).
    ntprob_arr: usize,
    /// The current step's transitions (`curr_tprob`).
    curr_tprob: Matrix<N>,
    /// The transition matrix raised to the powers 1 to `max(ntprob_arr, 1)`.
    tprob_arr: Vec<Matrix<N>>,
    init: StateProbs<N>,
    state: StateProbs<N>,
}

/// `HMM.c`'s `multiply_matrix`.
#[inline]
fn multiply_matrix<const N: usize>(a: &Matrix<N>, b: &Matrix<N>) -> Matrix<N> {
    let mut out = [[0.0; N]; N];
    for i in 0..N {
        for j in 0..N {
            let mut val = 0.0;
            for k in 0..N {
                val += a[i][k] * b[k][j];
            }
            out[i][j] = val;
        }
    }
    out
}

#[expect(
    clippy::needless_range_loop,
    reason = "the loops follow HMM.c's, whose order of operations decides the results"
)]
impl<const N: usize> Hmm<N> {
    /// Creates a model with the transition matrix `tprob` and `ntprob` precalculated powers of
    /// it, 0 for transitions independent of the distance (`hmm_init`).
    pub fn new(tprob: &Matrix<N>, ntprob: usize) -> Self {
        let state = StateProbs {
            snap_at_pos: 0,
            vit_prob: [0.0; N],
            fwd_prob: [0.0; N],
            bwd_prob: [0.0; N],
        };
        let mut hmm = Self {
            vprob: [0.0; N],
            vprob_tmp: [0.0; N],
            vpath: Vec::new(),
            bwd: [0.0; N],
            bwd_tmp: [0.0; N],
            fwd: Vec::new(),
            ntprob_arr: 0,
            curr_tprob: [[0.0; N]; N],
            tprob_arr: Vec::new(),
            init: state,
            state,
        };
        hmm.set_tprob(tprob, ntprob);
        hmm.init_states();
        hmm
    }

    /// Sets equal initial state probabilities (`hmm_init_states(hmm, NULL)`).
    fn init_states(&mut self) {
        self.init.snap_at_pos = 0;
        self.state.snap_at_pos = 0;
        self.init.vit_prob = [1.0 / N as f64; N];
        self.init.bwd_prob = [1.0; N];
        self.init.fwd_prob = self.init.vit_prob;
        self.state.vit_prob = self.init.vit_prob;
        self.state.fwd_prob = self.init.fwd_prob;
        self.state.bwd_prob = self.init.bwd_prob;
    }

    /// Sets the transition matrix and precalculates `ntprob` of its powers (`hmm_set_tprob`).
    pub fn set_tprob(&mut self, tprob: &Matrix<N>, ntprob: usize) {
        self.ntprob_arr = ntprob;
        self.tprob_arr.clear();
        self.tprob_arr.push(*tprob);
        for i in 1..ntprob.max(1) {
            let next = multiply_matrix(&self.tprob_arr[0], &self.tprob_arr[i - 1]);
            self.tprob_arr.push(next);
        }
    }

    /// Returns the transition matrix (the first of the precalculated ones, `hmm_get_tprob`).
    pub fn tprob(&self) -> &Matrix<N> {
        &self.tprob_arr[0]
    }

    /// Starts the next runs from a snapshot, or from the initial state when there is none or
    /// it is inactive (`hmm_restore`).
    pub fn restore(&mut self, snapshot: Option<&Snapshot<N>>) {
        match snapshot {
            Some(s) if s.snap_at_pos != 0 => {
                self.state.snap_at_pos = s.snap_at_pos;
                self.state.vit_prob = s.vit_prob;
                self.state.fwd_prob = s.fwd_prob;
            }
            _ => {
                self.state.snap_at_pos = 0;
                self.state.vit_prob = self.init.vit_prob;
                self.state.fwd_prob = self.init.fwd_prob;
            }
        }
    }

    /// Sets `curr_tprob` to the transitions over `pos_diff + 1` bases (`_set_tprob`).
    #[inline]
    fn set_curr_tprob(&mut self, pos_diff: usize) {
        if self.ntprob_arr == 0 {
            self.curr_tprob = self.tprob_arr[0];
            return;
        }
        self.curr_tprob = self.tprob_arr[pos_diff % self.ntprob_arr];
        let last = &self.tprob_arr[self.ntprob_arr - 1];
        for _ in 0..pos_diff / self.ntprob_arr {
            self.curr_tprob = multiply_matrix(last, &self.curr_tprob);
        }
    }

    /// Prepares the transitions of the step from `prev_pos` to `pos` (`pos >= prev_pos`),
    /// with the hook called as `set_tprob(hook_prev, hook_pos)`.
    #[inline]
    fn step_tprob(
        &mut self,
        prev_pos: u32,
        pos: u32,
        hook: &mut impl TransitionHook<N>,
        hook_args: (u32, u32),
    ) {
        // HMM.c computes this in `uint32_t` and asserts that it fits an `int`; positions are
        // sorted, so the unsigned difference is the distance (and 0 for a site at position 0
        // after one at -1, which wraps as in HMM.c).
        let pos_diff = if pos == prev_pos {
            0
        } else {
            pos.wrapping_sub(prev_pos).wrapping_sub(1) as usize
        };
        self.set_curr_tprob(pos_diff);
        hook.set_tprob(hook_args.0, hook_args.1, &mut self.curr_tprob);
    }

    /// Where the runs start: the snapshot's position, or the first site.
    fn start_pos(&self, sites: &[u32]) -> u32 {
        if self.state.snap_at_pos != 0 {
            self.state.snap_at_pos
        } else {
            sites[0]
        }
    }

    /// Runs the Viterbi algorithm over the sites, with the emission probabilities of each
    /// site's states, filling in `snapshot` at its position (`hmm_run_viterbi`).
    pub fn run_viterbi(
        &mut self,
        eprobs: &[[f64; N]],
        sites: &[u32],
        hook: &mut impl TransitionHook<N>,
        mut snapshot: Option<&mut Snapshot<N>>,
    ) {
        let n = sites.len();
        if self.vpath.len() < n {
            self.vpath.resize(n, [0; N]);
        }
        self.vprob = self.state.vit_prob;
        let mut prev_pos = self.start_pos(sites);

        for i in 0..n {
            self.step_tprob(prev_pos, sites[i], hook, (prev_pos, sites[i]));
            prev_pos = sites[i];
            let eprob = &eprobs[i];

            let mut vnorm = 0.0;
            for j in 0..N {
                let mut vmax = 0.0;
                let mut k_vmax = 0;
                for k in 0..N {
                    let pval = self.vprob[k] * self.curr_tprob[j][k];
                    if vmax < pval {
                        vmax = pval;
                        k_vmax = k;
                    }
                }
                self.vpath[i][j] = k_vmax as u8;
                self.vprob_tmp[j] = vmax * eprob[j];
                vnorm += self.vprob_tmp[j];
            }
            for v in &mut self.vprob_tmp {
                *v /= vnorm;
            }
            std::mem::swap(&mut self.vprob, &mut self.vprob_tmp);

            if let Some(s) = snapshot.as_deref_mut()
                && s.snap_at_pos == sites[i]
            {
                s.vit_prob = self.vprob;
            }
        }

        // The most likely last state, then the traceback.
        let mut iptr = 0;
        for i in 1..N {
            if self.vprob[iptr] < self.vprob[i] {
                iptr = i;
            }
        }
        for i in (0..n).rev() {
            let iptr_prev = self.vpath[i][iptr] as usize;
            self.vpath[i][0] = iptr as u8;
            iptr = iptr_prev;
        }
    }

    /// The state of site `i` on the Viterbi path of the last [`Self::run_viterbi`]
    /// (`HMM_VPATH`).
    #[inline]
    pub fn viterbi_state(&self, i: usize) -> u8 {
        self.vpath[i][0]
    }

    /// The posterior probability of `state` at site `i` after [`Self::run_fwd_bwd`]
    /// (`HMM_PPROB`).
    #[inline]
    pub fn posterior(&self, i: usize, state: usize) -> f64 {
        self.fwd[i + 1][state]
    }

    fn prepare_fwd_bwd(&mut self, n: usize) {
        if self.fwd.len() < n + 1 {
            self.fwd.resize(n + 1, [0.0; N]);
        }
        self.fwd[0] = self.state.fwd_prob;
        self.bwd = self.state.bwd_prob;
    }

    /// The forward pass shared by [`Self::run_fwd_bwd`] and [`Self::run_baum_welch`], filling
    /// in `snapshot` at its position.
    fn run_fwd(
        &mut self,
        eprobs: &[[f64; N]],
        sites: &[u32],
        hook: &mut impl TransitionHook<N>,
        mut snapshot: Option<&mut Snapshot<N>>,
    ) {
        let mut prev_pos = self.start_pos(sites);
        for i in 0..sites.len() {
            self.step_tprob(prev_pos, sites[i], hook, (prev_pos, sites[i]));
            prev_pos = sites[i];
            let eprob = &eprobs[i];
            let fwd_prev = self.fwd[i];
            let mut fwd = [0.0; N];

            let mut norm = 0.0;
            for j in 0..N {
                let mut pval = 0.0;
                for k in 0..N {
                    pval += fwd_prev[k] * self.curr_tprob[j][k];
                }
                fwd[j] = pval * eprob[j];
                norm += fwd[j];
            }
            for v in &mut fwd {
                *v /= norm;
            }
            self.fwd[i + 1] = fwd;

            if let Some(s) = snapshot.as_deref_mut()
                && s.snap_at_pos == sites[i]
            {
                s.fwd_prob = fwd;
            }
        }
    }

    /// One step of the backward pass: the new backward probabilities into `bwd_tmp`
    /// (normalized), from the site's emissions.
    #[inline]
    fn bwd_step(&mut self, eprob: &[f64; N]) {
        let mut bwd_norm = 0.0;
        for j in 0..N {
            let mut pval = 0.0;
            for k in 0..N {
                pval += self.bwd[k] * eprob[k] * self.curr_tprob[k][j];
            }
            self.bwd_tmp[j] = pval;
            bwd_norm += pval;
        }
        for v in &mut self.bwd_tmp {
            *v /= bwd_norm;
        }
    }

    /// Runs the forward-backward algorithm, filling in `snapshot` at its position
    /// (`hmm_run_fwd_bwd`).
    pub fn run_fwd_bwd(
        &mut self,
        eprobs: &[[f64; N]],
        sites: &[u32],
        hook: &mut impl TransitionHook<N>,
        snapshot: Option<&mut Snapshot<N>>,
    ) {
        let n = sites.len();
        self.prepare_fwd_bwd(n);
        self.run_fwd(eprobs, sites, hook, snapshot);

        let mut prev_pos = sites[n - 1];
        for i in 0..n {
            let site = sites[n - i - 1];
            self.step_tprob(site, prev_pos, hook, (site, prev_pos));
            prev_pos = site;
            self.bwd_step(&eprobs[n - i - 1]);

            // fwd now stores fwd*bwd.
            let fwd = &mut self.fwd[n - i];
            let mut norm = 0.0;
            for j in 0..N {
                fwd[j] *= self.bwd[j];
                norm += fwd[j];
            }
            for v in fwd.iter_mut() {
                *v /= norm;
            }
            std::mem::swap(&mut self.bwd, &mut self.bwd_tmp);
        }
    }

    /// Runs one iteration of the Baum-Welch algorithm and returns the new transition matrix;
    /// emission probabilities are not updated (`hmm_run_baum_welch`).
    pub fn run_baum_welch(
        &mut self,
        eprobs: &[[f64; N]],
        sites: &[u32],
        hook: &mut impl TransitionHook<N>,
    ) -> Matrix<N> {
        let n = sites.len();
        self.prepare_fwd_bwd(n);
        let mut tmp_xi = [[0.0; N]; N];
        let mut tmp_gamma = [0.0; N];
        let mut fwd_bwd = [0.0; N];

        self.run_fwd(eprobs, sites, hook, None);

        let mut prev_pos = sites[n - 1];
        for i in 0..n {
            let site = sites[n - i - 1];
            self.step_tprob(site, prev_pos, hook, (site, prev_pos));
            prev_pos = site;
            let eprob = &eprobs[n - i - 1];
            self.bwd_step(eprob);

            let fwd = &mut self.fwd[n - i];
            let mut norm = 0.0;
            for j in 0..N {
                fwd_bwd[j] = fwd[j] * self.bwd_tmp[j];
                norm += fwd_bwd[j];
            }
            for j in 0..N {
                fwd_bwd[j] /= norm;
                tmp_gamma[j] += fwd_bwd[j];
            }
            // HMM.c uses the base transition matrix here, not the step's.
            for j in 0..N {
                for k in 0..N {
                    tmp_xi[k][j] +=
                        fwd[j] * self.bwd[k] * self.tprob_arr[0][k][j] * eprob[k] / norm;
                }
            }
            *fwd = fwd_bwd;
            std::mem::swap(&mut self.bwd, &mut self.bwd_tmp);
        }
        for j in 0..N {
            let mut norm = 0.0;
            for k in 0..N {
                self.curr_tprob[k][j] = tmp_xi[k][j] / tmp_gamma[j];
                norm += self.curr_tprob[k][j];
            }
            for k in 0..N {
                self.curr_tprob[k][j] /= norm;
            }
        }
        self.curr_tprob
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: Matrix<2> = [[0.9, 0.2], [0.1, 0.8]];

    #[test]
    fn precalculated_powers_compose() {
        let hmm = Hmm::new(&T, 3);
        assert_eq!(hmm.tprob_arr[1], multiply_matrix(&T, &T));
        // Columns of a stochastic matrix's powers still sum to one.
        for m in &hmm.tprob_arr {
            assert!((m[0][0] + m[1][0] - 1.0).abs() < 1e-12, "{m:?}");
        }
    }

    #[test]
    fn far_steps_use_whole_blocks() {
        // A step over 7 bases with 3 precalculated powers is T^(7 % 3 + 1) * (T^3)^2 = T^8.
        let mut hmm = Hmm::new(&T, 3);
        hmm.set_curr_tprob(7);
        let mut want = T;
        for _ in 1..8 {
            want = multiply_matrix(&T, &want);
        }
        for (a, b) in hmm.curr_tprob.iter().flatten().zip(want.iter().flatten()) {
            assert!((a - b).abs() < 1e-12, "{:?} vs {want:?}", hmm.curr_tprob);
        }
    }

    #[test]
    fn runs_continue_from_a_snapshot() {
        let sites = [10, 20, 30, 40, 50];
        let eprobs = [[0.9, 0.1], [0.6, 0.4], [0.1, 0.9], [0.2, 0.8], [0.7, 0.3]];
        let mut full = Hmm::new(&T, 4);
        full.run_viterbi(&eprobs, &sites, &mut NoHook, None);
        full.run_fwd_bwd(&eprobs, &sites, &mut NoHook, None);

        // A first buffer of three sites, snapshotted at its second; the next buffer starts
        // there, as `flush_viterbi` does with overlapping buffers.
        let mut hmm = Hmm::new(&T, 4);
        let mut snap = Snapshot::new();
        snap.take_at(sites[1]);
        hmm.run_viterbi(&eprobs[..3], &sites[..3], &mut NoHook, Some(&mut snap));
        hmm.run_fwd_bwd(&eprobs[..3], &sites[..3], &mut NoHook, Some(&mut snap));
        hmm.restore(Some(&snap));
        hmm.run_viterbi(&eprobs[2..], &sites[2..], &mut NoHook, None);
        assert_eq!(hmm.vprob, full.vprob);
        hmm.run_fwd_bwd(&eprobs[2..], &sites[2..], &mut NoHook, None);
        // The last site has no later ones, so its posterior is the forward probability.
        assert_eq!(hmm.posterior(2, 1), full.posterior(4, 1));

        // An inactive snapshot starts from the initial state.
        snap.deactivate();
        hmm.restore(Some(&snap));
        hmm.run_viterbi(&eprobs, &sites, &mut NoHook, None);
        assert_eq!(hmm.vprob, full.vprob);
    }

    #[test]
    fn viterbi_follows_strong_emissions() {
        let mut hmm = Hmm::new(&T, 0);
        let sites = [10, 20, 30, 40];
        let eprobs = [[0.9, 0.1], [0.9, 0.1], [0.1, 0.9], [0.1, 0.9]];
        hmm.run_viterbi(&eprobs, &sites, &mut NoHook, None);
        let path: Vec<u8> = (0..4).map(|i| hmm.viterbi_state(i)).collect();
        assert_eq!(path, [0, 0, 1, 1]);
        hmm.run_fwd_bwd(&eprobs, &sites, &mut NoHook, None);
        for i in 0..4 {
            let sum = hmm.posterior(i, 0) + hmm.posterior(i, 1);
            assert!((sum - 1.0).abs() < 1e-12, "site {i}: {sum}");
        }
        assert!(hmm.posterior(0, 0) > 0.5 && hmm.posterior(3, 1) > 0.5);
    }
}
