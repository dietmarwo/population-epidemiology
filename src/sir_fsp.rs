//! Prepared finite state projection for the fixed SIR network.
//!
//! `crn_craft::fsp::solve` validates, parses and compiles a model and a run
//! document, discovers the reachable states and builds the uniformized
//! generator again for every rate binding. With positive normal rates, the
//! SIR reachable states, their discovery order, the generator edges and the
//! column order do not depend on the rate values. This module builds them
//! once. A candidate binds its two rates to that structure and runs the
//! solver's uniformization arithmetic in the same order, so means,
//! diagnostics and the FSP work count are bit-identical to the vendored
//! solver. Whenever a candidate reaches a condition that this kernel does not
//! reproduce, `solve` returns `None` and the caller runs the vendored solver.
use crn_craft::fsp::MAX_WORK;
use std::collections::HashMap;

/// The output budget the SIR run document declares.
pub const MAX_OUTPUT_BYTES: u64 = 8_388_608;
/// Species in declaration order: S, I and R.
const SPECIES: usize = 3;
/// Reactions in declaration order: infection and recovery.
const REACTIONS: usize = 2;

/// Rate-independent SIR projection structure.
pub struct Projection {
    /// Susceptible and infected counts per state, in CRN-craft discovery order.
    susceptible: Vec<f64>,
    infected: Vec<f64>,
    /// Target state of each state's infection and recovery transition.
    infection: Vec<Option<usize>>,
    recovery: Vec<Option<usize>>,
    /// Column-ordered generator: incoming entries of each target in
    /// ascending source order, with the diagonal at its own position and the
    /// escaping sink last (it has no incoming edges here).
    start: Vec<usize>,
    entries: Vec<(usize, Entry)>,
    edges: usize,
    times: Vec<f64>,
    tolerance: f64,
}

#[derive(Clone, Copy)]
enum Entry {
    Diagonal,
    Infection,
    Recovery,
}

/// Per-sample-time output of one rate binding.
pub struct Solution {
    /// Projected first moment of the infected count at each sample time.
    pub means: Vec<f64>,
    /// Physical leakage, series tail bound, mass defect and roundoff
    /// residual at each sample time, filled only when requested.
    pub errors: Vec<[f64; 4]>,
    /// The vendored solver's work count.
    pub work: u64,
}

impl Projection {
    /// Discover the states reachable from `initial` when both reactions have
    /// positive rates. The caller checks the documents with CRN-craft first.
    pub fn new(initial: [u32; SPECIES], times: &[f64], tolerance: f64) -> Self {
        let mut states = vec![initial.map(i64::from)];
        let mut index = HashMap::from([(states[0], 0)]);
        let mut infection = Vec::new();
        let mut recovery = Vec::new();
        let mut rows: Vec<Vec<(usize, Entry)>> = Vec::new();
        let mut i = 0;
        while i < states.len() {
            let [s, infected, r] = states[i];
            let mut targets = [None; REACTIONS];
            // Infection S + I -> 2I and recovery I -> R, in reaction order.
            let changes = [
                (s >= 1 && infected >= 1).then_some([s - 1, infected + 1, r]),
                (infected >= 1).then_some([s, infected - 1, r + 1]),
            ];
            for (target, next) in targets.iter_mut().zip(changes) {
                let Some(next) = next else { continue };
                let j = *index.entry(next).or_insert_with(|| {
                    states.push(next);
                    states.len() - 1
                });
                *target = Some(j);
            }
            infection.push(targets[0]);
            recovery.push(targets[1]);
            let mut row = vec![(i, Entry::Diagonal)];
            row.extend(targets[0].map(|j| (j, Entry::Infection)));
            row.extend(targets[1].map(|j| (j, Entry::Recovery)));
            rows.push(row);
            i += 1;
        }
        let n = states.len();
        let mut start = vec![0; n + 2];
        for row in &rows {
            for &(target, _) in row {
                start[target + 1] += 1;
            }
        }
        for c in 1..start.len() {
            start[c] += start[c - 1];
        }
        let mut fill = start.clone();
        let mut entries = vec![(0, Entry::Diagonal); start[n + 1]];
        for (source, row) in rows.iter().enumerate() {
            for &(target, kind) in row {
                entries[fill[target]] = (source, kind);
                fill[target] += 1;
            }
        }
        Self {
            susceptible: states.iter().map(|x| x[0] as f64).collect(),
            infected: states.iter().map(|x| x[1] as f64).collect(),
            infection,
            recovery,
            start,
            edges: entries.len() - n,
            entries,
            times: times.to_vec(),
            tolerance,
        }
    }

    /// Number of projected states.
    #[cfg(test)]
    pub fn states(&self) -> usize {
        self.infected.len()
    }

    /// Conservative compact-JSON size of the solver's full result, used to
    /// show that its output-size checks cannot fail for this structure.
    pub fn output_bound(&self) -> u128 {
        let (n, t) = (self.infected.len() as u128, self.times.len() as u128);
        // Fixed fields and hashes, each state's counts, and per snapshot its
        // keys plus at most 25 bytes per number.
        4096 + n * (SPECIES as u128 * 12 + 4) + t * (512 + 25 * (n + 1 + 2 * SPECIES as u128 + 6))
    }

    /// Solve one binding of the mass-action coefficients, or return `None`
    /// when the vendored solver must decide the outcome.
    pub fn solve(&self, infection: f64, gamma: f64, diagnostics: bool) -> Option<Solution> {
        if !infection.is_normal() || !gamma.is_normal() || infection < 0.0 || gamma < 0.0 {
            return None;
        }
        let n = self.infected.len();
        let mut infection_rate = vec![0.0; n];
        let mut recovery_rate = vec![0.0; n];
        let mut exits = Vec::with_capacity(n);
        let mut lambda: f64 = 0.0;
        for i in 0..n {
            // Mass-action pair and first-order kernels, with their products
            // in CRN-craft's order and their normal-range guards.
            let mut row = [(usize::MAX, 0.0); REACTIONS];
            if let Some(j) = self.infection[i] {
                let first = infection * self.susceptible[i];
                let rate = first * self.infected[i];
                if !first.is_normal() || !rate.is_normal() {
                    return None;
                }
                infection_rate[i] = rate;
                row[0] = (j, rate);
            }
            if let Some(j) = self.recovery[i] {
                let rate = gamma * self.infected[i];
                if !rate.is_normal() {
                    return None;
                }
                recovery_rate[i] = rate;
                row[1] = (j, rate);
            }
            // The solver sums each row in ascending target order.
            if row[1].0 < row[0].0 {
                row.swap(0, 1);
            }
            let exit: f64 = row
                .iter()
                .filter(|(j, _)| *j != usize::MAX)
                .map(|(_, rate)| rate)
                .sum();
            if !exit.is_finite() {
                return None;
            }
            lambda = lambda.max(exit);
            exits.push(exit);
        }
        let total_steps_f = self
            .times
            .windows(2)
            .map(|ts| ((ts[1] - ts[0]) * lambda / 16.0).ceil().max(1.0))
            .sum::<f64>();
        if !total_steps_f.is_finite() || total_steps_f > MAX_WORK as f64 {
            return None;
        }
        let tail_per_step = self.tolerance / total_steps_f;
        if tail_per_step < 1e-18 {
            return None;
        }
        let mut diagonal = exits;
        if lambda > 0.0 {
            for (i, d) in diagonal.iter_mut().enumerate() {
                *d = 1.0 - *d / lambda;
                for rate in [&mut infection_rate[i], &mut recovery_rate[i]] {
                    if *rate != 0.0 {
                        *rate /= lambda;
                        if *rate == 0.0 {
                            return None;
                        }
                    }
                }
            }
        }
        let columns: Vec<(usize, f64)> = self
            .entries
            .iter()
            .map(|&(source, kind)| {
                let weight = match kind {
                    Entry::Diagonal => diagonal[source],
                    Entry::Infection => infection_rate[source],
                    Entry::Recovery => recovery_rate[source],
                };
                (source, weight)
            })
            .collect();
        // Discovery visits every state and every reaction once.
        let mut work = (n * (1 + REACTIONS)) as u64;
        let mut term = vec![0.0; n + 1];
        let mut next = vec![0.0; n + 1];
        let mut result = vec![0.0; n + 1];
        let mut p = vec![0.0; n + 1];
        p[0] = 1.0;
        let mut means = Vec::with_capacity(self.times.len());
        let mut errors = Vec::with_capacity(if diagnostics { self.times.len() } else { 0 });
        let mut tail_bound = 0.0;
        let mut discarded = 0.0;
        let mut previous = 0.0;
        for &time in &self.times {
            let dt = time - previous;
            if lambda > 0.0 && dt > 0.0 {
                let steps = (dt * lambda / 16.0).ceil().max(1.0) as u64;
                let mu = lambda * (dt / steps as f64);
                for _ in 0..steps {
                    term.clone_from(&p);
                    let mut weight = libm::exp(-mu);
                    for (out, value) in result.iter_mut().zip(&term) {
                        *out = value * weight;
                    }
                    let mut k = 0;
                    loop {
                        work = work.saturating_add((n + self.edges) as u64);
                        if work > MAX_WORK {
                            return None;
                        }
                        let next_weight = weight * mu / (k + 1) as f64;
                        let ratio = mu / (k + 2) as f64;
                        if ratio < 1.0 && next_weight / (1.0 - ratio) <= tail_per_step {
                            tail_bound += next_weight / (1.0 - ratio);
                            discarded +=
                                omitted_weight(next_weight, mu, k + 1) * p.iter().sum::<f64>();
                            break;
                        }
                        multiply_accumulate(
                            &self.start,
                            &columns,
                            &term,
                            &mut next,
                            &mut result,
                            next_weight,
                        );
                        std::mem::swap(&mut term, &mut next);
                        weight = next_weight;
                        k += 1;
                    }
                    std::mem::swap(&mut p, &mut result);
                }
            }
            // Moments and checks visit every state once per species.
            work = work.saturating_add((n * (2 * SPECIES + 1)) as u64);
            if work > MAX_WORK {
                return None;
            }
            let mass: f64 = p.iter().sum();
            if p.iter().any(|x| !x.is_finite() || *x < 0.0) {
                return None;
            }
            means.push(self.infected.iter().zip(&p).map(|(x, q)| x * q).sum());
            if diagnostics {
                errors.push([p[n], tail_bound, 1.0 - mass, (1.0 - mass - discarded).abs()]);
            }
            previous = time;
        }
        Some(Solution {
            means,
            errors,
            work,
        })
    }
}

/// The vendored solver's omitted Poisson weight after term `first - 1`.
fn omitted_weight(first_weight: f64, mu: f64, first: u64) -> f64 {
    let mut term = first_weight;
    let mut sum = 0.0;
    let mut index = first;
    while term > 0.0 && term > sum * f64::EPSILON {
        sum += term;
        index += 1;
        term *= mu / index as f64;
    }
    sum
}

/// One series term, `next = term · P`, with the sink (index n) keeping its
/// mass, fused with `result += weight * next`, as in the vendored solver.
fn multiply_accumulate(
    start: &[usize],
    entries: &[(usize, f64)],
    term: &[f64],
    next: &mut [f64],
    result: &mut [f64],
    weight: f64,
) {
    let n = start.len() - 2;
    for (c, ((out, total), bounds)) in next
        .iter_mut()
        .zip(result.iter_mut())
        .zip(start.windows(2))
        .enumerate()
    {
        let mut sum = if c == n { term[n] } else { 0.0 };
        for &(i, rate) in &entries[bounds[0]..bounds[1]] {
            sum += term[i] * rate;
        }
        *out = sum;
        *total += weight * sum;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sir_documents;

    /// Compare one binding with the vendored solver bit for bit.
    fn compare(initial: [u32; 3], times: &[f64], tolerance: f64, beta: f64, gamma: f64) {
        let n = initial.iter().sum::<u32>();
        let infection = beta / f64::from(n);
        let (model, run) = sir_documents(initial, times, tolerance, infection, gamma, None);
        let control = crn_craft::simulation::Control::default();
        let reference = crn_craft::fsp::solve(model.as_bytes(), run.as_bytes(), &control);
        let projection = Projection::new(initial, times, tolerance);
        let fast = projection.solve(infection, gamma, true);
        let label = format!("{initial:?} {times:?} {tolerance} {beta} {gamma}");
        let Ok(out) = reference else {
            assert!(fast.is_none(), "{label}");
            return;
        };
        let fast = fast.unwrap_or_else(|| panic!("prepared path declined {label}"));
        assert_eq!(projection.states(), out.states.len(), "{label}");
        for (i, state) in out.states.iter().enumerate() {
            assert_eq!(projection.susceptible[i], state[0] as f64, "{label}");
            assert_eq!(projection.infected[i], state[1] as f64, "{label}");
        }
        assert_eq!(projection.edges, out.sparse_edges, "{label}");
        assert_eq!(fast.work, out.work, "{label}");
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let means: Vec<_> = out
            .snapshots
            .iter()
            .map(|s| s.projected_first_moment[1])
            .collect();
        assert_eq!(bits(&fast.means), bits(&means), "{label}");
        for (e, s) in fast.errors.iter().zip(&out.snapshots) {
            let expected = [
                s.physical_leakage_estimate,
                s.series_tail_upper_bound,
                s.mass_defect,
                s.roundoff_mass_residual,
            ];
            assert_eq!(bits(e), bits(&expected), "{label}");
        }
        assert_eq!(fast.errors.len(), out.snapshots.len());
        let lean = projection.solve(infection, gamma, false).unwrap();
        assert_eq!(bits(&lean.means), bits(&means), "{label}");
        assert!(lean.errors.is_empty());
        assert!(projection.output_bound() > serde_json::to_string(&out).unwrap().len() as u128);
    }

    #[test]
    fn prepared_projection_matches_the_vendored_solver_bit_for_bit() {
        let grids: [&[f64]; 4] = [
            &[0.0, 1.0, 2.0, 3.0, 5.0],
            &[0.0, 0.1],
            &[0.0, 1e-9, 0.5, 7.25, 20.0],
            &[0.0, 30.0, 120.0],
        ];
        let rates = [
            (0.7, 0.4),
            (5.0, 5.0),
            (1e-3, 2.0),
            (4.9, 1e-4),
            (0.3, 4.5),
            (1e-300, 1e-300),
        ];
        for initial in [
            [9, 2, 1],
            [1, 1, 0],
            [0, 3, 2],
            [5, 0, 5],
            [14, 1, 0],
            [3, 1, 2],
        ] {
            for times in grids {
                for tolerance in [1e-11, 1e-13, 1e-5] {
                    for &(beta, gamma) in &rates {
                        compare(initial, times, tolerance, beta, gamma);
                    }
                }
            }
        }
        // The largest supported population, including long multi-step intervals.
        compare([28, 2, 0], &[0.0, 30.0, 120.0], 1e-13, 5.0, 5.0);
        compare([20, 10, 0], &[0.0, 0.5, 1.0, 1.5, 2.0], 1e-12, 2.5, 0.3);
    }

    #[test]
    fn prepared_projection_declines_bindings_it_does_not_reproduce() {
        let projection = Projection::new([9, 2, 1], &[0.0, 1.0], 1e-11);
        for (infection, gamma) in [
            (0.0, 0.4),
            (0.1, 0.0),
            (f64::MIN_POSITIVE / 2.0, 0.4),
            (0.1, f64::NAN),
            (-0.1, 0.4),
        ] {
            assert!(projection.solve(infection, gamma, true).is_none());
        }
        assert!(projection.solve(0.1, 0.4, true).is_some());
    }
}
