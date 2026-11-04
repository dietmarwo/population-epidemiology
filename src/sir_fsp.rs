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
use crn_craft::fsp::{FspSnapshot, MAX_WORK};
use std::collections::HashMap;

/// The output budget the SIR run document declares.
pub const MAX_OUTPUT_BYTES: u64 = 8_388_608;
/// Species in declaration order: S, I and R.
const SPECIES: usize = 3;
/// Reactions in declaration order: infection and recovery.
const REACTIONS: usize = 2;

/// Rate-independent SIR projection structure.
pub struct Projection {
    /// Counts per state, in CRN-craft discovery order.
    counts: Vec<[i64; SPECIES]>,
    /// Susceptible and infected counts per state as floating-point values.
    susceptible: Vec<f64>,
    infected: Vec<f64>,
    /// Target state of each state's infection and recovery transition.
    infection: Vec<Option<usize>>,
    recovery: Vec<Option<usize>>,
    /// Column-ordered generator: the incoming entries of each target state
    /// in ascending source order, with the diagonal at its own position,
    /// padded after the real entries to three slots. The escaping sink has
    /// no incoming edges here.
    sources: Vec<[usize; 3]>,
    kinds: Vec<[Entry; 3]>,
    /// The same generator as a stencil on the (S, I) grid, when every
    /// column's ascending source order fits it.
    grid: Option<Grid>,
    edges: usize,
    times: Vec<f64>,
    tolerance: f64,
}

/// Zero-padded (S, I) grid storage for the series vectors. Every transition
/// lowers 2S + I by one, so discovery numbers states level by level and, as
/// checked at construction, each column (S, I) adds the recovery source
/// (S, I+1), the infection source (S+1, I-1) and its diagonal in that order.
struct Grid {
    stride: usize,
    /// Storage index of each state in discovery order.
    cells: Vec<usize>,
    /// Storage index of the escaping sink, the last element.
    sink: usize,
    /// Rows of reachable states: first storage index, length and first
    /// position in row-major state order.
    rows: Vec<(usize, usize, usize)>,
    /// In row-major order: the discovery index of each cell's state and of
    /// its recovery and infection sources, if those transitions exist.
    region: Vec<[Option<usize>; 3]>,
}

#[derive(Clone, Copy)]
enum Entry {
    Diagonal,
    Infection,
    Recovery,
    /// A zero weight after the real entries of a column.
    Padding,
}

/// How much of a solution to record.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Detail {
    /// Infected means only.
    Means,
    /// Means and per-snapshot error diagnostics.
    Diagnostics,
    /// Diagnostics and the vendored solver's full trace.
    Trace,
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
    /// The rest of the vendored solver's result in trace mode.
    pub trace: Option<Trace>,
}

/// Rate-dependent parts of `crn_craft::fsp::FspResult` other than hashes.
pub struct Trace {
    pub states: Vec<Vec<i64>>,
    pub snapshots: Vec<FspSnapshot>,
    pub uniformization_rate: f64,
    pub sparse_edges: usize,
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
        // A column holds its diagonal plus at most one infection source
        // (S + 1, I - 1) and one recovery source (I + 1, R - 1).
        assert_eq!(start[n], start[n + 1], "no transition reaches the sink");
        let mut sources = Vec::with_capacity(n);
        let mut kinds = Vec::with_capacity(n);
        for (c, bounds) in start[..=n].windows(2).enumerate() {
            let column = &entries[bounds[0]..bounds[1]];
            assert!(column.len() <= 3, "SIR columns have at most three entries");
            let mut slot_sources = [c; 3];
            let mut slot_kinds = [Entry::Padding; 3];
            for (k, &(source, kind)) in column.iter().enumerate() {
                slot_sources[k] = source;
                slot_kinds[k] = kind;
            }
            sources.push(slot_sources);
            kinds.push(slot_kinds);
        }
        let grid = Grid::new(&states, &sources, &kinds);
        Self {
            susceptible: states.iter().map(|x| x[0] as f64).collect(),
            infected: states.iter().map(|x| x[1] as f64).collect(),
            infection,
            recovery,
            sources,
            kinds,
            grid,
            edges: entries.len() - n,
            times: times.to_vec(),
            tolerance,
            counts: states,
        }
    }

    /// Whether the vendored solver's rate-independent output and workspace
    /// checks pass for this structure, so they cannot fail for a binding.
    pub fn fits_budgets(&self) -> bool {
        let n = self.infected.len() as u128;
        let t = self.times.len() as u128;
        let species = SPECIES as u128;
        let reserve = n * (species * 24 + 32) + n * t * 32 + t * (species * 48 + 400) + 16384;
        let memory = n * (species * 32 + 128)
            + self.edges as u128 * 64
            + n * (t * 8 + 64)
            + t * (species * 16 + 256);
        reserve <= u128::from(MAX_OUTPUT_BYTES)
            && memory <= 256 * 1024 * 1024
            && self.output_bound() <= u128::from(MAX_OUTPUT_BYTES)
    }

    /// Number of projected states.
    #[cfg(test)]
    pub fn states(&self) -> usize {
        self.infected.len()
    }

    /// Whether the grid stencil replaces the slot kernel.
    #[cfg(test)]
    pub fn uses_grid(&self) -> bool {
        self.grid.is_some()
    }

    /// The same projection with the slot kernel.
    #[cfg(test)]
    fn without_grid(mut self) -> Self {
        self.grid = None;
        self
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
    pub fn solve(&self, infection: f64, gamma: f64, detail: Detail) -> Option<Solution> {
        let diagnostics = detail != Detail::Means;
        let traced = detail == Detail::Trace;
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
        let kernel = match &self.grid {
            Some(grid) => {
                let weight = |state: Option<usize>, rates: &[f64]| state.map_or(0.0, |i| rates[i]);
                let mut stencil = [Vec::new(), Vec::new(), Vec::new()];
                for &[state, recovery, infection] in &grid.region {
                    stencil[0].push(weight(recovery, &recovery_rate));
                    stencil[1].push(weight(infection, &infection_rate));
                    stencil[2].push(weight(state, &diagonal));
                }
                Kernel::Grid(grid, stencil)
            }
            None => Kernel::Slots(
                self.kinds
                    .iter()
                    .zip(&self.sources)
                    .map(|(kinds, sources)| {
                        [0, 1, 2].map(|k| match kinds[k] {
                            Entry::Diagonal => diagonal[sources[k]],
                            Entry::Infection => infection_rate[sources[k]],
                            Entry::Recovery => recovery_rate[sources[k]],
                            Entry::Padding => 0.0,
                        })
                    })
                    .collect(),
            ),
        };
        let identity: Vec<usize>;
        let (cells, sink, storage) = match &self.grid {
            Some(grid) => (grid.cells.as_slice(), grid.sink, grid.sink + 1),
            None => {
                identity = (0..n).collect();
                (identity.as_slice(), n, n + 1)
            }
        };
        // Sums over the solver's vector: states in discovery order, then the sink.
        let total = |v: &[f64]| -> f64 { cells.iter().map(|&c| v[c]).chain([v[sink]]).sum() };
        // Discovery visits every state and every reaction once.
        let mut work = (n * (1 + REACTIONS)) as u64;
        let mut term = vec![0.0; storage];
        let mut next = vec![0.0; storage];
        let mut result = vec![0.0; storage];
        let mut p = vec![0.0; storage];
        p[cells[0]] = 1.0;
        let mut means = Vec::with_capacity(self.times.len());
        let mut errors = Vec::with_capacity(if diagnostics { self.times.len() } else { 0 });
        let mut snapshots = Vec::with_capacity(if traced { self.times.len() } else { 0 });
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
                            discarded += omitted_weight(next_weight, mu, k + 1) * total(&p);
                            break;
                        }
                        match &kernel {
                            Kernel::Slots(weights) => multiply_accumulate(
                                &self.sources,
                                weights,
                                &term,
                                &mut next,
                                &mut result,
                                next_weight,
                            ),
                            Kernel::Grid(grid, stencil) => {
                                grid.multiply_accumulate(
                                    stencil,
                                    &term,
                                    &mut next,
                                    &mut result,
                                    next_weight,
                                );
                            }
                        }
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
            let mass = total(&p);
            if cells
                .iter()
                .chain([&sink])
                .any(|&c| !p[c].is_finite() || p[c] < 0.0)
            {
                return None;
            }
            means.push(
                self.infected
                    .iter()
                    .zip(cells)
                    .map(|(x, &c)| x * p[c])
                    .sum(),
            );
            if traced {
                // The vendored solver's snapshot, field by field.
                snapshots.push(FspSnapshot {
                    time,
                    probabilities: cells.iter().map(|&c| p[c]).collect(),
                    physical_leakage_estimate: p[sink],
                    projected_first_moment: (0..SPECIES)
                        .map(|j| {
                            self.counts
                                .iter()
                                .zip(cells)
                                .map(|(x, &c)| x[j] as f64 * p[c])
                                .sum()
                        })
                        .collect(),
                    projected_second_raw_moment: (0..SPECIES)
                        .map(|j| {
                            self.counts
                                .iter()
                                .zip(cells)
                                .map(|(x, &c)| (x[j] as f64).powi(2) * p[c])
                                .sum()
                        })
                        .collect(),
                    series_tail_upper_bound: tail_bound,
                    discarded_weight_estimate: discarded,
                    mass_defect: 1.0 - mass,
                    roundoff_mass_residual: (1.0 - mass - discarded).abs(),
                });
            }
            if diagnostics {
                errors.push([
                    p[sink],
                    tail_bound,
                    1.0 - mass,
                    (1.0 - mass - discarded).abs(),
                ]);
            }
            previous = time;
        }
        let trace = traced.then(|| Trace {
            states: self.counts.iter().map(|x| x.to_vec()).collect(),
            snapshots,
            uniformization_rate: lambda,
            sparse_edges: self.edges,
        });
        Some(Solution {
            means,
            errors,
            work,
            trace,
        })
    }
}

/// Per-binding weights of the selected product kernel.
enum Kernel<'a> {
    Slots(Vec<[f64; 3]>),
    /// Recovery, infection and diagonal weights in row-major cell order.
    Grid(&'a Grid, [Vec<f64>; 3]),
}

impl Grid {
    /// Build the grid if every S row of states is a contiguous I range and
    /// every column's ascending source order is recovery, infection, diagonal.
    fn new(
        states: &[[i64; SPECIES]],
        sources: &[[usize; 3]],
        kinds: &[[Entry; 3]],
    ) -> Option<Self> {
        let s_low = states.iter().map(|x| x[0]).min()?;
        let s_high = states.iter().map(|x| x[0]).max()?;
        let i_high = states.iter().map(|x| x[1]).max()?;
        let stride = usize::try_from(i_high).ok()? + 3;
        let height = usize::try_from(s_high - s_low).ok()? + 2;
        let cell = |x: &[i64; SPECIES]| (x[0] - s_low) as usize * stride + x[1] as usize + 1;
        let mut by_row: Vec<Vec<usize>> = vec![Vec::new(); height];
        for (i, x) in states.iter().enumerate() {
            by_row[(x[0] - s_low) as usize].push(i);
        }
        let mut rows = Vec::new();
        let mut region = Vec::with_capacity(states.len());
        for members in &mut by_row {
            members.sort_by_key(|&i| states[i][1]);
            let Some(&first) = members.first() else {
                continue;
            };
            let contiguous = members
                .iter()
                .enumerate()
                .all(|(k, &i)| states[i][1] == states[first][1] + k as i64);
            if !contiguous {
                return None;
            }
            rows.push((cell(&states[first]), members.len(), region.len()));
            for &i in members.iter() {
                let mut order = Vec::new();
                let mut slots = [Some(i), None, None];
                for (&source, &kind) in sources[i].iter().zip(&kinds[i]) {
                    let (position, expected) = match kind {
                        Entry::Recovery => (0, [states[i][0], states[i][1] + 1]),
                        Entry::Infection => (1, [states[i][0] + 1, states[i][1] - 1]),
                        Entry::Diagonal => (2, [states[i][0], states[i][1]]),
                        Entry::Padding => continue,
                    };
                    if states[source][..2] != expected {
                        return None;
                    }
                    if position < 2 {
                        slots[position + 1] = Some(source);
                    }
                    order.push(position);
                }
                if !order.windows(2).all(|w| w[0] < w[1]) || order.last() != Some(&2) {
                    return None;
                }
                region.push(slots);
            }
        }
        Some(Self {
            stride,
            cells: states.iter().map(cell).collect(),
            sink: height * stride,
            rows,
            region,
        })
    }

    /// `next = term · P` and `result += weight * next` on the grid, with the
    /// sink keeping its mass. Each cell adds its recovery source, infection
    /// source and diagonal in the vendored solver's order; a missing
    /// transition has weight 0 and adds +0.0 to a non-negative sum.
    fn multiply_accumulate(
        &self,
        [recovery, infection, diagonal]: &[Vec<f64>; 3],
        term: &[f64],
        next: &mut [f64],
        result: &mut [f64],
        weight: f64,
    ) {
        let w = self.stride;
        for &(o, len, k) in &self.rows {
            let right = &term[o + 1..o + 1 + len];
            let lower = &term[o + w - 1..o + w - 1 + len];
            let here = &term[o..o + len];
            let (wr, wi, wd) = (
                &recovery[k..k + len],
                &infection[k..k + len],
                &diagonal[k..k + len],
            );
            let out = &mut next[o..o + len];
            let total = &mut result[o..o + len];
            for i in 0..len {
                let mut sum = right[i] * wr[i];
                sum += lower[i] * wi[i];
                sum += here[i] * wd[i];
                out[i] = sum;
                total[i] += weight * sum;
            }
        }
        next[self.sink] = term[self.sink];
        result[self.sink] += weight * term[self.sink];
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
/// The vendored solver starts each column at 0.0 and adds its entries in
/// ascending source order; `0.0 + x` is `x` and a padded slot adds +0.0 to a
/// non-negative sum, so each column's value is bitwise unchanged.
fn multiply_accumulate(
    sources: &[[usize; 3]],
    weights: &[[f64; 3]],
    term: &[f64],
    next: &mut [f64],
    result: &mut [f64],
    weight: f64,
) {
    let n = sources.len();
    for (((out, total), source), w) in next[..n]
        .iter_mut()
        .zip(&mut result[..n])
        .zip(sources)
        .zip(weights)
    {
        let mut sum = term[source[0]] * w[0];
        sum += term[source[1]] * w[1];
        sum += term[source[2]] * w[2];
        *out = sum;
        *total += weight * sum;
    }
    next[n] = term[n];
    result[n] += weight * term[n];
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sir_documents;

    /// Compare one binding with the vendored solver bit for bit, with the
    /// grid stencil and with the slot kernel.
    fn compare(initial: [u32; 3], times: &[f64], tolerance: f64, beta: f64, gamma: f64) {
        compare_kernel(initial, times, tolerance, beta, gamma, true);
        compare_kernel(initial, times, tolerance, beta, gamma, false);
    }

    fn compare_kernel(
        initial: [u32; 3],
        times: &[f64],
        tolerance: f64,
        beta: f64,
        gamma: f64,
        grid: bool,
    ) {
        let n = initial.iter().sum::<u32>();
        let infection = beta / f64::from(n);
        let (model, run) = sir_documents(initial, times, tolerance, infection, gamma, None);
        let control = crn_craft::simulation::Control::default();
        let reference = crn_craft::fsp::solve(model.as_bytes(), run.as_bytes(), &control);
        let projection = Projection::new(initial, times, tolerance);
        let projection = if grid {
            projection
        } else {
            projection.without_grid()
        };
        let fast = projection.solve(infection, gamma, Detail::Diagnostics);
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
        let lean = projection.solve(infection, gamma, Detail::Means).unwrap();
        assert_eq!(bits(&lean.means), bits(&means), "{label}");
        assert!(lean.errors.is_empty() && lean.trace.is_none());
        // The traced binding reproduces the solver's serialized result.
        let traced = projection.solve(infection, gamma, Detail::Trace).unwrap();
        assert_eq!(bits(&traced.means), bits(&means), "{label}");
        let trace = traced.trace.unwrap();
        let rebuilt = crn_craft::fsp::FspResult {
            schema: crn_craft::fsp::FspVersion::V01,
            model_sha256: crn_craft::sha256(model.as_bytes()),
            run_sha256: crn_craft::sha256(run.as_bytes()),
            species: ["S", "I", "R"].map(String::from).to_vec(),
            states: trace.states,
            snapshots: trace.snapshots,
            uniformization_rate: trace.uniformization_rate,
            sparse_edges: trace.sparse_edges,
            work: traced.work,
            numerical_error_certified: false,
        };
        assert_eq!(
            serde_json::to_string(&rebuilt).unwrap(),
            serde_json::to_string(&out).unwrap(),
            "{label}"
        );
        assert!(projection.fits_budgets());
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
        for initial in [
            [9, 2, 1],
            [1, 1, 0],
            [0, 3, 2],
            [14, 1, 0],
            [28, 2, 0],
            [20, 10, 0],
        ] {
            assert!(Projection::new(initial, &[0.0, 1.0], 1e-11).uses_grid());
        }
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
            assert!(projection.solve(infection, gamma, Detail::Trace).is_none());
        }
        assert!(projection.solve(0.1, 0.4, Detail::Trace).is_some());
    }
}
