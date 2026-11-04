//! Finite state projection (Munsky & Khammash 2006), row-oriented generator.
//! Positive sparse uniformization, split so lambda*dt <= 16. No normalization.
use crate::{
    CompiledNetwork, check_bytes,
    contracts::Method,
    diagnostic::Diagnostic,
    parse_model, parse_run, sha256,
    simulation::{Control, Poller},
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, time::Instant};

pub const MAX_STATES: usize = 100_000;
pub const MAX_EDGES: usize = 2_000_000;
pub const MAX_WORK: u64 = 200_000_000;
const MAX_BYTES: u128 = 256 * 1024 * 1024;
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum FspVersion {
    #[serde(rename = "crn-craft.fsp-result/0.1")]
    V01,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FspSnapshot {
    pub time: f64,
    /// Index order matches states, excludes the absorbing escaping-mass sink.
    pub probabilities: Vec<f64>,
    pub physical_leakage_estimate: f64,
    /// Unnormalized moments of the retained law, never of the full unbounded process.
    pub projected_first_moment: Vec<f64>,
    pub projected_second_raw_moment: Vec<f64>,
    /// Cumulative discarded positive Poisson weight, excluding roundoff.
    pub series_tail_upper_bound: f64,
    pub discarded_weight_estimate: f64,
    pub mass_defect: f64,
    /// Mass-balance residual only; not a rigorous floating-point error bound.
    pub roundoff_mass_residual: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FspResult {
    pub schema: FspVersion,
    pub model_sha256: String,
    pub run_sha256: String,
    pub species: Vec<String>,
    pub states: Vec<Vec<i64>>,
    pub snapshots: Vec<FspSnapshot>,
    pub uniformization_rate: f64,
    pub sparse_edges: usize,
    pub work: u64,
    pub numerical_error_certified: bool,
}
fn error(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E070",
        "/fsp",
        message,
        "Reduce the projection/work or revise bounds/tolerance. No boundary is reflected or normalized.",
    )
}
struct Work<'a> {
    poller: Poller<'a>,
    used: u64,
}
impl Work<'_> {
    fn add(&mut self, amount: usize) -> Result<(), Diagnostic> {
        self.used = self.used.saturating_add(amount as u64);
        if self.used > MAX_WORK {
            return Err(error("FSP work limit exceeded"));
        }
        if let Some(stop) = self.poller.poll() {
            return Err(error(format!("FSP interrupted: {stop:?}")));
        }
        Ok(())
    }
}

/// The Poisson weight omitted after term `first - 1`: the sum of the
/// weights from index `first` on, starting with `first_weight`, summed
/// explicitly rather than as `1 - kept weight`, which cancels to roundoff.
/// The series stops only once `mu / (first + 1) < 1`, so the terms decrease.
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

/// Compact JSON length plus the trailing newline the CLI writes.
fn serialized_len(result: &FspResult) -> u64 {
    struct Counter(u64);
    impl std::io::Write for Counter {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 += bytes.len() as u64;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut counter = Counter(1);
    serde_json::to_writer(&mut counter, result).expect("finite FSP result");
    counter.0
}

/// Normalized rows: the diagonal weight and the off-diagonal (target, weight)
/// edges of each projected state; `usize::MAX` targets the escaping sink.
type Rows = [(f64, Vec<(usize, f64)>)];

/// The same generator by column: the incoming (source, weight) entries of
/// each target, the sink last, each in ascending source order with the
/// diagonal at its own position. Gathering a column adds its contributions
/// in the order the row-oriented product scattered them.
struct Columns {
    start: Vec<usize>,
    entries: Vec<(usize, f64)>,
}

fn columns(rows: &Rows) -> Columns {
    let n = rows.len();
    let target = |j: usize| if j == usize::MAX { n } else { j };
    let mut start = vec![0; n + 2];
    for (i, (_, edges)) in rows.iter().enumerate() {
        start[i + 1] += 1;
        for &(j, _) in edges {
            start[target(j) + 1] += 1;
        }
    }
    for c in 1..start.len() {
        start[c] += start[c - 1];
    }
    let mut fill = start.clone();
    let mut entries = vec![(0, 0.0); start[n + 1]];
    for (i, (diagonal, edges)) in rows.iter().enumerate() {
        entries[fill[i]] = (i, *diagonal);
        fill[i] += 1;
        for &(j, weight) in edges {
            let c = target(j);
            entries[fill[c]] = (i, weight);
            fill[c] += 1;
        }
    }
    Columns { start, entries }
}

/// One series term, `next = term · P`, with the sink (index n) keeping its
/// mass, fused with the series accumulation `result += weight * next`. Each
/// entry of `next` is the same ascending-source sum as the row-oriented
/// product computed, so results are bitwise unchanged.
fn multiply_accumulate(
    columns: &Columns,
    term: &[f64],
    next: &mut [f64],
    result: &mut [f64],
    weight: f64,
) {
    let n = columns.start.len() - 2;
    for (c, ((out, total), bounds)) in next
        .iter_mut()
        .zip(result.iter_mut())
        .zip(columns.start.windows(2))
        .enumerate()
    {
        let mut sum = if c == n { term[n] } else { 0.0 };
        for &(i, rate) in &columns.entries[bounds[0]..bounds[1]] {
            sum += term[i] * rate;
        }
        *out = sum;
        *total += weight * sum;
    }
}

pub fn solve(
    model_bytes: &[u8],
    run_bytes: &[u8],
    control: &Control,
) -> Result<FspResult, Vec<Diagnostic>> {
    let checked = check_bytes(model_bytes, Some(run_bytes));
    if !checked.valid {
        return Err(checked.diagnostics);
    }
    let model = parse_model(model_bytes).map_err(|e| vec![e])?;
    let run = parse_run(run_bytes).map_err(|e| vec![e])?;
    let Method::Fsp {
        bounds,
        max_states,
        tolerance,
    } = &run.method
    else {
        return Err(vec![error("FSP requires an FSP run")]);
    };
    let network = CompiledNetwork::compile(&model)?;
    let mut work = Work {
        poller: Poller::new(control, Instant::now(), run.budget.wall_time_s),
        used: 0,
    };
    let action = (|| -> Result<FspResult, Diagnostic> {
        let bounds: Vec<_> = network.species().iter().map(|s| bounds[s]).collect();
        let capacity = (*max_states).min(MAX_STATES as u64) as usize;
        let mut states = vec![network.initial_counts().to_vec()];
        let mut index = BTreeMap::from([(states[0].clone(), 0)]);
        let mut rows = Vec::new();
        let mut edge_count = 0;
        let mut lambda: f64 = 0.0;
        let mut i = 0;
        while i < states.len() {
            work.add(1)?;
            let rates = network.propensities(&states[i])?;
            let mut row = BTreeMap::new();
            for (r, &rate) in rates.iter().enumerate() {
                work.add(1)?;
                if rate == 0.0 {
                    continue;
                }
                let mut next = states[i].clone();
                let mut escaped = false;
                for change in &network.reactions()[r].net_changes {
                    let value = i128::from(next[change.species]) + i128::from(change.coefficient);
                    if value < 0 {
                        return Err(error("negative reachable count"));
                    }
                    if value > i128::from(bounds[change.species]) {
                        escaped = true;
                        break;
                    }
                    next[change.species] = value as i64;
                }
                let target = if escaped {
                    usize::MAX
                } else if let Some(&j) = index.get(&next) {
                    j
                } else {
                    if states.len() >= capacity {
                        return Err(error(
                            "reachable-state cap exceeded; no incomplete generator is solved",
                        ));
                    }
                    network.validate_state(&next)?;
                    let j = states.len();
                    index.insert(next.clone(), j);
                    states.push(next);
                    j
                };
                // A physical null reaction contributes no generator exit or sink.
                if target != i {
                    *row.entry(target).or_insert(0.0) += rate;
                }
            }
            let exit: f64 = row.values().sum();
            if !exit.is_finite() {
                return Err(error("generator exit rate overflows"));
            }
            lambda = lambda.max(exit);
            edge_count += row.len();
            if edge_count > MAX_EDGES {
                return Err(error("FSP sparse-edge cap exceeded"));
            }
            // Retained snapshots: probabilities, two moment vectors and the record.
            let memory = states.len() as u128 * (bounds.len() as u128 * 32 + 128)
                + edge_count as u128 * 64
                + states.len() as u128 * (run.time.save_times.len() as u128 * 8 + 64)
                + run.time.save_times.len() as u128 * (bounds.len() as u128 * 16 + 256);
            if memory > MAX_BYTES {
                return Err(error("FSP workspace reservation exceeds 256 MiB"));
            }
            rows.push((exit, row.into_iter().collect::<Vec<_>>()));
            i += 1;
        }
        let n = states.len();
        // An early estimate; the serialized result is checked exactly at the end.
        let output_reserve = n as u128 * (bounds.len() as u128 * 24 + 32)
            + n as u128 * run.time.save_times.len() as u128 * 32
            + run.time.save_times.len() as u128 * (bounds.len() as u128 * 48 + 400)
            + 16384;
        if output_reserve > u128::from(run.budget.max_output_bytes) {
            return Err(error("FSP output reservation exceeds budget"));
        }
        let total_steps_f = run
            .time
            .save_times
            .windows(2)
            .map(|ts| ((ts[1] - ts[0]) * lambda / 16.0).ceil().max(1.0))
            .sum::<f64>();
        if !total_steps_f.is_finite() || total_steps_f > MAX_WORK as f64 {
            return Err(error("FSP uniformization step limit exceeded"));
        }
        let tail_per_step = *tolerance / total_steps_f;
        if tail_per_step < 1e-18 {
            return Err(error(
                "requested series tolerance is below the supported per-step range",
            ));
        }
        // Normalize each immutable row once; retain the exact per-entry division.
        if lambda > 0.0 {
            for (diagonal, edges) in &mut rows {
                *diagonal = 1.0 - *diagonal / lambda;
                for (_, weight) in edges {
                    *weight /= lambda;
                    if *weight == 0.0 {
                        return Err(error(
                            "positive FSP transition underflows during normalization",
                        ));
                    }
                }
            }
        }
        let columns = columns(&rows);
        drop(rows);
        let mut term = vec![0.0; n + 1];
        let mut next = vec![0.0; n + 1];
        let mut result = vec![0.0; n + 1];
        let mut p = vec![0.0; n + 1];
        p[0] = 1.0;
        let mut snapshots = Vec::new();
        let mut tail_bound = 0.0;
        let mut discarded = 0.0;
        let mut previous = 0.0;
        for &time in &run.time.save_times {
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
                        work.add(n + edge_count)?;
                        let next_weight = weight * mu / (k + 1) as f64;
                        let ratio = mu / (k + 2) as f64;
                        if ratio < 1.0 && next_weight / (1.0 - ratio) <= tail_per_step {
                            tail_bound += next_weight / (1.0 - ratio);
                            discarded +=
                                omitted_weight(next_weight, mu, k + 1) * p.iter().sum::<f64>();
                            break;
                        }
                        multiply_accumulate(&columns, &term, &mut next, &mut result, next_weight);
                        std::mem::swap(&mut term, &mut next);
                        weight = next_weight;
                        k += 1;
                    }
                    std::mem::swap(&mut p, &mut result);
                }
            }
            // Moments and checks visit every state once per species.
            work.add(n * (2 * bounds.len() + 1))?;
            let mass: f64 = p.iter().sum();
            if p.iter().any(|x| !x.is_finite() || *x < 0.0) {
                return Err(error("invalid probability after uniformization"));
            }
            snapshots.push(FspSnapshot {
                time,
                probabilities: p[..n].to_vec(),
                physical_leakage_estimate: p[n],
                projected_first_moment: (0..bounds.len())
                    .map(|j| states.iter().zip(&p).map(|(x, q)| x[j] as f64 * q).sum())
                    .collect(),
                projected_second_raw_moment: (0..bounds.len())
                    .map(|j| {
                        states
                            .iter()
                            .zip(&p)
                            .map(|(x, q)| (x[j] as f64).powi(2) * q)
                            .sum()
                    })
                    .collect(),
                series_tail_upper_bound: tail_bound,
                discarded_weight_estimate: discarded,
                mass_defect: 1.0 - mass,
                roundoff_mass_residual: (1.0 - mass - discarded).abs(),
            });
            previous = time;
        }
        let result = FspResult {
            schema: FspVersion::V01,
            model_sha256: sha256(model_bytes),
            run_sha256: sha256(run_bytes),
            species: network.species().to_vec(),
            states,
            snapshots,
            uniformization_rate: lambda,
            sparse_edges: edge_count,
            work: work.used,
            numerical_error_certified: false,
        };
        if serialized_len(&result) > run.budget.max_output_bytes {
            return Err(error("FSP result exceeds max_output_bytes"));
        }
        Ok(result)
    })();
    action.map_err(|e| vec![e])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The row-oriented scatter product used before P17 and P18.
    fn reference(rows: &Rows, term: &[f64], next: &mut [f64]) {
        let n = rows.len();
        next.fill(0.0);
        next[n] = term[n];
        for (i, (diagonal, edges)) in rows.iter().enumerate() {
            next[i] += term[i] * diagonal;
            for &(j, rate) in edges {
                next[if j == usize::MAX { n } else { j }] += term[i] * rate;
            }
        }
    }

    fn random_rows(n: usize, seed: &mut u64) -> Vec<(f64, Vec<(usize, f64)>)> {
        let mut uniform = || {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (*seed >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..n)
            .map(|_| {
                let mut edges: Vec<(usize, f64)> = (0..4)
                    .map(|_| {
                        let j = (uniform() * (n + 1) as f64) as usize;
                        (if j >= n { usize::MAX } else { j }, uniform() * 0.2)
                    })
                    .collect();
                edges.sort_by_key(|&(j, _)| j);
                edges.dedup_by_key(|&mut (j, _)| j);
                (uniform() * 0.2, edges)
            })
            .collect()
    }

    #[test]
    fn the_gathered_product_and_fused_sum_match_the_row_product_bit_for_bit() {
        let mut seed = 20261101;
        for n in [1, 2, 17, 300] {
            let rows = random_rows(n, &mut seed);
            let columns = columns(&rows);
            for density in [0.0, 0.1, 0.5, 1.0] {
                let mut term = vec![0.0; n + 1];
                for (i, value) in term.iter_mut().enumerate() {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let u = (seed >> 11) as f64 / (1u64 << 53) as f64;
                    if u < density {
                        // Include subnormal and tiny values.
                        *value = if i % 7 == 0 {
                            5e-324 * (i + 1) as f64
                        } else {
                            u * 1e-3
                        };
                    }
                }
                let weight = 0.37;
                let (mut fast, mut slow) = (vec![1.0; n + 1], vec![2.0; n + 1]);
                let mut fused: Vec<f64> = (0..=n).map(|i| i as f64 * 1e-4).collect();
                let mut separate = fused.clone();
                multiply_accumulate(&columns, &term, &mut fast, &mut fused, weight);
                reference(&rows, &term, &mut slow);
                for (out, value) in separate.iter_mut().zip(&slow) {
                    *out += weight * value;
                }
                let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                assert_eq!(bits(&fast), bits(&slow), "n {n}, density {density}");
                assert_eq!(bits(&fused), bits(&separate), "n {n}, density {density}");
            }
        }
    }
}
