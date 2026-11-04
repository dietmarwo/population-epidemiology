use crate::sir_fsp::{self, Detail, Projection, Solution};
use crate::{Evaluation, Work, infeasible, invalid};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json, value::RawValue};

pub const CAPABILITIES: &[&str] = &[
    "sir_via_crn_fsp",
    "finite_neutral_two_mutant_moran",
    "exact_two_state_rate_calendar",
];
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Input {
    Sir {
        susceptible: u32,
        infected: u32,
        recovered: u32,
        sample_times: Vec<f64>,
        observations: Vec<Option<f64>>,
        train_count: usize,
        tolerance: f64,
    },
    Moran {
        population: usize,
        mutant_a: usize,
        mutant_b: usize,
        updates: usize,
    },
    Calendar {
        initial_on_probability: f64,
        segments: Vec<Segment>,
    },
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Segment {
    duration: f64,
    on_rate: f64,
    off_rate: f64,
}
pub struct Model {
    input: Input,
    moran: Option<Moran>,
    calendar: Vec<[f64; 2]>,
    /// Calendar metrics depend only on the model, so they are serialized once.
    calendar_metrics: Option<Box<RawValue>>,
    /// Prepared SIR projection, absent when CRN-craft rejects the documents.
    sir: Option<Projection>,
}
struct Moran {
    states: Vec<(usize, usize)>,
    /// Row stride of the zero-padded (N+3) x (N+3) law grid. State (a, b)
    /// lives at `(a + 1) * stride + b + 1`.
    stride: usize,
    /// First triangular index of each row a.
    rows: Vec<usize>,
    /// Weights by slot and triangular target index. A target (a, b) gets
    /// mass from (a-1, b), (a-1, b+1), (a, b-1), three self-replacement
    /// entries, (a, b+1), (a+1, b-1) and (a+1, b), in that order, which is
    /// ascending source order. Missing transitions have weight 0.
    weights: [Vec<f64>; MORAN_SLOTS],
    initial: usize,
}

/// Neighbor and self-replacement slots of the Moran stencil.
const MORAN_SLOTS: usize = 9;

pub fn prepare(input: Input) -> Result<Model, String> {
    let mut moran = None;
    let mut calendar = Vec::new();
    let mut sir = None;
    match &input {
        Input::Sir {
            susceptible,
            infected,
            recovered,
            sample_times,
            observations,
            train_count,
            tolerance,
        } => {
            let n = u64::from(*susceptible) + u64::from(*infected) + u64::from(*recovered);
            let dimensions = (1..=30).contains(&n)
                && (2..=64).contains(&sample_times.len())
                && sample_times.len() == observations.len()
                && (1..sample_times.len()).contains(train_count);
            let times = sample_times.iter().all(|t| (0.0..=120.0).contains(t))
                && sample_times.windows(2).all(|w| w[0] < w[1]);
            if !(dimensions
                && times
                && observations
                    .iter()
                    .flatten()
                    .all(|o| (0.0..=n as f64).contains(o))
                && observations[..*train_count].iter().any(Option::is_some)
                && (1e-13..=1e-5).contains(tolerance))
            {
                return Err("invalid SIR population, time grid, observations or tolerance".into());
            }
            // Positive normal rates leave CRN-craft's validation independent of
            // their values, so one checked binding qualifies the structure.
            let initial = [*susceptible, *infected, *recovered];
            let (model, run) = sir_documents(initial, sample_times, *tolerance, 1.0, 1.0, None);
            if crn_craft::check_bytes(model.as_bytes(), Some(run.as_bytes())).valid {
                let projection = Projection::new(initial, sample_times, *tolerance);
                if projection.fits_budgets() {
                    sir = Some(projection);
                }
            }
        }
        Input::Moran {
            population,
            mutant_a,
            mutant_b,
            updates,
        } => {
            if !(2..=32).contains(population)
                || mutant_a
                    .checked_add(*mutant_b)
                    .is_none_or(|s| s > *population)
                || *updates > 10_000
            {
                return Err("invalid Moran state or update count".into());
            }
            moran = Some(prepare_moran(*population, *mutant_a, *mutant_b));
        }
        Input::Calendar {
            initial_on_probability,
            segments,
        } => {
            if !(0.0..=1.0).contains(initial_on_probability)
                || !(1..=1024).contains(&segments.len())
                || segments.iter().any(|s| {
                    !(1e-10..=120.0).contains(&s.duration)
                        || ![s.on_rate, s.off_rate]
                            .iter()
                            .all(|v| (0.0..=1000.0).contains(v))
                })
            {
                return Err("invalid calendar probability, duration or rates".into());
            }
            let mut probability = *initial_on_probability;
            let mut time = 0.0;
            calendar.reserve(segments.len());
            for s in segments {
                let rate = s.on_rate + s.off_rate;
                if rate > 0.0 {
                    let equilibrium = s.on_rate / rate;
                    probability += (equilibrium - probability) * (-(-rate * s.duration).exp_m1());
                }
                time += s.duration;
                calendar.push([time, probability]);
            }
        }
    }
    let calendar_metrics = (!calendar.is_empty()).then(|| {
        serde_json::value::to_raw_value(&json!({"boundary_time_and_on_probability":calendar}))
            .expect("finite calendar boundaries")
    });
    Ok(Model {
        input,
        moran,
        calendar,
        calendar_metrics,
        sir,
    })
}

/// The CRN-craft model and FSP run documents for one SIR rate binding.
pub(crate) fn sir_documents(
    [susceptible, infected, recovered]: [u32; 3],
    sample_times: &[f64],
    tolerance: f64,
    infection: f64,
    gamma: f64,
    wall: Option<f64>,
) -> (String, String) {
    let n = susceptible + infected + recovered;
    let model = json!({"schema":"crn-craft.model/0.1","name":"small_sir","convention":"stochastic_counts","omega":1.0,
        "parameters":{"infection":infection,"gamma":gamma},
        "species":[{"name":"S","initial_count":susceptible},{"name":"I","initial_count":infected},{"name":"R","initial_count":recovered}],
        "reactions":[{"name":"infection","reactants":{"S":1,"I":1},"products":{"I":2},
            "rate":{"kind":"mass_action","coefficient":"infection"}},
            {"name":"recovery","reactants":{"I":1},"products":{"R":1},
            "rate":{"kind":"mass_action","coefficient":"gamma"}}]});
    let bytes = model.to_string();
    let mut run = json!({"schema":"crn-craft.run/0.1","model":{"path":"sir.model.json","sha256":crn_craft::sha256(bytes.as_bytes())},
        "method":{"kind":"fsp","bounds":{"S":n,"I":n,"R":n},"max_states":(n+1)*(n+2)/2,"tolerance":tolerance},
        "time":{"end":sample_times.last(),"save_times":sample_times},
        "outputs":{"statistics":["mean","variance"],"keep_trajectories":0},
        "budget":{"threads":1,"max_output_bytes":sir_fsp::MAX_OUTPUT_BYTES}});
    if let Some(wall) = wall {
        run["budget"]["wall_time_s"] = json!(wall);
    }
    (bytes, run.to_string())
}

/// SIR metrics with fields in the sorted key order of a `Value` map.
#[derive(Serialize)]
struct SirMetrics<'a> {
    crn_model_sha256: Option<&'a str>,
    crn_run_sha256: Option<&'a str>,
    error: Vec<ErrorMetrics>,
    fsp: Option<SortedFsp<'a>>,
    fsp_work: u64,
    holdout_mse: Option<f64>,
    infected_means: &'a [f64],
    numerical_error_certified: bool,
}

#[derive(Serialize)]
struct ErrorMetrics {
    mass_defect: f64,
    physical_leakage: f64,
    roundoff_mass_residual: f64,
    series_tail_upper_bound: f64,
}

/// `crn_craft::fsp::FspResult` with fields in sorted key order.
#[derive(Serialize)]
struct SortedFsp<'a> {
    model_sha256: &'a str,
    numerical_error_certified: bool,
    run_sha256: &'a str,
    schema: &'a crn_craft::fsp::FspVersion,
    snapshots: Vec<SortedSnapshot<'a>>,
    sparse_edges: usize,
    species: &'a [String],
    states: &'a [Vec<i64>],
    uniformization_rate: f64,
    work: u64,
}

/// `crn_craft::fsp::FspSnapshot` with fields in sorted key order.
#[derive(Serialize)]
struct SortedSnapshot<'a> {
    discarded_weight_estimate: f64,
    mass_defect: f64,
    physical_leakage_estimate: f64,
    probabilities: &'a [f64],
    projected_first_moment: &'a [f64],
    projected_second_raw_moment: &'a [f64],
    roundoff_mass_residual: f64,
    series_tail_upper_bound: f64,
    time: f64,
}

/// Serialize SIR metrics to the same bytes as the equivalent `json!` value.
fn sir_metrics(solved: &Solved, holdout: Option<f64>) -> Box<RawValue> {
    let fsp = solved.trace.as_ref().map(|t| SortedFsp {
        model_sha256: &t.model_sha256,
        numerical_error_certified: t.numerical_error_certified,
        run_sha256: &t.run_sha256,
        schema: &t.schema,
        snapshots: t
            .snapshots
            .iter()
            .map(|s| SortedSnapshot {
                discarded_weight_estimate: s.discarded_weight_estimate,
                mass_defect: s.mass_defect,
                physical_leakage_estimate: s.physical_leakage_estimate,
                probabilities: &s.probabilities,
                projected_first_moment: &s.projected_first_moment,
                projected_second_raw_moment: &s.projected_second_raw_moment,
                roundoff_mass_residual: s.roundoff_mass_residual,
                series_tail_upper_bound: s.series_tail_upper_bound,
                time: s.time,
            })
            .collect(),
        sparse_edges: t.sparse_edges,
        species: &t.species,
        states: &t.states,
        uniformization_rate: t.uniformization_rate,
        work: t.work,
    });
    let metrics = SirMetrics {
        crn_model_sha256: solved.model_sha256.as_deref(),
        crn_run_sha256: solved.run_sha256.as_deref(),
        error: solved
            .solution
            .errors
            .iter()
            .map(|e| ErrorMetrics {
                physical_leakage: e[0],
                series_tail_upper_bound: e[1],
                mass_defect: e[2],
                roundoff_mass_residual: e[3],
            })
            .collect(),
        fsp,
        fsp_work: solved.solution.work,
        holdout_mse: holdout,
        infected_means: &solved.solution.means,
        numerical_error_certified: solved.certified,
    };
    serde_json::value::to_raw_value(&metrics).expect("finite SIR metrics")
}

/// Means, diagnostics, work and provenance of one solved SIR binding.
struct Solved {
    solution: Solution,
    model_sha256: Option<String>,
    run_sha256: Option<String>,
    certified: bool,
    trace: Option<crn_craft::fsp::FspResult>,
}

pub fn evaluate(model: &Model, parameters: &[f64], trace: bool, work: &mut Work<'_>) -> Evaluation {
    match &model.input {
        Input::Sir {
            susceptible,
            infected,
            recovered,
            sample_times,
            observations,
            train_count,
            tolerance,
        } => {
            let &[beta, gamma] = parameters else {
                return Err(invalid("expected beta and gamma per time unit"));
            };
            if ![beta, gamma].iter().all(|v| (0.0..=5.0).contains(v)) {
                return Err(invalid("rates must be in [0,5]"));
            }
            work.step()?;
            if work.used == work.budget.max_steps {
                work.step()?;
            }
            let n = susceptible + infected + recovered;
            let initial = [*susceptible, *infected, *recovered];
            let infection = beta / f64::from(n);
            let detail = if trace {
                Detail::Trace
            } else if work.metrics_enabled {
                Detail::Diagnostics
            } else {
                Detail::Means
            };
            let prepared = model
                .sir
                .as_ref()
                .filter(|_| !work.monitored)
                .and_then(|p| p.solve(infection, gamma, detail));
            let solved = match prepared {
                Some(mut solution) => {
                    let hashes = (work.metrics_enabled || trace).then(|| {
                        let (model_bytes, run_bytes) = sir_documents(
                            initial,
                            sample_times,
                            *tolerance,
                            infection,
                            gamma,
                            None,
                        );
                        (
                            crn_craft::sha256(model_bytes.as_bytes()),
                            crn_craft::sha256(run_bytes.as_bytes()),
                        )
                    });
                    let trace = solution.trace.take().map(|t| crn_craft::fsp::FspResult {
                        schema: crn_craft::fsp::FspVersion::V01,
                        model_sha256: hashes.as_ref().expect("trace hashes").0.clone(),
                        run_sha256: hashes.as_ref().expect("trace hashes").1.clone(),
                        species: ["S", "I", "R"].map(String::from).to_vec(),
                        states: t.states,
                        snapshots: t.snapshots,
                        uniformization_rate: t.uniformization_rate,
                        sparse_edges: t.sparse_edges,
                        work: solution.work,
                        numerical_error_certified: false,
                    });
                    let (model_sha256, run_sha256) = hashes.unzip();
                    Ok(Solved {
                        solution,
                        model_sha256,
                        run_sha256,
                        // The vendored solver never certifies its numerical error.
                        certified: false,
                        trace,
                    })
                }
                None => {
                    let wall = work.deadline.map(|d| {
                        d.saturating_duration_since(std::time::Instant::now())
                            .as_secs_f64()
                            .max(1e-6)
                    });
                    let (bytes, run) =
                        sir_documents(initial, sample_times, *tolerance, infection, gamma, wall);
                    let mut control = crn_craft::simulation::Control::default();
                    control.cancel_file = work.budget.cancel_file.clone();
                    crn_craft::fsp::solve(bytes.as_bytes(), run.as_bytes(), &control).map(|out| {
                        let solution = Solution {
                            means: out
                                .snapshots
                                .iter()
                                .map(|s| s.projected_first_moment[1])
                                .collect(),
                            errors: out
                                .snapshots
                                .iter()
                                .map(|s| {
                                    [
                                        s.physical_leakage_estimate,
                                        s.series_tail_upper_bound,
                                        s.mass_defect,
                                        s.roundoff_mass_residual,
                                    ]
                                })
                                .collect(),
                            work: out.work,
                            trace: None,
                        };
                        Solved {
                            solution,
                            model_sha256: work.metrics_enabled.then(|| out.model_sha256.clone()),
                            run_sha256: work.metrics_enabled.then(|| out.run_sha256.clone()),
                            certified: out.numerical_error_certified,
                            trace: trace.then_some(out),
                        }
                    })
                }
            };
            work.step()?;
            let solved =
                solved.map_err(|d| infeasible(serde_json::to_string(&d).unwrap_or_default()))?;
            let mut loss = [0.0; 2];
            let mut count = [0; 2];
            for (i, &mean) in solved.solution.means.iter().enumerate() {
                if let Some(o) = observations[i] {
                    let k = usize::from(i >= *train_count);
                    loss[k] += (mean - o).powi(2);
                    count[k] += 1;
                }
            }
            let train = loss[0] / f64::from(count[0]);
            let holdout = (count[1] > 0).then(|| loss[1] / f64::from(count[1]));
            Ok((train, work.raw_metrics(|| sir_metrics(&solved, holdout))))
        }
        Input::Moran {
            population,
            updates,
            ..
        } => {
            if !parameters.is_empty() {
                return Err(invalid("neutral Moran has no candidate parameters"));
            }
            moran(
                model.moran.as_ref().expect("prepared transitions"),
                *population,
                *updates,
                trace,
                work,
            )
        }
        Input::Calendar { .. } => {
            if !parameters.is_empty() {
                return Err(invalid("calendar rates come from the model"));
            }
            let boundaries = &model.calendar;
            work.advance(boundaries.len() as u64)?;
            Ok((
                boundaries.last().expect("nonempty calendar")[1],
                work.raw_metrics(|| {
                    model
                        .calendar_metrics
                        .clone()
                        .expect("prepared calendar metrics")
                }),
            ))
        }
    }
}

/// Outgoing (target, probability) transitions of one Moran state.
type Transitions = Vec<(usize, f64)>;

/// Triangular (a, b) states and their outgoing transitions in parent/death
/// order, with self replacement repeated.
fn moran_transitions(n: usize) -> (Vec<(usize, usize)>, Vec<Transitions>) {
    let states: Vec<_> = (0..=n)
        .flat_map(|a| (0..=n - a).map(move |b| (a, b)))
        .collect();
    let index = |a: usize, b: usize| a * (n + 1) - a * a.saturating_sub(1) / 2 + b;
    let denom = (n * n) as f64;
    let rows = states
        .iter()
        .map(|&(a, b)| {
            let counts = [n - a - b, a, b];
            let mut row = Vec::with_capacity(9);
            for parent in 0..3 {
                for death in 0..3 {
                    if counts[parent] == 0 || counts[death] == 0 {
                        continue;
                    }
                    let mut next = counts;
                    next[parent] += 1;
                    next[death] -= 1;
                    row.push((
                        index(next[1], next[2]),
                        (counts[parent] * counts[death]) as f64 / denom,
                    ));
                }
            }
            row
        })
        .collect();
    (states, rows)
}

fn prepare_moran(n: usize, a0: usize, b0: usize) -> Moran {
    let (states, transitions) = moran_transitions(n);
    let index = |a: usize, b: usize| a * (n + 1) - a * a.saturating_sub(1) / 2 + b;
    // Place every transition in its target's stencil slot. Visiting sources
    // in ascending order puts self-replacement entries in parent/death order.
    let mut weights: [Vec<f64>; MORAN_SLOTS] = std::array::from_fn(|_| vec![0.0; states.len()]);
    let mut selves = vec![0; states.len()];
    for (from, row) in transitions.iter().enumerate() {
        let (fa, fb) = states[from];
        for &(to, p) in row {
            let (ta, tb) = states[to];
            let slot = match (fa as isize - ta as isize, fb as isize - tb as isize) {
                (-1, 0) => 0,
                (-1, 1) => 1,
                (0, -1) => 2,
                (0, 0) => {
                    selves[to] += 1;
                    2 + selves[to]
                }
                (0, 1) => 6,
                (1, -1) => 7,
                (1, 0) => 8,
                _ => unreachable!("a Moran update moves one individual"),
            };
            assert!(weights[slot][to] == 0.0, "one transition per stencil slot");
            weights[slot][to] = p;
        }
    }
    Moran {
        states,
        stride: n + 3,
        rows: (0..=n).map(|a| index(a, 0)).collect(),
        weights,
        initial: index(a0, b0),
    }
}

fn moran(
    topology: &Moran,
    n: usize,
    updates: usize,
    trace: bool,
    work: &mut Work<'_>,
) -> Evaluation {
    let Moran {
        states,
        stride,
        rows,
        weights,
        initial,
    } = topology;
    let stride = *stride;
    // Charge the same logical row visits as the cold construction and one
    // row visit per state and update, as in the scatter update. Without a
    // deadline or cancel file, the work limit is the only possible stop, so
    // charging every visit up front gives the same status and step count.
    let monitored = work.monitored;
    if monitored {
        for _ in states {
            work.step()?;
        }
    } else {
        work.advance(states.len() as u64 * (updates as u64 + 1))?;
    }
    let cell = |a: usize, b: usize| (a + 1) * stride + b + 1;
    let mut law = vec![0.0; stride * stride];
    let mut next = vec![0.0; stride * stride];
    let (a0, b0) = states[*initial];
    law[cell(a0, b0)] = 1.0;
    for _ in 0..updates {
        for (a, &first) in rows.iter().enumerate() {
            let len = n + 1 - a;
            if monitored {
                for _ in 0..len {
                    work.step()?;
                }
            }
            let o = cell(a, 0);
            // Neighbor rows around (a, b) for b in 0..len; padding cells are 0.
            let up = &law[o - stride..o - stride + len + 1];
            let mid = &law[o - 1..o + len + 1];
            let down = &law[o + stride - 1..o + stride + len];
            let w: [&[f64]; MORAN_SLOTS] = std::array::from_fn(|k| &weights[k][first..first + len]);
            let out = &mut next[o..o + len];
            for b in 0..len {
                // The gather started at 0.0; 0.0 + x is x for x >= 0, and a
                // zero-weight slot adds +0.0 to a non-negative sum.
                let mut sum = up[b] * w[0][b];
                sum += up[b + 1] * w[1][b];
                sum += mid[b] * w[2][b];
                sum += mid[b + 1] * w[3][b];
                sum += mid[b + 1] * w[4][b];
                sum += mid[b + 1] * w[5][b];
                sum += mid[b + 2] * w[6][b];
                sum += down[b] * w[7][b];
                sum += down[b + 1] * w[8][b];
                out[b] = sum;
            }
        }
        std::mem::swap(&mut law, &mut next);
    }
    let law: Vec<f64> = states.iter().map(|&(a, b)| law[cell(a, b)]).collect();
    let mut a = vec![0.0; n + 1];
    let mut b = vec![0.0; n + 1];
    for (&(x, y), &mass) in states.iter().zip(&law) {
        a[x] += mass;
        b[y] += mass;
    }
    let mass: f64 = law.iter().sum();
    let joint: Option<Value> = trace.then(|| json!({"states":states,"probabilities":law}));
    Ok((
        a[n] + b[n],
        work.metrics(
            || json!({"marginal_a":a,"marginal_b":b,"mass":mass,"joint":joint,"updates":updates}),
        ),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The `json!` construction used through round 6.
    fn sir_metrics_value(solved: &Solved, holdout: Option<f64>) -> Value {
        let errors: Vec<_> = solved
            .solution
            .errors
            .iter()
            .map(|e| {
                json!({"physical_leakage":e[0],"series_tail_upper_bound":e[1],
                       "mass_defect":e[2],"roundoff_mass_residual":e[3]})
            })
            .collect();
        json!({"infected_means":solved.solution.means,"holdout_mse":holdout,"error":errors,
        "fsp_work":solved.solution.work,"numerical_error_certified":solved.certified,
        "crn_model_sha256":solved.model_sha256,"crn_run_sha256":solved.run_sha256,
        "fsp":solved.trace})
    }

    #[test]
    fn typed_sir_metrics_serialize_like_the_json_value() {
        let initial = [9, 2, 1];
        let times = [0.0, 1.0, 2.0, 3.0, 5.0];
        let projection = Projection::new(initial, &times, 1e-11);
        for (beta, gamma) in [(0.7, 0.4), (5.0, 0.01), (0.001, 4.0)] {
            let infection = beta / 12.0;
            let (model, run) = sir_documents(initial, &times, 1e-11, infection, gamma, None);
            let control = crn_craft::simulation::Control::default();
            let out = crn_craft::fsp::solve(model.as_bytes(), run.as_bytes(), &control).unwrap();
            for traced in [false, true] {
                let solution = projection
                    .solve(infection, gamma, Detail::Diagnostics)
                    .unwrap();
                let solved = Solved {
                    solution,
                    model_sha256: Some(out.model_sha256.clone()),
                    run_sha256: (!traced).then(|| out.run_sha256.clone()),
                    certified: traced,
                    trace: traced.then(|| out.clone()),
                };
                for holdout in [None, Some(0.125), Some(1e-300)] {
                    assert_eq!(
                        sir_metrics(&solved, holdout).get(),
                        sir_metrics_value(&solved, holdout).to_string()
                    );
                }
            }
        }
    }

    /// The scatter update used through 0.2.0.
    fn scatter(n: usize, a0: usize, b0: usize, updates: usize) -> Vec<f64> {
        let (states, rows) = moran_transitions(n);
        let index = |a: usize, b: usize| a * (n + 1) - a * a.saturating_sub(1) / 2 + b;
        let mut law = vec![0.0; states.len()];
        let mut next = vec![0.0; states.len()];
        law[index(a0, b0)] = 1.0;
        for _ in 0..updates {
            next.fill(0.0);
            for (&mass, row) in law.iter().zip(&rows) {
                if mass == 0.0 {
                    continue;
                }
                for &(to, p) in row {
                    next[to] += mass * p;
                }
            }
            std::mem::swap(&mut law, &mut next);
        }
        law
    }

    #[test]
    fn gathered_moran_update_matches_the_scatter_update_bit_for_bit() {
        let budget = crate::Budget::default();
        for (n, a0, b0, updates) in [
            (2, 1, 1, 0),
            (2, 0, 0, 3),
            (3, 1, 1, 3),
            (7, 2, 4, 40),
            (18, 4, 3, 65),
            (32, 0, 32, 5),
            (32, 10, 11, 300),
            (5, 5, 0, 30),
            (6, 0, 0, 12),
            (11, 3, 8, 77),
            (24, 1, 1, 150),
        ] {
            let mut work = crate::Work {
                budget: &budget,
                deadline: None,
                used: 0,
                monitored: false,
                metrics_enabled: true,
            };
            let (_, metrics) =
                moran(&prepare_moran(n, a0, b0), n, updates, true, &mut work).unwrap();
            let metrics: Value = serde_json::from_str(metrics.unwrap().get()).unwrap();
            let law: Vec<f64> =
                serde_json::from_value(metrics["joint"]["probabilities"].clone()).unwrap();
            let reference = scatter(n, a0, b0, updates);
            let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(&law), bits(&reference), "{n} {a0} {b0} {updates}");
            let states = (n + 1) * (n + 2) / 2;
            assert_eq!(work.used, (states * (updates + 1)) as u64);
        }
    }
}
