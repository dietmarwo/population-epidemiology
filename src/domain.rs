use crate::sir_fsp::{self, Projection, Solution};
use crate::{Evaluation, Work, infeasible, invalid};
use serde::Deserialize;
use serde_json::{Value, json};

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
    /// Prepared SIR projection, absent when CRN-craft rejects the documents.
    sir: Option<Projection>,
}
struct Moran {
    states: Vec<(usize, usize)>,
    /// Incoming transitions of each target state, as (source, probability)
    /// in ascending source order and then the source's parent/death order.
    offsets: Vec<usize>,
    incoming: Vec<(usize, f64)>,
    initial: usize,
}

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
                if projection.output_bound() <= u128::from(sir_fsp::MAX_OUTPUT_BYTES) {
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
    Ok(Model {
        input,
        moran,
        calendar,
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
            let prepared = model
                .sir
                .as_ref()
                .filter(|_| !work.monitored && !trace)
                .and_then(|p| p.solve(infection, gamma, work.metrics_enabled));
            let solved = match prepared {
                Some(solution) => {
                    let (model_sha256, run_sha256) = if work.metrics_enabled {
                        let (model, run) = sir_documents(
                            initial,
                            sample_times,
                            *tolerance,
                            infection,
                            gamma,
                            None,
                        );
                        (
                            Some(crn_craft::sha256(model.as_bytes())),
                            Some(crn_craft::sha256(run.as_bytes())),
                        )
                    } else {
                        (None, None)
                    };
                    Ok(Solved {
                        solution,
                        model_sha256,
                        run_sha256,
                        // The vendored solver never certifies its numerical error.
                        certified: false,
                        trace: None,
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
            Ok((
                train,
                work.metrics(|| {
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
                }),
            ))
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
                work.metrics(|| json!({"boundary_time_and_on_probability":boundaries})),
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
    let (states, rows) = moran_transitions(n);
    let index = |a: usize, b: usize| a * (n + 1) - a * a.saturating_sub(1) / 2 + b;
    // Transpose by target. Visiting sources in ascending order keeps every
    // target's additions in the order the scatter update performed them.
    let mut offsets = vec![0; states.len() + 1];
    for &(to, _) in rows.iter().flatten() {
        offsets[to + 1] += 1;
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    let mut fill = offsets.clone();
    let mut incoming = vec![(0, 0.0); offsets[states.len()]];
    for (from, row) in rows.iter().enumerate() {
        for &(to, p) in row {
            incoming[fill[to]] = (from, p);
            fill[to] += 1;
        }
    }
    Moran {
        states,
        offsets,
        incoming,
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
        offsets,
        incoming,
        initial,
    } = topology;
    // Charge the same logical row visits as the cold construction.
    for _ in states {
        work.step()?;
    }
    let mut law = vec![0.0; states.len()];
    let mut next = vec![0.0; states.len()];
    law[*initial] = 1.0;
    for _ in 0..updates {
        // One logical row visit per state, as in the scatter update. A
        // zero-mass source adds +0.0, which leaves the sum unchanged.
        for (out, offset) in next.iter_mut().zip(offsets.windows(2)) {
            work.step()?;
            let mut sum = 0.0;
            for &(from, p) in &incoming[offset[0]..offset[1]] {
                sum += law[from] * p;
            }
            *out = sum;
        }
        std::mem::swap(&mut law, &mut next);
    }
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
