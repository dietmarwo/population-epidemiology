use crate::{Candidate, Evaluation, Work, infeasible, invalid};
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
}
struct Moran {
    states: Vec<(usize, usize)>,
    offsets: Vec<usize>,
    transitions: Vec<(usize, f64)>,
    initial: usize,
}

pub fn prepare(input: Input) -> Result<Model, String> {
    let mut moran = None;
    let mut calendar = Vec::new();
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
    })
}

pub fn evaluate(model: &Model, c: &Candidate, work: &mut Work<'_>) -> Evaluation {
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
            let &[beta, gamma] = c.parameters.as_slice() else {
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
            let model = json!({"schema":"crn-craft.model/0.1","name":"small_sir","convention":"stochastic_counts","omega":1.0,
                "parameters":{"infection":beta/f64::from(n),"gamma":gamma},
                "species":[{"name":"S","initial_count":susceptible},{"name":"I","initial_count":infected},{"name":"R","initial_count":recovered}],
                "reactions":[{"name":"infection","reactants":{"S":1,"I":1},"products":{"I":2},
                    "rate":{"kind":"mass_action","coefficient":"infection"}},
                    {"name":"recovery","reactants":{"I":1},"products":{"R":1},
                    "rate":{"kind":"mass_action","coefficient":"gamma"}}]});
            let bytes = model.to_string();
            let wall = work.deadline.map(|d| {
                d.saturating_duration_since(std::time::Instant::now())
                    .as_secs_f64()
                    .max(1e-6)
            });
            let mut run = json!({"schema":"crn-craft.run/0.1","model":{"path":"sir.model.json","sha256":crn_craft::sha256(bytes.as_bytes())},
                "method":{"kind":"fsp","bounds":{"S":n,"I":n,"R":n},"max_states":(n+1)*(n+2)/2,"tolerance":tolerance},
                "time":{"end":sample_times.last(),"save_times":sample_times},
                "outputs":{"statistics":["mean","variance"],"keep_trajectories":0},
                "budget":{"threads":1,"max_output_bytes":8_388_608}});
            if let Some(wall) = wall {
                run["budget"]["wall_time_s"] = json!(wall);
            }
            let mut control = crn_craft::simulation::Control::default();
            control.cancel_file = work.budget.cancel_file.clone();
            let solved =
                crn_craft::fsp::solve(bytes.as_bytes(), run.to_string().as_bytes(), &control);
            work.step()?;
            let out =
                solved.map_err(|d| infeasible(serde_json::to_string(&d).unwrap_or_default()))?;
            let mut loss = [0.0; 2];
            let mut count = [0; 2];
            let mut means = Vec::new();
            let mut errors = Vec::new();
            for (i, s) in out.snapshots.iter().enumerate() {
                let mean = s.projected_first_moment[1];
                if work.metrics_enabled {
                    means.push(mean);
                    errors.push(json!({"physical_leakage":s.physical_leakage_estimate,
                    "series_tail_upper_bound":s.series_tail_upper_bound,"mass_defect":s.mass_defect,
                    "roundoff_mass_residual":s.roundoff_mass_residual}));
                }
                if let Some(o) = observations[i] {
                    let k = usize::from(i >= *train_count);
                    loss[k] += (mean - o).powi(2);
                    count[k] += 1;
                }
            }
            let train = loss[0] / f64::from(count[0]);
            let holdout = (count[1] > 0).then(|| loss[1] / f64::from(count[1]));
            let fsp_work = out.work;
            let model_hash = work.metrics_enabled.then(|| out.model_sha256.clone());
            let run_hash = work.metrics_enabled.then(|| out.run_sha256.clone());
            let certified = out.numerical_error_certified;
            Ok((
                train,
                work.metrics(|| json!({"infected_means":means,"holdout_mse":holdout,"error":errors,
                "fsp_work":fsp_work,"numerical_error_certified":certified,"crn_model_sha256":model_hash,
                "crn_run_sha256":run_hash,"fsp":c.trace.then_some(out)})),
            ))
        }
        Input::Moran {
            population,
            updates,
            ..
        } => {
            if !c.parameters.is_empty() {
                return Err(invalid("neutral Moran has no candidate parameters"));
            }
            moran(
                model.moran.as_ref().expect("prepared transitions"),
                *population,
                *updates,
                c.trace,
                work,
            )
        }
        Input::Calendar { .. } => {
            if !c.parameters.is_empty() {
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

fn prepare_moran(n: usize, a0: usize, b0: usize) -> Moran {
    let states: Vec<_> = (0..=n)
        .flat_map(|a| (0..=n - a).map(move |b| (a, b)))
        .collect();
    let index = |a: usize, b: usize| a * (n + 1) - a * a.saturating_sub(1) / 2 + b;
    let denom = (n * n) as f64;
    let mut offsets = Vec::with_capacity(states.len() + 1);
    let mut transitions = Vec::with_capacity(states.len() * 9);
    for &(a, b) in &states {
        let counts = [n - a - b, a, b];
        offsets.push(transitions.len());
        for parent in 0..3 {
            for death in 0..3 {
                if counts[parent] == 0 || counts[death] == 0 {
                    continue;
                }
                let mut next = counts;
                next[parent] += 1;
                next[death] -= 1;
                transitions.push((
                    index(next[1], next[2]),
                    (counts[parent] * counts[death]) as f64 / denom,
                ));
            }
        }
    }
    offsets.push(transitions.len());
    Moran {
        states,
        offsets,
        transitions,
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
        transitions,
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
        next.fill(0.0);
        for (&mass, offset) in law.iter().zip(offsets.windows(2)) {
            work.step()?;
            if mass == 0.0 {
                continue;
            }
            for &(to, p) in &transitions[offset[0]..offset[1]] {
                next[to] += mass * p;
            }
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
