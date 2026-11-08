//! Where a prepared call's time goes, measured through the public API.
//!
//! `cargo run --release --example breakdown -- models/example.json` prints
//! JSON with median microseconds per call over 41 alternating rounds. For each workload it times
//! candidate decoding alone, the typed matrix API (evaluation without
//! JSON), lean `scores_json` and full `evaluate_json`. The differences give
//! the JSON boundary and the cost of full-record metrics. It also times model
//! preparation, the CRN-craft solver fallback for one SIR candidate (a
//! cancel-file budget routes SIR through `crn_craft::fsp::solve`) and worker
//! scaling for the 96-candidate SIR batch.
use population_epidemiology::{Candidate, Prepared};
use serde_json::{Value, json};
use std::{hint::black_box, time::Instant};

/// Alternating measurement rounds per case; the median is reported.
const ROUNDS: usize = 41;

/// Median microseconds per call over alternating rounds of `repeats` calls.
fn medians(rounds: usize, repeats: u32, cases: &mut [&mut dyn FnMut()]) -> Vec<f64> {
    let mut samples = vec![Vec::with_capacity(rounds); cases.len()];
    for case in cases.iter_mut() {
        case();
    }
    for _ in 0..rounds {
        for (case, times) in cases.iter_mut().zip(&mut samples) {
            let start = Instant::now();
            for _ in 0..repeats {
                case();
            }
            times.push(start.elapsed().as_secs_f64() * 1e6 / f64::from(repeats));
        }
    }
    samples
        .into_iter()
        .map(|mut v| {
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        })
        .collect()
}

fn candidates(count: usize, parameters: &[f64]) -> (String, Vec<f64>) {
    let rows: Vec<Vec<f64>> = (0..count)
        .map(|i| {
            let mut row = parameters.to_vec();
            if let Some(first) = row.first_mut() {
                *first *= 1.0 + (i % 13) as f64 * 0.003;
            }
            row
        })
        .collect();
    let json: Vec<_> = rows
        .iter()
        .enumerate()
        .map(|(i, r)| json!({"id": format!("candidate-{i:04}"), "parameters": r}))
        .collect();
    (Value::from(json).to_string(), rows.concat())
}

fn phases(name: &str, model: &Value, count: usize, parameters: &[f64], repeats: u32) -> Value {
    let text = model.to_string();
    let prepare = medians(
        ROUNDS,
        repeats.max(5),
        &mut [&mut || {
            black_box(Prepared::new(&text).unwrap());
        }],
    )[0];
    let prepared = Prepared::new(&text).unwrap();
    let (batch, flat) = candidates(count, parameters);
    let budget = r#"{"workers":1}"#;
    let times = medians(
        ROUNDS,
        repeats,
        &mut [
            &mut || {
                black_box(serde_json::from_str::<Vec<Candidate>>(&batch).unwrap());
            },
            &mut || {
                black_box(prepared.score_matrix(count, &flat, budget).unwrap());
            },
            &mut || {
                black_box(prepared.scores_json(&batch, budget).unwrap());
            },
            &mut || {
                black_box(prepared.evaluate_json(&batch, budget).unwrap());
            },
        ],
    );
    let [decode, matrix, lean, full] = [times[0], times[1], times[2], times[3]];
    json!({
        "workload": name, "candidates": count, "prepare_us": prepare,
        "decode_us": decode, "matrix_us": matrix, "lean_us": lean, "full_us": full,
        "matrix_us_per_candidate": matrix / count as f64,
        "lean_us_per_candidate": lean / count as f64,
        "full_us_per_candidate": full / count as f64,
        "lean_json_boundary_us": lean - matrix,
        "encode_and_checks_us": lean - matrix - decode,
        "full_metrics_us": full - lean,
        "full_output_bytes": prepared.evaluate_json(&batch, budget).unwrap().len(),
    })
}

fn main() {
    let path = std::env::args().nth(1).expect("example JSON path");
    let example: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let sir = &example["model"];
    let fit: Vec<f64> = serde_json::from_value(example["parameters"].clone()).unwrap();
    let moran = json!({"kind":"moran","population":18,"mutant_a":4,"mutant_b":3,"updates":65});
    let calendar = json!({"kind":"calendar","initial_on_probability":0.3,
        "segments":vec![json!({"duration":0.25,"on_rate":0.9,"off_rate":0.2}); 350]});
    let workloads = vec![
        phases("one (SIR)", sir, 1, &fit, 200),
        phases("serial (SIR)", sir, 96, &fit, 4),
        phases("invalid (SIR, 63 values)", sir, 1024, &[1.0; 63], 4),
        phases("domain (Moran)", &moran, 48, &[], 4),
        phases("domain2 (calendar)", &calendar, 48, &[], 40),
    ];

    // One SIR candidate through the prepared projection and through the
    // CRN-craft solver, which a cancel-file budget selects.
    let prepared = Prepared::new(&sir.to_string()).unwrap();
    let (one, _) = candidates(1, &fit);
    let absent = std::env::temp_dir().join(format!("breakdown-absent-{}", std::process::id()));
    let crn_solver = json!({"cancel_file": absent}).to_string();
    let mut traced: Value = serde_json::from_str(&one).unwrap();
    traced[0]["trace"] = json!(true);
    let traced = traced.to_string();
    let fallback = medians(
        ROUNDS,
        100,
        &mut [
            &mut || {
                black_box(prepared.scores_json(&one, "{}").unwrap());
            },
            &mut || {
                black_box(prepared.scores_json(&one, &crn_solver).unwrap());
            },
            &mut || {
                black_box(prepared.evaluate_json(&traced, "{}").unwrap());
            },
            &mut || {
                black_box(prepared.evaluate_json(&traced, &crn_solver).unwrap());
            },
        ],
    );

    // Worker scaling for the 96-candidate SIR batch.
    let (batch, flat) = candidates(96, &fit);
    let mut scaling = Vec::new();
    for workers in [1, 2, 4] {
        let budget = json!({"workers": workers}).to_string();
        let t = medians(
            ROUNDS,
            4,
            &mut [
                &mut || {
                    black_box(prepared.score_matrix(96, &flat, &budget).unwrap());
                },
                &mut || {
                    black_box(prepared.scores_json(&batch, &budget).unwrap());
                },
                &mut || {
                    black_box(prepared.evaluate_json(&batch, &budget).unwrap());
                },
            ],
        );
        scaling
            .push(json!({"workers": workers, "matrix_us": t[0], "lean_us": t[1], "full_us": t[2]}));
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "unit": "median microseconds per call",
            "available_parallelism": std::thread::available_parallelism().map(usize::from).ok(),
            "workloads": workloads,
            "sir_one_candidate": {
                "prepared_lean_us": fallback[0], "crn_solver_lean_us": fallback[1],
                "prepared_traced_full_us": fallback[2], "crn_solver_traced_full_us": fallback[3],
            },
            "sir_96_worker_scaling": scaling,
        }))
        .unwrap()
    );
}
