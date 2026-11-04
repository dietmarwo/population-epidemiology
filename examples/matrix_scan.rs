//! Scan SIR rates through the typed score-matrix API.
//!
//! `cargo run --release --example matrix_scan -- models/example.json` scores
//! a 21 x 21 grid of (beta, gamma) rows without candidate or record JSON
//! and prints the best complete row, then replays it as a full JSON record.
use population_epidemiology::Prepared;
use serde_json::{Value, json};

fn main() -> Result<(), String> {
    let path = std::env::args()
        .nth(1)
        .ok_or("example JSON path required")?;
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let example: Value = serde_json::from_str(&text).map_err(|e| e.to_string())?;
    let prepared = Prepared::new(&example["model"].to_string())?;
    // Row-major parameters: two values (beta, gamma) per row.
    let mut parameters = Vec::new();
    for i in 0..21 {
        for j in 0..21 {
            parameters.extend([0.3 + 0.04 * f64::from(i), 0.1 + 0.03 * f64::from(j)]);
        }
    }
    let rows = parameters.len() / 2;
    let scores = prepared.score_matrix(rows, &parameters, r#"{"workers":2}"#)?;
    let (best, score) = scores
        .iter()
        .enumerate()
        .filter_map(|(i, row)| row.score.map(|s| (i, s)))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .ok_or("no complete row")?;
    let complete = scores.iter().filter(|r| r.status == "complete").count();
    let x = &parameters[2 * best..2 * best + 2];
    println!(
        "{complete} of {rows} rows complete; best beta {:.2} gamma {:.2} score {score:.3e}",
        x[0], x[1]
    );
    let replay =
        prepared.evaluate_json(&json!([{"id": "best", "parameters": x}]).to_string(), "{}")?;
    let record: Value = serde_json::from_str(&replay).map_err(|e| e.to_string())?;
    println!(
        "replay holdout MSE {}",
        record["records"][0]["metrics"]["holdout_mse"]
    );
    Ok(())
}
