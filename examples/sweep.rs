use serde_json::{Value, json};
use std::io;

fn main() -> Result<(), io::Error> {
    let file = std::env::args_os()
        .nth(1)
        .ok_or_else(|| io::Error::other("example.json path required"))?;
    let example: Value = serde_json::from_str(&std::fs::read_to_string(file)?)?;
    let parameters: Vec<f64> = serde_json::from_value(example["parameters"].clone())?;
    let prepared = population_epidemiology::Prepared::new(&example["model"].to_string())
        .map_err(io::Error::other)?;
    let candidates: Vec<_> = [0.75, 0.9, 1.0, 1.1, 1.25].iter().enumerate()
        .map(|(i, scale)| json!({"id":i.to_string(),"parameters":parameters.iter().map(|p|p*scale).collect::<Vec<_>>()}))
        .collect();
    let result: Value = serde_json::from_str(
        &prepared
            .evaluate_json(&serde_json::to_string(&candidates)?, r#"{"workers":2}"#)
            .map_err(io::Error::other)?,
    )?;
    let best = result["records"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|row| row["status"] == "complete")
        .min_by(|a, b| {
            a["score"]
                .as_f64()
                .unwrap()
                .total_cmp(&b["score"].as_f64().unwrap())
        })
        .ok_or_else(|| io::Error::other("no feasible candidate"))?;
    println!("{best}");
    Ok(())
}
