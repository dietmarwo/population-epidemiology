use population_epidemiology::Prepared;
use serde_json::{Value, json};
use std::{hint::black_box, time::Instant};

fn main() {
    let mut args = std::env::args().skip(1);
    let example = std::fs::read_to_string(args.next().expect("example JSON")).unwrap();
    let example: Value = serde_json::from_str(&example).unwrap();
    let mode = args.next().expect("workload");
    let mut model = example["model"].clone();
    #[allow(unused_mut)]
    let mut parameters = example["parameters"].clone();
    let (count, workers, traced, repeats) = match mode.as_str() {
        "one" => (1, 1, false, 40),
        "serial" => (96, 1, false, 4),
        "parallel" => (96, 4, false, 4),
        "trace" => (12, 1, true, 6),
        "invalid" => (1024, 1, false, 6),
        "domain" => {
            model = json!({"kind":"moran","population":18,"mutant_a":4,"mutant_b":3,"updates":65});
            parameters = json!([]);
            (48, 1, false, 5)
        }
        "domain2" => {
            model = json!({"kind":"calendar","initial_on_probability":0.3,"segments":vec![json!({"duration":0.25,"on_rate":0.9,"off_rate":0.2}); 350]});
            parameters = json!([]);
            (48, 1, false, 5)
        }
        _ => panic!("unknown workload"),
    };
    let candidates: Vec<_> = (0..count)
        .map(|i| {
            let mut values = parameters.clone();
            if mode == "invalid" {
                values = json!(vec![1.0; 63]);
            } else if let Some(first) = values.as_array_mut().unwrap().first_mut() {
                *first = json!(first.as_f64().unwrap() * (1.0 + (i % 13) as f64 * 0.003));
            }
            json!({"id": format!("candidate-{i:04}"), "parameters": values, "trace": traced})
        })
        .collect();
    let candidates = serde_json::to_string(&candidates).unwrap();
    let budget = json!({"workers":workers}).to_string();
    let evaluator = Prepared::new(&model.to_string()).unwrap();
    let lean = std::env::var_os("CRAFT_SCORES").is_some();
    let evaluate = || {
        if lean {
            evaluator.scores_json(&candidates, &budget)
        } else {
            evaluator.evaluate_json(&candidates, &budget)
        }
    };
    let first = evaluate().unwrap();
    if args.next().as_deref() == Some("probe") {
        println!("{first}");
        return;
    }
    black_box(first);
    let start = Instant::now();
    for _ in 0..repeats {
        black_box(evaluate().unwrap());
    }
    println!("{}", start.elapsed().as_secs_f64() / f64::from(repeats));
}
