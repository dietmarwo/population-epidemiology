mod domain;
mod sir_fsp;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json, value::RawValue};
use sha2::{Digest, Sha256};
use std::{
    collections::HashSet,
    path::PathBuf,
    sync::{
        OnceLock,
        atomic::{AtomicUsize, Ordering},
    },
    time::{Duration, Instant},
};

pub type Evaluation = Result<(f64, Value), Failure>;

#[derive(Debug, Clone)]
pub struct Failure {
    status: &'static str,
    detail: String,
    metrics: Option<Value>,
}
pub fn invalid(detail: impl Into<String>) -> Failure {
    Failure {
        status: "invalid",
        detail: detail.into(),
        metrics: None,
    }
}
pub fn infeasible(detail: impl Into<String>) -> Failure {
    Failure {
        status: "infeasible",
        detail: detail.into(),
        metrics: None,
    }
}

impl Failure {
    pub fn with_metrics(mut self, metrics: Value) -> Self {
        self.metrics = Some(metrics);
        self
    }
}

#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Candidate {
    pub id: String,
    pub parameters: Vec<f64>,
    #[serde(default)]
    pub trace: bool,
}

#[derive(Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Budget {
    workers: usize,
    max_steps: u64,
    wall_ms: Option<u64>,
    cancel_file: Option<PathBuf>,
}
impl Default for Budget {
    fn default() -> Self {
        Self {
            workers: 1,
            max_steps: 1_000_000,
            wall_ms: None,
            cancel_file: None,
        }
    }
}

pub struct Work<'a> {
    budget: &'a Budget,
    deadline: Option<Instant>,
    used: u64,
    monitored: bool,
    metrics_enabled: bool,
}
impl Work<'_> {
    pub fn metrics(&self, make: impl FnOnce() -> Value) -> Value {
        if self.metrics_enabled {
            make()
        } else {
            Value::Null
        }
    }

    #[inline]
    pub fn step(&mut self) -> Result<(), Failure> {
        if self.used < self.budget.max_steps && !self.monitored {
            self.used += 1;
            return Ok(());
        }
        self.checked_step()
    }

    /// Charge logical steps without bypassing monitoring or changing the count on failure.
    fn advance(&mut self, count: u64) -> Result<(), Failure> {
        if self.monitored {
            return (0..count).try_for_each(|_| self.step());
        }
        let available = self.budget.max_steps - self.used;
        self.used += count.min(available);
        if count > available {
            self.checked_step()
        } else {
            Ok(())
        }
    }

    #[cold]
    fn checked_step(&mut self) -> Result<(), Failure> {
        if self.budget.cancel_file.as_ref().is_some_and(|p| p.exists()) {
            return Err(Failure {
                status: "cancelled",
                detail: "cancel file exists".into(),
                metrics: None,
            });
        }
        if self.deadline.is_some_and(|d| Instant::now() >= d) {
            return Err(Failure {
                status: "wall_limit",
                detail: "batch deadline reached".into(),
                metrics: None,
            });
        }
        if self.used == self.budget.max_steps {
            return Err(Failure {
                status: "work_limit",
                detail: "candidate step budget reached".into(),
                metrics: None,
            });
        }
        self.used += 1;
        Ok(())
    }
}

/// One untraced candidate's outcome from `Prepared::score_matrix`.
#[derive(Clone, Debug, PartialEq)]
pub struct MatrixScore {
    /// Same strings as the JSON records.
    pub status: &'static str,
    pub steps: u64,
    /// `Some` only for `"complete"`.
    pub score: Option<f64>,
}

#[derive(Serialize)]
struct Record<'a> {
    id: &'a str,
    status: &'static str,
    steps: u64,
    score: Option<f64>,
    metrics: Option<Value>,
    detail: Option<String>,
}

#[derive(Serialize)]
struct Batch<'a> {
    model_sha256: &'a str,
    records: Vec<Record<'a>>,
    schema: &'static str,
    version: &'static str,
}

pub struct Prepared {
    model: domain::Model,
    model_sha256: String,
}
impl Prepared {
    pub fn new(model: &str) -> Result<Self, String> {
        let model_sha256 = format!("{:x}", Sha256::digest(model.as_bytes()));
        let model = domain::prepare(serde_json::from_str(model).map_err(|e| e.to_string())?)?;
        Ok(Self {
            model,
            model_sha256,
        })
    }

    pub fn evaluate_json(&self, candidates: &str, budget: &str) -> Result<String, String> {
        let candidates: Vec<Candidate> =
            serde_json::from_str(candidates).map_err(|e| e.to_string())?;
        let budget = serde_json::from_str(budget).map_err(|e| e.to_string())?;
        self.evaluate_encoded(candidates, budget, true)
    }

    /// Evaluate untraced candidates with full feasibility checks and null metrics.
    pub fn scores_json(&self, candidates: &str, budget: &str) -> Result<String, String> {
        let candidates: Vec<Candidate> =
            serde_json::from_str(candidates).map_err(|e| e.to_string())?;
        let budget = serde_json::from_str(budget).map_err(|e| e.to_string())?;
        if candidates.iter().any(|c| c.trace) {
            return Err("scores mode accepts untraced candidates only".into());
        }
        self.evaluate_encoded(candidates, budget, false)
    }

    /// Evaluate `rows` untraced candidates whose parameters are stored row-major in
    /// `parameters`. Each row has `parameters.len() / rows` values; dimension 0 is
    /// allowed. Budget JSON, limits, worker dispatch, deadlines, cancellation and
    /// failure classification are the same as `scores_json`.
    pub fn score_matrix(
        &self,
        rows: usize,
        parameters: &[f64],
        budget: &str,
    ) -> Result<Vec<MatrixScore>, String> {
        if !(1..=4096).contains(&rows) {
            return Err("matrix rows must be in 1..=4096".into());
        }
        if !parameters.len().is_multiple_of(rows) {
            return Err("parameter count is not a multiple of the row count".into());
        }
        let dimension = parameters.len() / rows;
        if dimension > 64 {
            return Err("parameter vector exceeds limit".into());
        }
        let budget: Budget = serde_json::from_str(budget).map_err(|e| e.to_string())?;
        check_budget(&budget, rows)?;
        let deadline = deadline(&budget);
        Ok(dispatch(rows, budget.workers, |index| {
            let row = &parameters[index * dimension..(index + 1) * dimension];
            let (result, steps) = self.outcome(row, false, &budget, deadline, false);
            match result {
                Ok((score, _)) => MatrixScore {
                    status: "complete",
                    steps,
                    score: Some(score),
                },
                Err(e) => MatrixScore {
                    status: e.status,
                    steps,
                    score: None,
                },
            }
        }))
    }

    /// SHA-256 of the exact model bytes passed to `Prepared::new`.
    pub fn model_sha256(&self) -> &str {
        &self.model_sha256
    }

    fn evaluate_encoded(
        &self,
        candidates: Vec<Candidate>,
        budget: Budget,
        metrics_enabled: bool,
    ) -> Result<String, String> {
        serde_json::to_string(&self.evaluate(&candidates, budget, metrics_enabled)?)
            .map_err(|e| e.to_string())
    }

    /// Run one candidate and return its checked result and logical step count.
    fn outcome(
        &self,
        parameters: &[f64],
        trace: bool,
        budget: &Budget,
        deadline: Option<Instant>,
        metrics_enabled: bool,
    ) -> (Evaluation, u64) {
        let mut work = Work {
            budget,
            deadline,
            used: 0,
            metrics_enabled,
            monitored: deadline.is_some() || budget.cancel_file.is_some(),
        };
        let result = domain::evaluate(&self.model, parameters, trace, &mut work);
        let result = result.and_then(|(score, metrics)| {
            if score.is_finite() {
                Ok((score, metrics))
            } else {
                Err(invalid("nonfinite objective"))
            }
        });
        (result, work.used)
    }

    fn record<'a>(
        &self,
        c: &'a Candidate,
        budget: &Budget,
        deadline: Option<Instant>,
        metrics_enabled: bool,
    ) -> Record<'a> {
        let (result, steps) =
            self.outcome(&c.parameters, c.trace, budget, deadline, metrics_enabled);
        match result {
            Ok((score, metrics)) => Record {
                id: &c.id,
                status: "complete",
                steps,
                score: Some(score),
                metrics: Some(metrics),
                detail: None,
            },
            Err(e) => Record {
                id: &c.id,
                status: e.status,
                steps,
                score: None,
                metrics: if metrics_enabled { e.metrics } else { None },
                detail: Some(e.detail),
            },
        }
    }

    fn evaluate<'a>(
        &'a self,
        candidates: &'a [Candidate],
        budget: Budget,
        metrics_enabled: bool,
    ) -> Result<Batch<'a>, String> {
        check_budget(&budget, candidates.len())?;
        let mut seen = HashSet::with_capacity(candidates.len());
        for c in candidates {
            if c.id.is_empty() || c.id.len() > 128 || !seen.insert(&c.id) || c.parameters.len() > 64
            {
                return Err(
                    "invalid or duplicate candidate id, or parameter vector exceeds limit".into(),
                );
            }
        }
        if candidates.iter().filter(|c| c.trace).count() > 16 {
            return Err("at most 16 traced candidates per batch".into());
        }
        let deadline = deadline(&budget);
        let records = dispatch(candidates.len(), budget.workers, |index| {
            self.record(&candidates[index], &budget, deadline, metrics_enabled)
        });
        Ok(Batch {
            model_sha256: &self.model_sha256,
            records,
            schema: concat!(env!("CARGO_PKG_NAME"), ".result/0.1"),
            version: env!("CARGO_PKG_VERSION"),
        })
    }
}

fn check_budget(budget: &Budget, count: usize) -> Result<(), String> {
    if !(1..=64).contains(&budget.workers)
        || !(1..=100_000_000).contains(&budget.max_steps)
        || budget.wall_ms.is_some_and(|n| n == 0 || n > 86_400_000)
        || count == 0
        || count > 4096
    {
        return Err("invalid budget or batch size".into());
    }
    Ok(())
}

fn deadline(budget: &Budget) -> Option<Instant> {
    budget
        .wall_ms
        .map(|n| Instant::now() + Duration::from_millis(n))
}

/// Evaluate `count` indexes inline or on scoped workers and return results in index order.
fn dispatch<T: Send + Sync>(
    count: usize,
    workers: usize,
    evaluate: impl Fn(usize) -> T + Sync,
) -> Vec<T> {
    if workers == 1 || count == 1 {
        return (0..count).map(evaluate).collect();
    }
    let next = AtomicUsize::new(0);
    let slots: Vec<_> = (0..count).map(|_| OnceLock::new()).collect();
    std::thread::scope(|scope| {
        for _ in 0..workers.min(count) {
            scope.spawn(|| {
                loop {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    if index >= count {
                        break;
                    }
                    let result = evaluate(index);
                    assert!(slots[index].set(result).is_ok(), "duplicate work index");
                }
            });
        }
    });
    slots
        .into_iter()
        .map(|slot| slot.into_inner().expect("missing evaluation"))
        .collect()
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema: String,
    model: Box<RawValue>,
    candidates: Vec<Candidate>,
    #[serde(default)]
    budget: Budget,
}

pub fn run(text: &str) -> Result<Value, String> {
    let req: Request = serde_json::from_str(text).map_err(|e| e.to_string())?;
    if req.schema != concat!(env!("CARGO_PKG_NAME"), ".request/0.1") {
        return Err("unsupported request schema".into());
    }
    let prepared = Prepared::new(req.model.get())?;
    serde_json::to_value(prepared.evaluate(&req.candidates, req.budget, true)?)
        .map_err(|e| e.to_string())
}

pub fn capabilities() -> Value {
    json!({"schema": concat!(env!("CARGO_PKG_NAME"), ".capabilities/0.1"),
           "version": env!("CARGO_PKG_VERSION"), "models": domain::CAPABILITIES,
           "batch": true, "worker_limit": 64})
}

#[cfg(feature = "python")]
mod python {
    use super::*;
    use pyo3::{exceptions::PyValueError, prelude::*};

    #[pyclass(name = "Prepared")]
    struct PyPrepared(Prepared);

    #[pymethods]
    impl PyPrepared {
        #[new]
        fn new(model: &str) -> PyResult<Self> {
            Prepared::new(model)
                .map(Self)
                .map_err(PyValueError::new_err)
        }
        fn scores(&self, py: Python<'_>, candidates: String, budget: String) -> PyResult<String> {
            py.detach(|| self.0.scores_json(&candidates, &budget))
                .map_err(PyValueError::new_err)
        }
        fn batch(&self, py: Python<'_>, candidates: String, budget: String) -> PyResult<String> {
            py.detach(|| self.0.evaluate_json(&candidates, &budget))
                .map_err(PyValueError::new_err)
        }
        /// Score untraced parameter rows. Returns scores (NaN unless complete),
        /// statuses and logical steps in row order.
        fn score_matrix(
            &self,
            py: Python<'_>,
            rows: Vec<Vec<f64>>,
            budget: String,
        ) -> PyResult<(Vec<f64>, Vec<&'static str>, Vec<u64>)> {
            let dimension = rows.first().map_or(0, Vec::len);
            if rows.iter().any(|row| row.len() != dimension) {
                return Err(PyValueError::new_err("ragged parameter rows"));
            }
            let count = rows.len();
            let flat = rows.concat();
            let rows = py
                .detach(|| self.0.score_matrix(count, &flat, &budget))
                .map_err(PyValueError::new_err)?;
            let mut scores = Vec::with_capacity(count);
            let mut statuses = Vec::with_capacity(count);
            let mut steps = Vec::with_capacity(count);
            for row in rows {
                scores.push(row.score.unwrap_or(f64::NAN));
                statuses.push(row.status);
                steps.push(row.steps);
            }
            Ok((scores, statuses, steps))
        }
    }

    #[pymodule]
    fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<PyPrepared>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> Value {
        serde_json::from_str(include_str!("../models/example.json")).unwrap()
    }
    #[test]
    fn prepared_workers_return_same_completed_records() {
        let example = fixture();
        let prepared = Prepared::new(&example["model"].to_string()).unwrap();
        let candidates = json!([
            {"id":"first","parameters":example["parameters"]},
            {"id":"second","parameters":example["parameters"]},
        ]);
        let serial = prepared
            .evaluate_json(&candidates.to_string(), r#"{"workers":1}"#)
            .unwrap();
        let parallel = prepared
            .evaluate_json(&candidates.to_string(), r#"{"workers":4}"#)
            .unwrap();
        assert_eq!(serial, parallel);
        assert!(serial.contains(r#""status":"complete""#));
    }
    #[test]
    fn scores_keep_order_status_steps_and_objective() {
        let example = fixture();
        let prepared = Prepared::new(&example["model"].to_string()).unwrap();
        let mut candidates = json!([
            {"id":"valid","parameters":example["parameters"]},
            {"id":"invalid","parameters": vec![1.0; 63]},
        ]);
        for budget in ["{}", r#"{"workers":4}"#, r#"{"max_steps":1}"#] {
            let mut full: Value = serde_json::from_str(
                &prepared
                    .evaluate_json(&candidates.to_string(), budget)
                    .unwrap(),
            )
            .unwrap();
            for row in full["records"].as_array_mut().unwrap() {
                row["metrics"] = Value::Null;
            }
            let scores: Value = serde_json::from_str(
                &prepared
                    .scores_json(&candidates.to_string(), budget)
                    .unwrap(),
            )
            .unwrap();
            assert_eq!(full, scores);
        }
        candidates[0]["trace"] = json!(true);
        assert!(prepared.scores_json(&candidates.to_string(), "{}").is_err());
    }

    /// Compare every matrix row with the `scores_json` record for the same vector.
    fn assert_matrix_matches_scores(model: &Value, rows: &[Vec<f64>], budget: &str) {
        let prepared = Prepared::new(&model.to_string()).unwrap();
        let candidates: Vec<_> = rows
            .iter()
            .enumerate()
            .map(|(i, row)| json!({"id": format!("row-{i}"), "parameters": row}))
            .collect();
        let lean: Value = serde_json::from_str(
            &prepared
                .scores_json(&Value::from(candidates).to_string(), budget)
                .unwrap(),
        )
        .unwrap();
        let flat = rows.concat();
        let matrix = prepared.score_matrix(rows.len(), &flat, budget).unwrap();
        let records = lean["records"].as_array().unwrap();
        assert_eq!(matrix.len(), records.len());
        for (row, record) in matrix.iter().zip(records) {
            assert_eq!(row.status, record["status"].as_str().unwrap(), "{budget}");
            assert_eq!(row.steps, record["steps"].as_u64().unwrap(), "{budget}");
            assert_eq!(
                row.score.map(f64::to_bits),
                record["score"].as_f64().map(f64::to_bits),
                "{budget}"
            );
            assert_eq!(row.score.is_some(), row.status == "complete");
        }
    }

    #[test]
    fn score_matrix_matches_scores_json_for_sir_rows_and_budgets() {
        let example = fixture();
        let base: Vec<f64> = serde_json::from_value(example["parameters"].clone()).unwrap();
        let rows = vec![
            base.clone(),
            vec![base[0] * 1.2, base[1] * 0.9],
            vec![7.0, 0.4],
            vec![-0.1, 0.4],
            vec![0.0, 0.4],
            vec![0.7, 0.0],
        ];
        let cancel = std::env::temp_dir().join(format!("pe-matrix-cancel-{}", std::process::id()));
        std::fs::write(&cancel, b"").unwrap();
        let cancel_budget = json!({"cancel_file": cancel}).to_string();
        for budget in [
            "{}",
            r#"{"workers":4}"#,
            r#"{"max_steps":1}"#,
            r#"{"max_steps":2}"#,
            r#"{"workers":3,"max_steps":1}"#,
            r#"{"wall_ms":60000}"#,
            &cancel_budget,
        ] {
            assert_matrix_matches_scores(&example["model"], &rows, budget);
        }
        std::fs::remove_file(&cancel).unwrap();
        // Wrong dimension is a per-row invalid status, as in the JSON path.
        assert_matrix_matches_scores(&example["model"], &[vec![0.5; 3], vec![0.5; 3]], "{}");
        let prepared = Prepared::new(&example["model"].to_string()).unwrap();
        let absent = json!({"cancel_file": cancel.with_extension("absent")}).to_string();
        let unaffected = prepared.score_matrix(1, &base, &absent).unwrap();
        assert_eq!(unaffected[0].status, "complete");
    }

    #[test]
    fn score_matrix_matches_scores_json_for_dimension_zero_models() {
        let moran = json!({"kind":"moran","population":6,"mutant_a":2,"mutant_b":1,"updates":5});
        let calendar = json!({"kind":"calendar","initial_on_probability":0.25,
            "segments":[{"duration":1.0,"on_rate":0.5,"off_rate":0.25},
                        {"duration":0.5,"on_rate":0.0,"off_rate":2.0}]});
        for model in [moran, calendar] {
            for budget in [
                "{}",
                r#"{"workers":2}"#,
                r#"{"max_steps":1}"#,
                r#"{"max_steps":2}"#,
            ] {
                assert_matrix_matches_scores(&model, &[vec![], vec![], vec![]], budget);
                assert_matrix_matches_scores(&model, &[vec![1.0], vec![2.0]], budget);
            }
        }
    }

    #[test]
    fn score_matrix_rejects_invalid_shapes_and_budgets() {
        let example = fixture();
        let prepared = Prepared::new(&example["model"].to_string()).unwrap();
        assert!(prepared.score_matrix(0, &[], "{}").is_err());
        assert!(prepared.score_matrix(4097, &[], "{}").is_err());
        assert!(prepared.score_matrix(4096, &[], "{}").is_ok());
        assert!(prepared.score_matrix(2, &[1.0, 2.0, 3.0], "{}").is_err());
        assert!(prepared.score_matrix(1, &[0.5; 65], "{}").is_err());
        assert_eq!(
            prepared.score_matrix(1, &[0.5; 64], "{}").unwrap()[0].status,
            "invalid"
        );
        for budget in [
            "",
            "null",
            r#"{"workers":0}"#,
            r#"{"workers":65}"#,
            r#"{"max_steps":0}"#,
            r#"{"wall_ms":0}"#,
            r#"{"unknown":1}"#,
        ] {
            assert!(
                prepared.score_matrix(1, &[0.7, 0.4], budget).is_err(),
                "{budget}"
            );
        }
        // Values that JSON cannot carry are classified per row.
        let rows = prepared
            .score_matrix(2, &[f64::NAN, 0.4, 0.7, f64::INFINITY], "{}")
            .unwrap();
        assert!(
            rows.iter()
                .all(|r| r.status == "invalid" && r.score.is_none())
        );
    }

    #[test]
    fn prepared_sir_records_equal_the_vendored_solver_records() {
        // A cancel-file budget whose file never exists routes every SIR
        // candidate through the vendored solver with identical documents.
        let absent = std::env::temp_dir().join(format!("pe-absent-{}", std::process::id()));
        let reference = json!({"cancel_file": absent}).to_string();
        let mut example = fixture();
        let models = [
            example["model"].clone(),
            json!({"kind":"sir","susceptible":4,"infected":1,"recovered":0,
                   "sample_times":[0.0,0.5,4.0],"observations":[1.0,null,0.5],
                   "train_count":2,"tolerance":1e-13}),
            json!({"kind":"sir","susceptible":3,"infected":0,"recovered":1,
                   "sample_times":[0.0,2.0],"observations":[0.0,0.0],
                   "train_count":1,"tolerance":1e-9}),
            // CRN-craft requires the save grid to start at 0.
            json!({"kind":"sir","susceptible":3,"infected":2,"recovered":1,
                   "sample_times":[1.0,2.0],"observations":[1.0,1.0],
                   "train_count":1,"tolerance":1e-9}),
        ];
        let candidates = json!([
            {"id":"fit","parameters":[0.7,0.4]},
            {"id":"fast","parameters":[5.0,5.0]},
            {"id":"zero-beta","parameters":[0.0,0.4]},
            {"id":"zero-gamma","parameters":[0.7,0.0]},
            {"id":"tiny","parameters":[1e-310,0.4]},
            {"id":"range","parameters":[5.5,0.4]},
        ])
        .to_string();
        for model in models {
            example["model"] = model;
            let prepared = Prepared::new(&example["model"].to_string()).unwrap();
            for metrics in [true, false] {
                let run = |budget: &str| {
                    if metrics {
                        prepared.evaluate_json(&candidates, budget).unwrap()
                    } else {
                        prepared.scores_json(&candidates, budget).unwrap()
                    }
                };
                assert_eq!(run("{}"), run(&reference));
                assert_eq!(run(r#"{"workers":3}"#), run(&reference));
            }
        }
        let rejected = Prepared::new(
            &json!({"kind":"sir","susceptible":3,"infected":2,"recovered":1,
                    "sample_times":[1.0,2.0],"observations":[1.0,1.0],
                    "train_count":1,"tolerance":1e-9})
            .to_string(),
        )
        .unwrap();
        let rows = rejected.score_matrix(1, &[0.7, 0.4], "{}").unwrap();
        assert_eq!(rows[0].status, "infeasible");
    }

    #[test]
    fn bulk_moran_charges_match_per_row_charges_at_every_limit() {
        let model = json!({"kind":"moran","population":4,"mutant_a":1,"mutant_b":2,"updates":6});
        let prepared = Prepared::new(&model.to_string()).unwrap();
        let candidates = json!([{"id":"m","parameters":[]}]).to_string();
        let absent = std::env::temp_dir().join(format!("pe-moran-absent-{}", std::process::id()));
        // 15 states, charged once for construction and once per update.
        let total = 15 * 7;
        for max_steps in [1, 14, 15, 16, 29, 30, 31, total - 1, total, total + 1] {
            let unmonitored = json!({"max_steps": max_steps}).to_string();
            let monitored = json!({"max_steps": max_steps, "cancel_file": absent}).to_string();
            let record = prepared.evaluate_json(&candidates, &unmonitored).unwrap();
            assert_eq!(
                record,
                prepared.evaluate_json(&candidates, &monitored).unwrap()
            );
            let row: Value = serde_json::from_str(&record).unwrap();
            let expected = if max_steps >= total {
                "complete"
            } else {
                "work_limit"
            };
            assert_eq!(row["records"][0]["status"], expected);
            assert_eq!(row["records"][0]["steps"], max_steps.min(total));
        }
    }

    #[test]
    fn model_and_candidate_validation_are_independent() {
        let mut example = fixture();
        example["model"]["unexpected"] = json!(true);
        assert!(Prepared::new(&example["model"].to_string()).is_err());
        let example = fixture();
        let prepared = Prepared::new(&example["model"].to_string()).unwrap();
        let duplicate = json!([
            {"id":"a","parameters":example["parameters"]},
            {"id":"a","parameters":example["parameters"]},
        ]);
        assert!(
            prepared
                .evaluate_json(&duplicate.to_string(), "{}")
                .is_err()
        );
    }
}
