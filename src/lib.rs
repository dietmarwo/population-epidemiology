mod domain;

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

    fn evaluate_encoded(
        &self,
        candidates: Vec<Candidate>,
        budget: Budget,
        metrics_enabled: bool,
    ) -> Result<String, String> {
        serde_json::to_string(&self.evaluate(&candidates, budget, metrics_enabled)?)
            .map_err(|e| e.to_string())
    }

    fn record<'a>(
        &self,
        c: &'a Candidate,
        budget: &Budget,
        deadline: Option<Instant>,
        metrics_enabled: bool,
    ) -> Record<'a> {
        let mut work = Work {
            budget,
            deadline,
            used: 0,
            metrics_enabled,
            monitored: deadline.is_some() || budget.cancel_file.is_some(),
        };
        let result = domain::evaluate(&self.model, c, &mut work);
        let result = result.and_then(|(score, metrics)| {
            if score.is_finite() {
                Ok((score, metrics))
            } else {
                Err(invalid("nonfinite objective"))
            }
        });
        match result {
            Ok((score, metrics)) => Record {
                id: &c.id,
                status: "complete",
                steps: work.used,
                score: Some(score),
                metrics: Some(metrics),
                detail: None,
            },
            Err(e) => Record {
                id: &c.id,
                status: e.status,
                steps: work.used,
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
        if !(1..=64).contains(&budget.workers)
            || !(1..=100_000_000).contains(&budget.max_steps)
            || budget.wall_ms.is_some_and(|n| n == 0 || n > 86_400_000)
            || candidates.is_empty()
            || candidates.len() > 4096
        {
            return Err("invalid budget or batch size".into());
        }
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
        let deadline = budget
            .wall_ms
            .map(|n| Instant::now() + Duration::from_millis(n));
        let records: Vec<_> = if budget.workers == 1 || candidates.len() == 1 {
            candidates
                .iter()
                .map(|c| self.record(c, &budget, deadline, metrics_enabled))
                .collect()
        } else {
            let next = AtomicUsize::new(0);
            let slots: Vec<_> = (0..candidates.len()).map(|_| OnceLock::new()).collect();
            std::thread::scope(|scope| {
                for _ in 0..budget.workers.min(candidates.len()) {
                    scope.spawn(|| {
                        loop {
                            let index = next.fetch_add(1, Ordering::Relaxed);
                            let Some(c) = candidates.get(index) else {
                                break;
                            };
                            let record = self.record(c, &budget, deadline, metrics_enabled);
                            assert!(slots[index].set(record).is_ok(), "duplicate work index");
                        }
                    });
                }
            });
            slots
                .into_iter()
                .map(|slot| slot.into_inner().expect("missing evaluation"))
                .collect()
        };
        Ok(Batch {
            model_sha256: &self.model_sha256,
            records,
            schema: concat!(env!("CARGO_PKG_NAME"), ".result/0.1"),
            version: env!("CARGO_PKG_VERSION"),
        })
    }
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
