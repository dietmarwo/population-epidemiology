//! Versioned data contracts. Validation is a separate, mandatory step.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ModelVersion {
    #[serde(rename = "crn-craft.model/0.1")]
    V01,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum RunVersion {
    #[serde(rename = "crn-craft.run/0.1")]
    V01,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ResultVersion {
    #[serde(rename = "crn-craft.result/0.1")]
    V01,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Convention {
    StochasticCounts,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Model {
    pub schema: ModelVersion,
    pub name: String,
    pub convention: Convention,
    /// Counts per chosen concentration unit. Does not rescale SSA coefficients.
    #[schemars(schema_with = "positive_number")]
    pub omega: f64,
    pub parameters: BTreeMap<String, f64>,
    pub species: Vec<Species>,
    pub reactions: Vec<Reaction>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Species {
    pub name: String,
    #[schemars(range(min = 0))]
    pub initial_count: i64,
    /// A safety stop, never a reflecting boundary or an FSP truncation bound.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub population_limit: Option<i64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Reaction {
    pub name: String,
    pub reactants: BTreeMap<String, i64>,
    pub products: BTreeMap<String, i64>,
    pub rate: Rate,
}

/// A literal or the exact name of a model parameter; no string evaluation.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum Scalar {
    Literal(f64),
    Parameter(String),
}

/// Custom laws are complete propensities, not mass-action multipliers.
/// All laws require sufficient reactants before firing (implemented in M1/M2).
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Rate {
    MassAction {
        coefficient: Scalar,
    },
    HillActivation {
        regulator: String,
        vmax: Scalar,
        basal: Scalar,
        k: Scalar,
        n: Scalar,
    },
    HillRepression {
        regulator: String,
        vmax: Scalar,
        basal: Scalar,
        k: Scalar,
        n: Scalar,
    },
    MichaelisMenten {
        substrate: String,
        vmax: Scalar,
        km: Scalar,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelRef {
    /// Portable, relative to the run file's directory; no parent traversal.
    pub path: String,
    /// SHA-256 of the exact model input bytes, lowercase hexadecimal.
    #[schemars(regex(pattern = "^[0-9a-f]{64}$"))]
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    Ssa {
        aggregator: Aggregator,
    },
    Fsp {
        bounds: BTreeMap<String, i64>,
        max_states: u64,
        tolerance: f64,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Aggregator {
    Direct,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Run {
    pub schema: RunVersion,
    pub model: ModelRef,
    pub method: Method,
    pub time: TimeGrid,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensemble: Option<Ensemble>,
    pub outputs: Outputs,
    pub budget: Budget,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TimeGrid {
    #[schemars(schema_with = "positive_number")]
    pub end: f64,
    /// Explicit sorted grid including 0 and end; no accumulated step rounding.
    pub save_times: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Ensemble {
    #[schemars(range(min = 1))]
    pub replicates: u64,
    pub root_seed: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Statistic {
    Mean,
    Variance,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Outputs {
    pub statistics: Vec<Statistic>,
    pub keep_trajectories: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    /// Required for SSA; omitted for FSP, which has max_states.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_events_per_replicate: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_s: Option<f64>,
    #[schemars(range(min = 1, max = 1024))]
    pub threads: u32,
    #[schemars(range(min = 1))]
    pub max_output_bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    Completed,
    EventLimit,
    WallTimeLimit,
    PopulationLimit,
    OutputLimit,
    Cancelled,
    Failed,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Assessment {
    ExactEnsemble,
    FspApproximation,
    Incomplete,
    FailedChecks,
}

/// Reserved result envelope; M0 never fabricates a simulation result.
/// Full statistics/archive contracts will be added with a new version in M2/M3.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunResult {
    pub schema: ResultVersion,
    pub model_sha256: String,
    pub run_sha256: String,
    pub termination: Termination,
    pub assessment: Assessment,
    pub requested_replicates: u64,
    pub completed_replicates: u64,
    pub events: u64,
    pub evidence: Vec<Evidence>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Evidence {
    pub check: String,
    pub passed: bool,
    pub detail: String,
}

fn positive_number(generator: &mut schemars::r#gen::SchemaGenerator) -> schemars::schema::Schema {
    let mut schema = f64::json_schema(generator).into_object();
    schema.number().exclusive_minimum = Some(0.0);
    schema.into()
}
