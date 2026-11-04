//! Versioned, JSON-friendly analysis and named-state input.
use crate::{
    conservation::{ConservationReport, analyze_conservation},
    diagnostic::Diagnostic,
    network::{CompiledNetwork, ReactionStructure},
    parse_model, sha256, strict_json,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum StateVersion {
    #[serde(rename = "crn-craft.state/0.1")]
    V01,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum AnalysisVersion {
    #[serde(rename = "crn-craft.analysis/0.1")]
    V01,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NamedState {
    pub schema: StateVersion,
    pub counts: BTreeMap<String, i64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NetworkAnalysis {
    pub species: Vec<String>,
    pub initial_counts: Vec<i64>,
    pub evaluation_counts: Vec<i64>,
    pub reactions: Vec<ReactionStructure>,
    /// Sorted propensity invalidations for each fired reaction index.
    pub dependencies: Vec<Vec<usize>>,
    pub propensities: Vec<f64>,
    pub conservation: ConservationReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnalysisReport {
    pub schema: AnalysisVersion,
    pub success: bool,
    pub simulation_available: bool,
    pub model_sha256: Option<String>,
    pub state_sha256: Option<String>,
    pub analysis: Option<NetworkAnalysis>,
    pub diagnostics: Vec<Diagnostic>,
}

impl AnalysisReport {
    pub fn failure(error: Diagnostic) -> Self {
        Self {
            schema: AnalysisVersion::V01,
            success: false,
            simulation_available: false,
            model_sha256: None,
            state_sha256: None,
            analysis: None,
            diagnostics: vec![error],
        }
    }
}

pub fn analyze_bytes(model_bytes: &[u8], state_bytes: Option<&[u8]>) -> AnalysisReport {
    let mut report = AnalysisReport {
        schema: AnalysisVersion::V01,
        success: false,
        simulation_available: false,
        model_sha256: Some(sha256(model_bytes)),
        state_sha256: state_bytes.map(sha256),
        analysis: None,
        diagnostics: Vec::new(),
    };
    let result = (|| {
        let model = parse_model(model_bytes).map_err(|d| vec![d])?;
        let network = CompiledNetwork::compile(&model)?;
        let state = if let Some(bytes) = state_bytes {
            let state: NamedState = strict_json::parse(bytes).map_err(|d| vec![d])?;
            network.named_state(&state.counts).map_err(|d| vec![d])?
        } else {
            network.initial_counts().to_vec()
        };
        let propensities = network.propensities(&state).map_err(|d| vec![d])?;
        let conservation = analyze_conservation(&network);
        let bounded_paths: Vec<String> = conservation
            .initial_class_upper_bounds
            .iter()
            .enumerate()
            .filter(|(_, bound)| bound.is_some())
            .map(|(i, _)| format!("/species/{i}"))
            .collect();
        report.diagnostics = network
            .warnings()
            .iter()
            .filter(|d| d.code != "CRN-W001" || !bounded_paths.contains(&d.path))
            .cloned()
            .collect();
        report
            .diagnostics
            .extend(conservation.diagnostics.iter().cloned());
        Ok(NetworkAnalysis {
            species: network.species().to_vec(),
            initial_counts: network.initial_counts().to_vec(),
            evaluation_counts: state,
            reactions: network.reactions().to_vec(),
            dependencies: network.dependencies().to_vec(),
            propensities,
            conservation,
        })
    })();
    match result {
        Ok(analysis) => {
            report.analysis = Some(analysis);
            report.success = true;
        }
        Err(diagnostics) => {
            report.diagnostics = diagnostics;
        }
    }
    report
}
