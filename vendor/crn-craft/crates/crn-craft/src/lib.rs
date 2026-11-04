#![forbid(unsafe_code)]
//! CRN-craft M2: strict inputs, network analysis and Direct SSA.
pub mod analysis;
pub mod analysis_v1;
pub mod archive;
pub mod codegen;
pub mod cone;
pub mod conservation;
pub mod contracts;
pub mod diagnostic;
pub mod export;
pub mod fsp;
pub mod kinetics;
pub mod linkage;
pub mod network;
#[doc(hidden)]
pub mod next_reaction;
pub mod parameters;
mod rates;
pub mod sampling;
pub mod simulation;
pub mod ssa;
pub mod strict_json;
pub mod tau;
pub mod validate;

pub use analysis::{AnalysisReport, analyze_bytes};
pub use network::CompiledNetwork;

use contracts::{Model, Run, RunResult};
use diagnostic::{Diagnostic, Severity};
use schemars::{JsonSchema, schema_for};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::io::Read;
use std::path::Path;

#[derive(Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckReport {
    pub schema: String,
    pub valid: bool,
    pub scope: String,
    pub execution_available: bool,
    pub model_sha256: Option<String>,
    pub run_sha256: Option<String>,
    pub diagnostics: Vec<Diagnostic>,
}

impl CheckReport {
    fn new() -> Self {
        Self {
            schema: "crn-craft.check/0.1".into(),
            valid: false,
            scope: "structural_and_domain_validation".into(),
            execution_available: false,
            model_sha256: None,
            run_sha256: None,
            diagnostics: Vec::new(),
        }
    }
    pub fn from_error(diagnostic: Diagnostic) -> Self {
        let mut report = Self::new();
        report.diagnostics.push(diagnostic);
        report
    }
    fn finish(&mut self) {
        self.valid = !self
            .diagnostics
            .iter()
            .any(|item| item.severity == Severity::Error);
    }
}

pub fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

pub fn parse_model(bytes: &[u8]) -> Result<Model, Diagnostic> {
    strict_json::parse(bytes)
}
pub fn parse_run(bytes: &[u8]) -> Result<Run, Diagnostic> {
    strict_json::parse(bytes)
}

pub fn analyze_files(model_path: &Path, state_path: Option<&Path>) -> AnalysisReport {
    let model = match read_document(model_path) {
        Ok(bytes) => bytes,
        Err(error) => return AnalysisReport::failure(error),
    };
    let state = match state_path.map(read_document).transpose() {
        Ok(bytes) => bytes,
        Err(error) => return AnalysisReport::failure(error),
    };
    analyze_bytes(&model, state.as_deref())
}

/// Legacy M1a output and warning semantics, retained for existing consumers.
pub fn analyze_legacy_files(
    model_path: &Path,
    state_path: Option<&Path>,
) -> analysis_v1::AnalysisReport {
    let model = match read_document(model_path) {
        Ok(bytes) => bytes,
        Err(error) => return analysis_v1::AnalysisReport::failure(error),
    };
    let state = match state_path.map(read_document).transpose() {
        Ok(bytes) => bytes,
        Err(error) => return analysis_v1::AnalysisReport::failure(error),
    };
    analysis_v1::analyze_bytes(&model, state.as_deref())
}

pub fn export_gillespy2_file(model_path: &Path) -> export::ExportReport {
    match read_document(model_path) {
        Ok(bytes) => export::export_gillespy2(&bytes),
        Err(error) => export::ExportReport::failure(error),
    }
}

/// Validate in-memory input bytes. Only check_files verifies model.path identity.
pub fn check_bytes(model_bytes: &[u8], run_bytes: Option<&[u8]>) -> CheckReport {
    let mut report = CheckReport::new();
    report.model_sha256 = Some(sha256(model_bytes));
    report.run_sha256 = run_bytes.map(sha256);
    let model = parse_model(model_bytes);
    if let Ok(model) = &model {
        report.diagnostics.extend(validate::validate_model(model));
    }
    if let Some(bytes) = run_bytes {
        match parse_run(bytes) {
            Ok(run) => {
                if let Ok(model) = &model {
                    report.diagnostics.extend(validate::validate_run(
                        model,
                        &run,
                        report.model_sha256.as_deref().unwrap_or_default(),
                    ));
                }
            }
            Err(error) => report.diagnostics.push(error),
        }
    }
    if let Err(error) = model {
        report.diagnostics.insert(0, error);
    }
    report.finish();
    report
}

fn read_document(path: &Path) -> Result<Vec<u8>, Diagnostic> {
    let io_error = |error: std::io::Error| {
        Diagnostic::error(
            "CRN-E030",
            "",
            format!("{}: {error}", path.display()),
            "Provide a readable regular JSON file.",
        )
    };
    if !std::fs::metadata(path).map_err(io_error)?.is_file() {
        return Err(Diagnostic::error(
            "CRN-E030",
            "",
            format!("{} is not a regular file", path.display()),
            "Provide a regular JSON file.",
        ));
    }
    let file = std::fs::File::open(path).map_err(io_error)?;
    if !file.metadata().map_err(io_error)?.is_file() {
        return Err(Diagnostic::error(
            "CRN-E030",
            "",
            format!("{} is not a regular file", path.display()),
            "Provide a regular JSON file.",
        ));
    }
    let mut bytes = Vec::new();
    file.take(strict_json::MAX_DOCUMENT_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(io_error)?;
    if bytes.len() > strict_json::MAX_DOCUMENT_BYTES {
        return Err(Diagnostic::error(
            "CRN-E004",
            "",
            "document exceeds 16 MiB",
            "Reduce the input file.",
        ));
    }
    Ok(bytes)
}

/// Resolve model.path relative to RUN and require it to denote the supplied MODEL.
/// A symlink is allowed only if its canonical target stays inside RUN's directory.
pub fn check_files(model_path: &Path, run_path: Option<&Path>) -> CheckReport {
    let model_bytes = match read_document(model_path) {
        Ok(x) => x,
        Err(e) => return CheckReport::from_error(e),
    };
    let run_bytes = match run_path.map(read_document).transpose() {
        Ok(x) => x,
        Err(e) => return CheckReport::from_error(e),
    };
    let mut report = check_bytes(&model_bytes, run_bytes.as_deref());
    if let (Some(path), Some(bytes)) = (run_path, run_bytes.as_deref()) {
        if let Ok(run) = parse_run(bytes) {
            if validate::valid_relative_path(&run.model.path) {
                let parent = path
                    .parent()
                    .filter(|x| !x.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."));
                let resolved = (
                    parent.canonicalize(),
                    parent.join(&run.model.path).canonicalize(),
                    model_path.canonicalize(),
                );
                let matches = match resolved {
                    (Ok(base), Ok(referenced), Ok(supplied)) => {
                        referenced.starts_with(&base) && referenced == supplied
                    }
                    _ => false,
                };
                if !matches {
                    report.diagnostics.push(Diagnostic::error("CRN-E022", "/model/path", "model reference does not resolve to the supplied model inside the run directory", "Place the model below the run directory and pass that same file to check."));
                }
            }
        }
    }
    report.finish();
    report
}

pub fn schema(kind: &str) -> Option<Value> {
    match kind {
        "model" => serde_json::to_value(schema_for!(Model)).ok(),
        "run" => serde_json::to_value(schema_for!(Run)).ok(),
        "result" => serde_json::to_value(schema_for!(simulation::SimulationResult)).ok(),
        "result-v1" => serde_json::to_value(schema_for!(RunResult)).ok(),
        "archive" => serde_json::to_value(schema_for!(archive::Manifest)).ok(),
        "archive-check" => serde_json::to_value(schema_for!(archive::ArchiveReport)).ok(),
        "check" => serde_json::to_value(schema_for!(CheckReport)).ok(),
        "state" => serde_json::to_value(schema_for!(analysis::NamedState)).ok(),
        "analysis" => serde_json::to_value(schema_for!(AnalysisReport)).ok(),
        "analysis-v1" => serde_json::to_value(schema_for!(analysis_v1::AnalysisReport)).ok(),
        "export" => serde_json::to_value(schema_for!(export::ExportReport)).ok(),
        "tau-request" => serde_json::to_value(schema_for!(tau::Request)).ok(),
        "tau-result" => serde_json::to_value(schema_for!(tau::TauResult)).ok(),
        "fsp-result" => serde_json::to_value(schema_for!(fsp::FspResult)).ok(),
        "ssa-run" => serde_json::to_value(schema_for!(ssa::Request)).ok(),
        "ssa-result" => serde_json::to_value(schema_for!(ssa::SsaResult)).ok(),
        _ => None,
    }
}

pub fn describe() -> Value {
    json!({
        "schema": "crn-craft.describe/0.8", "version": env!("CARGO_PKG_VERSION"), "milestone": "SSA-acceleration",
        "commands": ["describe", "check", "analyze", "analyze-v1", "export", "simulate", "verify", "simulate-cached", "fsp", "tau", "ssa"], "simulation_available": true,
        "aggregators": [{"name": "direct", "status": "implemented", "implementation_milestone": "M2"}],
        "methods": [{"name": "ssa", "status": "implemented"}, {"name": "fsp", "status": "implemented"}, {"name": "tau", "status": "experimental_approximation"}],
        "sampler": sampling::SAMPLER,
        "implementation_variants": ["direct", "dependency_cached_direct"],
        "alternative_exact_ssa": {"request":"crn-craft.ssa-run/0.3", "result":"crn-craft.ssa-result/0.3", "accepted_requests":["crn-craft.ssa-run/0.1","crn-craft.ssa-run/0.2","crn-craft.ssa-run/0.3","crn-craft.ssa-run/0.4"], "next_reaction_request":"crn-craft.ssa-run/0.4", "next_reaction_result":"crn-craft.ssa-result/0.4", "aggregators":["direct","cached_direct","direct_cr","sorting_direct","rssa_cr","next_reaction"], "waiting_times":["inverse","ziggurat"], "fast":{"default":true, "default_by_request":{"crn-craft.ssa-run/0.1":"unavailable","crn-craft.ssa-run/0.2":false,"crn-craft.ssa-run/0.3":true}, "checked_policy":"method.fast=false: ChaCha8 words, lost-increment and 53-bit resolution checks, stop polls before every random word", "emulates":"JumpProcesses.jl 10.0.1: Xoshiro256++, running totals without loss or resolution checks, batched stop polls"}, "archive_format":null},
        "optional_integrations": {"python":"bindings/python", "jobs":"integrations/jobs", "ode_bdf_sensitivities":"integrations/diffsol"},
        "rate_laws": ["mass_action", "hill_activation", "hill_repression", "michaelis_menten"],
        "reference_exporters": [{"target":"gillespy2", "version": export::GILLESPY2_VERSION, "subset":"elementary_mass_action", "maximum_reactant_order":2}],
        "rate_evaluation_available": true, "observables": ["python_trajectory_autocorrelation_spectrum", "python_interval_censored_first_passage"],
        "model_cards": ["docs/model-cards/rate-laws.md", "docs/model-cards/direct.md"],
        "limits": {"document_bytes": strict_json::MAX_DOCUMENT_BYTES, "json_depth": strict_json::MAX_JSON_DEPTH,
            "species": validate::MAX_SPECIES, "parameters": validate::MAX_PARAMETERS,
            "reactions": validate::MAX_REACTIONS, "save_times": validate::MAX_SAVE_TIMES,
            "dependency_edges": network::MAX_DEPENDENCY_EDGES,
            "conservation_species": conservation::MAX_CONSERVATION_SPECIES,
            "conservation_reactions": conservation::MAX_CONSERVATION_REACTIONS,
            "cone_species": cone::MAX_CONE_SPECIES, "cone_reactions": cone::MAX_CONE_REACTIONS,
            "cone_rays": cone::MAX_CONE_RAYS, "cone_work": cone::MAX_CONE_WORK,
            "simulation_replicates": simulation::MAX_REPLICATES, "simulation_workers": simulation::MAX_WORKERS,
            "simulation_workspace_bytes": simulation::MAX_WORKSPACE_BYTES, "archive_bytes": simulation::MAX_RESULT_BYTES},
        "schemas": {"model": schema("model"), "run": schema("run"), "result": schema("result"), "result-v1":schema("result-v1"), "archive":schema("archive"), "archive-check":schema("archive-check"), "check": schema("check"), "state": schema("state"), "analysis": schema("analysis"), "analysis-v1": schema("analysis-v1"), "export": schema("export")}
    })
}
