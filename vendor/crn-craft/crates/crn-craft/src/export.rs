//! Explicit elementary stochastic-count exporter. No target runtime is linked.
use crate::{
    contracts::{Rate, Scalar},
    diagnostic::{Diagnostic, Severity},
    parse_model, sha256,
    validate::validate_model,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const GILLESPY2_VERSION: &str = "1.8.3";
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ExportVersion {
    #[serde(rename = "crn-craft.export/0.1")]
    V01,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PythonExport {
    pub target: String,
    pub target_version: String,
    pub source: String,
    pub species_aliases: BTreeMap<String, String>,
    pub reaction_aliases: BTreeMap<String, String>,
    pub notes: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExportReport {
    pub schema: ExportVersion,
    pub success: bool,
    pub simulation_available: bool,
    pub model_sha256: Option<String>,
    pub artifact: Option<PythonExport>,
    pub diagnostics: Vec<Diagnostic>,
}
impl ExportReport {
    pub fn failure(error: Diagnostic) -> Self {
        Self {
            schema: ExportVersion::V01,
            success: false,
            simulation_available: false,
            model_sha256: None,
            artifact: None,
            diagnostics: vec![error],
        }
    }
}
fn unsupported(path: &str, message: &str) -> Diagnostic {
    Diagnostic::error(
        "CRN-E050",
        path,
        message,
        "This exporter accepts elementary mass action without population limits; keep unsupported semantics in the source model.",
    )
}
pub fn export_gillespy2(model_bytes: &[u8]) -> ExportReport {
    let mut report = ExportReport {
        schema: ExportVersion::V01,
        success: false,
        simulation_available: false,
        model_sha256: Some(sha256(model_bytes)),
        artifact: None,
        diagnostics: Vec::new(),
    };
    let model = match parse_model(model_bytes) {
        Ok(model) => model,
        Err(error) => {
            report.diagnostics.push(error);
            return report;
        }
    };
    report.diagnostics = validate_model(&model);
    if report
        .diagnostics
        .iter()
        .any(|d| d.severity == Severity::Error)
    {
        return report;
    }
    for (i, species) in model.species.iter().enumerate() {
        if species.population_limit.is_some() {
            report.diagnostics.push(unsupported(
                &format!("/species/{i}/population_limit"),
                "GillesPy2 cannot preserve the CRN-craft population safety stop",
            ));
        }
        if species.initial_count > 9_007_199_254_740_992 {
            report.diagnostics.push(unsupported(
                &format!("/species/{i}/initial_count"),
                "initial counts above 2^53 cannot be represented exactly by the target solver",
            ));
        }
    }
    let mut coefficients = Vec::new();
    for (i, reaction) in model.reactions.iter().enumerate() {
        let Rate::MassAction { coefficient } = &reaction.rate else {
            report.diagnostics.push(unsupported(
                &format!("/reactions/{i}/rate"),
                "custom propensity export is not implemented",
            ));
            continue;
        };
        if reaction.reactants.values().sum::<i64>() > 2 {
            report.diagnostics.push(unsupported(
                &format!("/reactions/{i}/reactants"),
                "GillesPy2 elementary mass action supports total reactant order at most two",
            ));
            continue;
        }
        let coefficient = match coefficient {
            Scalar::Literal(x) => *x,
            Scalar::Parameter(name) => model.parameters[name],
        };
        // GillesPy2's homodimer uses k*x*(x-1)/volume, without the /2 factorial.
        let translated = if reaction.reactants.values().any(|&v| v == 2) {
            coefficient / 2.0
        } else {
            coefficient
        };
        if reaction.reactants.values().any(|&v| v == 2) && translated * 2.0 != coefficient {
            report.diagnostics.push(Diagnostic::error("CRN-E051", &format!("/reactions/{i}/rate"),
                "homodimer coefficient cannot be halved exactly in f64",
                "Rescale the reference fixture; do not silently round or disable the exported channel."));
        }
        coefficients.push(translated);
    }
    if report
        .diagnostics
        .iter()
        .any(|d| d.severity == Severity::Error)
    {
        return report;
    }
    let aliases: BTreeMap<_, _> = model
        .species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.clone(), format!("S{i}")))
        .collect();
    let indices: BTreeMap<_, _> = model
        .species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();
    let mut source = format!(
        "# Generated by CRN-craft. Stochastic-count semantics; reference SSA only.\n# Model SHA-256: {}\nimport gillespy2 as gp\n\ndef build_model():\n    if gp.__version__ != {GILLESPY2_VERSION:?}:\n        raise RuntimeError('Requires GillesPy2 {GILLESPY2_VERSION}')\n    model = gp.Model(name={:?}, population=True, volume=1.0)\n    species = [\n",
        sha256(model_bytes),
        model.name
    );
    for (i, s) in model.species.iter().enumerate() {
        source.push_str(&format!(
            "        gp.Species(name='S{i}', initial_value={}, mode='discrete'),\n",
            s.initial_count
        ));
    }
    source.push_str("    ]\n    model.add_species(species)\n");
    let side = |terms: &BTreeMap<String, i64>| {
        let text: Vec<_> = terms
            .iter()
            .map(|(name, c)| format!("species[{}]: {c}", indices[name.as_str()]))
            .collect();
        format!("{{{}}}", text.join(", "))
    };
    for (i, reaction) in model.reactions.iter().enumerate() {
        source.push_str(&format!("    rate = gp.Parameter(name='K{i}', expression={:?})\n    model.add_parameter(rate)\n    model.add_reaction(gp.Reaction(name='R{i}', reactants={}, products={}, rate=rate))\n",
            coefficients[i], side(&reaction.reactants), side(&reaction.products)));
    }
    source.push_str("    return model\n");
    report.artifact = Some(PythonExport { target: "gillespy2".into(), target_version: GILLESPY2_VERSION.into(), source,
        species_aliases: aliases, reaction_aliases: model.reactions.iter().enumerate().map(|(i, r)| (r.name.clone(), format!("R{i}"))).collect(),
        notes: vec![
            "Target volume is fixed at 1: source coefficients already use stochastic counts; Omega does not rescale them.".into(),
            "Homodimer rate constants are divided by 2 to preserve the binomial factorial. Parameters are resolved into one constant per channel.".into(),
            "Use a stochastic solver. Target ODE semantics, live budgets, i64 range/underflow diagnostics and cross-runtime bit identity are not preserved.".into(),
            "The target uses floating-point populations; even supported initial counts may exceed its exact range during a reference run.".into(),
        ] });
    report.success = true;
    report
}
