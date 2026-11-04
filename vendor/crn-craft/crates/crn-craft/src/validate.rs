use crate::contracts::*;
use crate::diagnostic::{Diagnostic, pointer};
use std::collections::{BTreeMap, BTreeSet};

pub const MAX_SPECIES: usize = 4096;
pub const MAX_PARAMETERS: usize = 16384;
pub const MAX_REACTIONS: usize = 65536;
pub const MAX_SAVE_TIMES: usize = 1_000_000;

fn identifier(name: &str) -> bool {
    let mut bytes = name.bytes();
    name.len() <= 64
        && bytes
            .next()
            .is_some_and(|x| x.is_ascii_alphabetic() || x == b'_')
        && bytes.all(|x| x.is_ascii_alphanumeric() || x == b'_')
}

fn name_check(name: &str, path: &str, out: &mut Vec<Diagnostic>) {
    if !identifier(name) {
        out.push(Diagnostic::error(
            "CRN-E010",
            path,
            "invalid identifier",
            "Use 1–64 ASCII letters, digits or underscores, starting with a letter or underscore.",
        ));
    }
}

fn species_ref(name: &str, path: &str, names: &BTreeSet<&str>, out: &mut Vec<Diagnostic>) {
    if !names.contains(name) {
        out.push(Diagnostic::error(
            "CRN-E012",
            path,
            format!("unknown species {name:?}"),
            "Declare the species or correct its name.",
        ));
    }
}

fn scalar(
    value: &Scalar,
    path: &str,
    positive: bool,
    parameters: &BTreeMap<String, f64>,
    out: &mut Vec<Diagnostic>,
) -> Option<f64> {
    let number = match value {
        Scalar::Literal(x) => *x,
        Scalar::Parameter(name) => match parameters.get(name) {
            Some(x) => *x,
            None => {
                out.push(Diagnostic::error(
                    "CRN-E014",
                    path,
                    format!("unknown parameter {name:?}"),
                    "Declare the parameter or provide a numeric literal.",
                ));
                return None;
            }
        },
    };
    if !number.is_finite() || number < 0.0 || (positive && number == 0.0) {
        out.push(Diagnostic::error(
            "CRN-E015",
            path,
            "rate argument outside its finite numeric domain",
            if positive {
                "Use a finite, strictly positive value."
            } else {
                "Use a finite, non-negative value."
            },
        ));
        None
    } else {
        Some(number)
    }
}

/// Structural and domain checks only. Does not evaluate a propensity or prove boundedness.
pub fn validate_model(model: &Model) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    name_check(&model.name, "/name", &mut out);
    if !model.omega.is_finite() || model.omega <= 0.0 {
        out.push(Diagnostic::error(
            "CRN-E015",
            "/omega",
            "invalid system size",
            "Use a finite positive counts-per-concentration scale.",
        ));
    }
    if model.species.is_empty()
        || model.species.len() > MAX_SPECIES
        || model.parameters.len() > MAX_PARAMETERS
        || model.reactions.len() > MAX_REACTIONS
    {
        out.push(Diagnostic::error(
            "CRN-E004",
            "",
            "model collection limit exceeded or no species",
            "Use 1–4096 species, at most 16384 parameters and 65536 reactions.",
        ));
    }
    let mut names = BTreeSet::new();
    for (i, item) in model.species.iter().enumerate() {
        let base = format!("/species/{i}");
        name_check(&item.name, &pointer(&base, "name"), &mut out);
        if !names.insert(item.name.as_str()) {
            out.push(Diagnostic::error(
                "CRN-E011",
                &pointer(&base, "name"),
                "duplicate species name",
                "Give each species a unique name.",
            ));
        }
        if item.initial_count < 0 {
            out.push(Diagnostic::error(
                "CRN-E016",
                &pointer(&base, "initial_count"),
                "negative initial population",
                "Use a non-negative i64 count.",
            ));
        }
        if item
            .population_limit
            .is_some_and(|limit| limit < 0 || limit < item.initial_count)
        {
            out.push(Diagnostic::error(
                "CRN-E016",
                &pointer(&base, "population_limit"),
                "population limit below the initial population",
                "Use a non-negative limit at least as large as initial_count.",
            ));
        }
    }
    for (name, value) in &model.parameters {
        let path = pointer("/parameters", name);
        name_check(name, &path, &mut out);
        if names.contains(name.as_str()) {
            out.push(Diagnostic::error(
                "CRN-E011",
                &path,
                "parameter and species share a name",
                "Keep parameter and species namespaces disjoint.",
            ));
        }
        if !value.is_finite() {
            out.push(Diagnostic::error(
                "CRN-E015",
                &path,
                "non-finite parameter",
                "Use a finite number.",
            ));
        }
    }
    let mut reaction_names = BTreeSet::new();
    let mut increased = BTreeSet::new();
    let mut decreased = BTreeSet::new();
    for (i, reaction) in model.reactions.iter().enumerate() {
        let base = format!("/reactions/{i}");
        name_check(&reaction.name, &pointer(&base, "name"), &mut out);
        if !reaction_names.insert(&reaction.name) {
            out.push(Diagnostic::error(
                "CRN-E011",
                &pointer(&base, "name"),
                "duplicate reaction name",
                "Give each reaction a unique name.",
            ));
        }
        for (label, side) in [
            ("reactants", &reaction.reactants),
            ("products", &reaction.products),
        ] {
            for (name, value) in side {
                let path = pointer(&pointer(&base, label), name);
                species_ref(name, &path, &names, &mut out);
                if *value <= 0 || *value > i64::from(i32::MAX) {
                    out.push(Diagnostic::error(
                        "CRN-E013",
                        &path,
                        "invalid stoichiometry",
                        "Use positive integers up to 2147483647; omit zero entries.",
                    ));
                }
            }
        }
        if reaction.reactants == reaction.products {
            out.push(Diagnostic::warning(
                "CRN-W003",
                &base,
                "reaction has no population effect",
                "Remove it unless event counting has a stated purpose.",
            ));
        }
        for name in reaction.reactants.keys().chain(reaction.products.keys()) {
            let left = reaction.reactants.get(name).copied().unwrap_or(0);
            let right = reaction.products.get(name).copied().unwrap_or(0);
            if right > left {
                increased.insert(name.as_str());
            }
            if right < left {
                decreased.insert(name.as_str());
            }
        }
        let base = pointer(&base, "rate");
        match &reaction.rate {
            Rate::MassAction { coefficient } => {
                if scalar(
                    coefficient,
                    &pointer(&base, "coefficient"),
                    false,
                    &model.parameters,
                    &mut out,
                ) == Some(0.0)
                {
                    out.push(Diagnostic::warning(
                        "CRN-W002",
                        &pointer(&base, "coefficient"),
                        "reaction coefficient is zero",
                        "Confirm that this disabled reaction is intentional.",
                    ));
                }
            }
            Rate::HillActivation {
                regulator,
                vmax,
                basal,
                k,
                n,
            }
            | Rate::HillRepression {
                regulator,
                vmax,
                basal,
                k,
                n,
            } => {
                species_ref(regulator, &pointer(&base, "regulator"), &names, &mut out);
                for (key, value, positive) in [
                    ("vmax", vmax, false),
                    ("basal", basal, false),
                    ("k", k, true),
                    ("n", n, true),
                ] {
                    scalar(
                        value,
                        &pointer(&base, key),
                        positive,
                        &model.parameters,
                        &mut out,
                    );
                }
            }
            Rate::MichaelisMenten {
                substrate,
                vmax,
                km,
            } => {
                species_ref(substrate, &pointer(&base, "substrate"), &names, &mut out);
                scalar(
                    vmax,
                    &pointer(&base, "vmax"),
                    false,
                    &model.parameters,
                    &mut out,
                );
                scalar(km, &pointer(&base, "km"), true, &model.parameters, &mut out);
            }
        }
    }
    for (i, item) in model.species.iter().enumerate() {
        if increased.contains(item.name.as_str()) && !decreased.contains(item.name.as_str()) {
            out.push(Diagnostic::warning(
                "CRN-W001",
                &format!("/species/{i}"),
                "species has a net-producing reaction and no net-consuming reaction",
                "Review growth and budgets; this structural warning does not prove unboundedness.",
            ));
        }
    }
    out
}

pub fn valid_relative_path(path: &str) -> bool {
    !path.is_empty()
        && !path.contains(['\\', ':', '\0'])
        && path
            .split('/')
            .all(|part| !part.is_empty() && part != "." && part != "..")
}

pub fn validate_run(model: &Model, run: &Run, model_hash: &str) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    if !valid_relative_path(&run.model.path) {
        out.push(Diagnostic::error("CRN-E022", "/model/path", "non-portable or non-relative model path", "Use a slash-separated path below the run directory without '.', '..', empty components, drive prefixes or backslashes."));
    }
    if run.model.sha256.len() != 64
        || !run
            .model
            .sha256
            .bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
    {
        out.push(Diagnostic::error(
            "CRN-E021",
            "/model/sha256",
            "invalid SHA-256 encoding",
            "Use exactly 64 lowercase hexadecimal characters.",
        ));
    } else if run.model.sha256 != model_hash {
        out.push(Diagnostic::error(
            "CRN-E021",
            "/model/sha256",
            "model byte hash mismatch",
            "Recompute SHA-256 from the exact model bytes after reviewing the change.",
        ));
    }
    let times = &run.time.save_times;
    if !run.time.end.is_finite() || run.time.end <= 0.0 {
        out.push(Diagnostic::error(
            "CRN-E020",
            "/time/end",
            "invalid end time",
            "Use a finite positive time in seconds.",
        ));
    }
    if times.len() < 2
        || times.len() > MAX_SAVE_TIMES
        || times.first() != Some(&0.0)
        || times.last() != Some(&run.time.end)
    {
        out.push(Diagnostic::error(
            "CRN-E020",
            "/time/save_times",
            "invalid save grid endpoints or size",
            "Supply 2–1000000 increasing times, starting at 0 and ending exactly at end.",
        ));
    }
    for (i, time) in times.iter().enumerate() {
        if !time.is_finite()
            || *time < 0.0
            || *time > run.time.end
            || (i > 0 && *time <= times[i - 1])
        {
            out.push(Diagnostic::error(
                "CRN-E020",
                &format!("/time/save_times/{i}"),
                "invalid or non-increasing save time",
                "Use finite, strictly increasing values within [0, end].",
            ));
        }
    }
    if run.budget.threads == 0 || run.budget.threads > 1024 {
        out.push(Diagnostic::error(
            "CRN-E020",
            "/budget/threads",
            "invalid thread allowance",
            "Use 1–1024 threads.",
        ));
    }
    if run.budget.max_output_bytes == 0 {
        out.push(Diagnostic::error(
            "CRN-E020",
            "/budget/max_output_bytes",
            "zero output budget",
            "Provide a positive byte budget.",
        ));
    }
    if run
        .budget
        .wall_time_s
        .is_some_and(|time| !time.is_finite() || time <= 0.0)
    {
        out.push(Diagnostic::error(
            "CRN-E020",
            "/budget/wall_time_s",
            "invalid wall-time budget",
            "Omit the limit or use a finite positive number of seconds.",
        ));
    }
    let mut seen = BTreeSet::new();
    for (i, statistic) in run.outputs.statistics.iter().enumerate() {
        if !seen.insert(statistic) {
            out.push(Diagnostic::error(
                "CRN-E020",
                &format!("/outputs/statistics/{i}"),
                "duplicate statistic",
                "Request each statistic once.",
            ));
        }
    }
    if run.outputs.statistics.is_empty() && run.outputs.keep_trajectories == 0 {
        out.push(Diagnostic::error(
            "CRN-E020",
            "/outputs",
            "no output requested",
            "Request statistics or retained trajectories.",
        ));
    }
    // A conservative minimum payload, not a claim about JSON/archive size.
    let values_per_sample = (model.species.len() as u128)
        * (run.outputs.statistics.len() as u128 + u128::from(run.outputs.keep_trajectories));
    let min_bytes = (times.len() as u128) * (values_per_sample + 1) * 8;
    if min_bytes > u128::from(run.budget.max_output_bytes) {
        out.push(Diagnostic::error(
            "CRN-E020",
            "/budget/max_output_bytes",
            format!("output budget below minimum numeric payload ({min_bytes} bytes)"),
            "Reduce retained output or increase the budget; serialization needs additional space.",
        ));
    }
    match &run.method {
        Method::Ssa { .. } => {
            match &run.ensemble {
                Some(ensemble) if ensemble.replicates > 0 => {
                    if run.outputs.keep_trajectories > ensemble.replicates {
                        out.push(Diagnostic::error(
                            "CRN-E020",
                            "/outputs/keep_trajectories",
                            "more kept trajectories than replicates",
                            "Keep at most the requested replicate count.",
                        ));
                    }
                    if ensemble.replicates < 2 && seen.contains(&Statistic::Variance) {
                        out.push(Diagnostic::error(
                            "CRN-E020",
                            "/ensemble/replicates",
                            "sample variance needs at least two replicates",
                            "Use two or more replicates or omit variance.",
                        ));
                    }
                }
                _ => out.push(Diagnostic::error(
                    "CRN-E020",
                    "/ensemble",
                    "SSA requires a positive replicate count",
                    "Supply ensemble.replicates >= 1 and root_seed.",
                )),
            }
            if run.budget.max_events_per_replicate.is_none_or(|n| n == 0) {
                out.push(Diagnostic::error(
                    "CRN-E020",
                    "/budget/max_events_per_replicate",
                    "SSA requires a positive event cap",
                    "Supply a positive event cap per replicate.",
                ));
            }
        }
        Method::Fsp {
            bounds,
            max_states,
            tolerance,
        } => {
            if run.ensemble.is_some()
                || run.budget.max_events_per_replicate.is_some()
                || run.outputs.keep_trajectories != 0
            {
                out.push(Diagnostic::error(
                    "CRN-E020",
                    "/method",
                    "FSP does not use ensembles, event caps or trajectories",
                    "Omit ensemble and max_events_per_replicate; set keep_trajectories to 0.",
                ));
            }
            if *max_states == 0 || !tolerance.is_finite() || *tolerance <= 0.0 || *tolerance >= 1.0
            {
                out.push(Diagnostic::error(
                    "CRN-E023",
                    "/method",
                    "invalid FSP capacity or tolerance",
                    "Use max_states >= 1 and 0 < tolerance < 1.",
                ));
            }
            let names = model.species.iter().map(|x| x.name.as_str()).collect();
            for name in bounds.keys() {
                species_ref(name, &pointer("/method/bounds", name), &names, &mut out);
            }
            for item in &model.species {
                if bounds
                    .get(&item.name)
                    .is_none_or(|bound| *bound < 0 || *bound < item.initial_count)
                {
                    out.push(Diagnostic::error("CRN-E023", &pointer("/method/bounds", &item.name), "missing FSP bound or initial state outside truncation", "Give every species a finite integer truncation bound >= initial_count; escaping probability goes to a sink."));
                }
            }
        }
    }
    out
}
