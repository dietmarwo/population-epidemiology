//! Six fixed files, checksummed manifest, completion marker written last.
//! Durable-write pattern inspected in optimization-jobs storage.rs at cc6b03d.
use crate::{
    diagnostic::Diagnostic,
    parse_model, parse_run, sampling, sha256,
    simulation::{self, Assessment, Control, Provenance, SimulationResult, Stop},
    strict_json,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeSet,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::Path,
};

const PAYLOAD_NAMES: [&str; 4] = ["model.json", "run.json", "result.json", "provenance.json"];
const ALL_NAMES: [&str; 6] = [
    "model.json",
    "run.json",
    "result.json",
    "provenance.json",
    "manifest.json",
    "COMPLETE",
];
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum ArchiveVersion {
    #[serde(rename = "crn-craft.archive/0.1")]
    V01,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ArchiveFile {
    pub path: String,
    pub bytes: u64,
    pub sha256: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub schema: ArchiveVersion,
    pub files: Vec<ArchiveFile>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ArchiveReport {
    pub schema: String,
    pub valid: bool,
    pub archive_complete: bool,
    pub simulation_complete: bool,
    pub termination: Option<Stop>,
    pub assessment: Option<Assessment>,
    pub model_sha256: Option<String>,
    pub run_sha256: Option<String>,
    pub diagnostics: Vec<Diagnostic>,
}
impl ArchiveReport {
    pub fn failure(diagnostics: Vec<Diagnostic>) -> Self {
        Self {
            schema: "crn-craft.archive-check/0.1".into(),
            valid: false,
            archive_complete: false,
            simulation_complete: false,
            termination: None,
            assessment: None,
            model_sha256: None,
            run_sha256: None,
            diagnostics,
        }
    }
    fn success(result: &SimulationResult) -> Self {
        Self {
            schema: "crn-craft.archive-check/0.1".into(),
            valid: true,
            archive_complete: true,
            simulation_complete: result.termination == Stop::Completed,
            termination: Some(result.termination),
            assessment: Some(result.assessment),
            model_sha256: Some(result.model_sha256.clone()),
            run_sha256: Some(result.run_sha256.clone()),
            diagnostics: Vec::new(),
        }
    }
}

fn error(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E065",
        "/archive",
        message,
        "Use a new output directory; preserve all files and verify the completion marker and checksums.",
    )
}
fn io_error(e: std::io::Error) -> Diagnostic {
    error(e.to_string())
}

struct LimitedBuffer {
    bytes: Vec<u8>,
    limit: u64,
}
impl Write for LimitedBuffer {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if self.bytes.len() as u64 + bytes.len() as u64 > self.limit {
            return Err(std::io::Error::other(
                "serialized output byte budget exhausted",
            ));
        }
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encoded<T: Serialize>(value: &T, limit: u64) -> Result<Vec<u8>, Diagnostic> {
    let mut buffer = LimitedBuffer {
        bytes: Vec::new(),
        limit,
    };
    serde_json::to_writer(&mut buffer, value).map_err(|e| error(e.to_string()))?;
    buffer.write_all(b"\n").map_err(io_error)?;
    Ok(buffer.bytes)
}
fn write_new(path: &Path, bytes: &[u8]) -> Result<(), Diagnostic> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(io_error)?;
    file.write_all(bytes).map_err(io_error)?;
    file.sync_all().map_err(io_error)
}
fn read_regular(path: &Path, limit: u64) -> Result<Vec<u8>, Diagnostic> {
    let metadata = fs::symlink_metadata(path).map_err(io_error)?;
    if !metadata.is_file() || metadata.len() > limit {
        return Err(error("not a bounded regular archive file"));
    }
    let file = File::open(path).map_err(io_error)?;
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(io_error)?;
    if bytes.len() as u64 > limit {
        return Err(error("archive file grew beyond its bound"));
    }
    Ok(bytes)
}
fn executable_hash() -> Option<String> {
    let mut file = File::open(std::env::current_exe().ok()?).ok()?;
    let mut digest = Sha256::new();
    let mut buffer = [0; 65536];
    loop {
        let n = file.read(&mut buffer).ok()?;
        if n == 0 {
            break;
        }
        digest.update(&buffer[..n]);
    }
    Some(format!("{:x}", digest.finalize()))
}

/// Preserves the exact consumed bytes; run path identity is checked before use.
pub fn simulate_files(
    model_path: &Path,
    run_path: &Path,
    output: &Path,
    control: &Control,
) -> ArchiveReport {
    simulate_using(model_path, run_path, output, control, false)
}
/// Optional dependency-cached execution of the same Direct stochastic contract.
pub fn simulate_cached_files(
    model_path: &Path,
    run_path: &Path,
    output: &Path,
    control: &Control,
) -> ArchiveReport {
    simulate_using(model_path, run_path, output, control, true)
}
fn simulate_using(
    model_path: &Path,
    run_path: &Path,
    output: &Path,
    control: &Control,
    cached: bool,
) -> ArchiveReport {
    let action = (|| -> Result<ArchiveReport, Vec<Diagnostic>> {
        let model = crate::read_document(model_path).map_err(|e| vec![e])?;
        let run = crate::read_document(run_path).map_err(|e| vec![e])?;
        let prepared = simulation::prepare(&model, &run)?;
        let parent = run_path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let identity = (
            parent.canonicalize(),
            parent.join(&prepared.run.model.path).canonicalize(),
            model_path.canonicalize(),
        );
        if !matches!(identity, (Ok(ref base), Ok(ref referenced), Ok(ref supplied)) if referenced.starts_with(base) && referenced == supplied)
        {
            return Err(vec![Diagnostic::error(
                "CRN-E022",
                "/model/path",
                "run model reference does not identify the supplied model inside the run directory",
                "Use the referenced model below the run directory.",
            )]);
        }
        fs::create_dir(output).map_err(|e| vec![io_error(e)])?;
        write_new(&output.join("model.json"), &model).map_err(|e| vec![e])?;
        write_new(&output.join("run.json"), &run).map_err(|e| vec![e])?;
        let mut simulation = if cached {
            simulation::execute_cached(&prepared, control)
        } else {
            simulation::execute(&prepared, control)
        }
        .map_err(|e| vec![e])?;
        simulation.provenance.executable_sha256 = executable_hash();
        let provenance = encoded(&simulation.provenance, 8192).map_err(|e| vec![e])?;
        let total_limit = prepared
            .run
            .budget
            .max_output_bytes
            .min(simulation::MAX_RESULT_BYTES);
        let result_limit = total_limit
            .saturating_sub(model.len() as u64 + run.len() as u64 + provenance.len() as u64 + 8192);
        let result = match encoded(&simulation.result, result_limit) {
            Ok(bytes) => bytes,
            Err(_) => {
                // Retain full accounting; never present truncated statistics as exact.
                simulation.result.termination = Stop::OutputLimit;
                simulation.result.assessment = Assessment::Incomplete;
                simulation.result.statistics = None;
                simulation.result.trajectories.clear();
                simulation.result.diagnostics.push(Diagnostic::error("CRN-E064", "/outputs", "serialized numerical payload exceeded the output byte budget", "Increase max_output_bytes; trajectories/statistics were removed, replicate accounting was retained."));
                encoded(&simulation.result, result_limit).map_err(|e| vec![e])?
            }
        };
        validate_result(&model, &run, &simulation.result).map_err(|e| vec![e])?;
        write_new(&output.join("result.json"), &result).map_err(|e| vec![e])?;
        write_new(&output.join("provenance.json"), &provenance).map_err(|e| vec![e])?;
        let files = PAYLOAD_NAMES
            .iter()
            .zip([&model, &run, &result, &provenance])
            .map(|(name, bytes)| ArchiveFile {
                path: (*name).into(),
                bytes: bytes.len() as u64,
                sha256: sha256(bytes),
            })
            .collect();
        let manifest = encoded(
            &Manifest {
                schema: ArchiveVersion::V01,
                files,
            },
            8192,
        )
        .map_err(|e| vec![e])?;
        let marker = format!("{}\n", sha256(&manifest));
        if model.len() as u64
            + run.len() as u64
            + result.len() as u64
            + provenance.len() as u64
            + manifest.len() as u64
            + marker.len() as u64
            > total_limit
        {
            return Err(vec![error("archive byte budget exceeded")]);
        }
        write_new(&output.join("manifest.json"), &manifest).map_err(|e| vec![e])?;
        write_new(&output.join("COMPLETE"), marker.as_bytes()).map_err(|e| vec![e])?;
        #[cfg(unix)]
        File::open(output)
            .and_then(|f| f.sync_all())
            .map_err(|e| vec![io_error(e)])?;
        Ok(ArchiveReport::success(&simulation.result))
    })();
    action.unwrap_or_else(ArchiveReport::failure)
}

/// Integrity and relational checks, not cryptographic authenticity or proof of
/// stochastic correctness. Never executes supplied code or follows manifest paths.
pub fn verify_archive(path: &Path) -> ArchiveReport {
    let action = (|| -> Result<ArchiveReport, Diagnostic> {
        if !fs::symlink_metadata(path).map_err(io_error)?.is_dir() {
            return Err(error("archive root is not a real directory"));
        }
        let mut names = BTreeSet::new();
        let mut total_bytes = 0_u64;
        for entry in fs::read_dir(path).map_err(io_error)? {
            let entry = entry.map_err(io_error)?;
            let name = entry
                .file_name()
                .into_string()
                .map_err(|_| error("non-UTF8 archive name"))?;
            if !ALL_NAMES.contains(&name.as_str())
                || !entry.file_type().map_err(io_error)?.is_file()
            {
                return Err(error("unexpected file, directory or symlink in archive"));
            }
            total_bytes = total_bytes
                .checked_add(entry.metadata().map_err(io_error)?.len())
                .ok_or_else(|| error("archive length overflow"))?;
            names.insert(name);
        }
        if names.len() != ALL_NAMES.len() || total_bytes > simulation::MAX_RESULT_BYTES {
            return Err(error("incomplete or oversized archive"));
        }
        let manifest_bytes = read_regular(&path.join("manifest.json"), 8192)?;
        let marker = read_regular(&path.join("COMPLETE"), 65)?;
        if marker != format!("{}\n", sha256(&manifest_bytes)).as_bytes() {
            return Err(error("completion marker does not match manifest"));
        }
        let manifest: Manifest = strict_json::parse(&manifest_bytes)?;
        let paths: BTreeSet<_> = manifest.files.iter().map(|f| f.path.as_str()).collect();
        if manifest.files.len() != 4 || paths != PAYLOAD_NAMES.into_iter().collect() {
            return Err(error(
                "manifest paths are not exactly the four unique payload names",
            ));
        }
        let mut payloads = std::collections::BTreeMap::new();
        for entry in manifest.files {
            let limit = match entry.path.as_str() {
                "result.json" => simulation::MAX_RESULT_BYTES,
                "provenance.json" => 8192,
                _ => strict_json::MAX_DOCUMENT_BYTES as u64,
            };
            let bytes = read_regular(&path.join(&entry.path), limit)?;
            if entry.bytes != bytes.len() as u64 || entry.sha256 != sha256(&bytes) {
                return Err(error("payload checksum or length mismatch"));
            }
            payloads.insert(entry.path, bytes);
        }
        let result: SimulationResult = strict_json::parse_bounded(
            &payloads["result.json"],
            simulation::MAX_RESULT_BYTES as usize,
        )?;
        let provenance: Provenance = strict_json::parse(&payloads["provenance.json"])?;
        let run = parse_run(&payloads["run.json"])?;
        if total_bytes > run.budget.max_output_bytes
            || !provenance.elapsed_seconds.is_finite()
            || provenance.elapsed_seconds < 0.0
            || provenance.workers == 0
            || provenance.workers > (run.budget.threads as usize).min(simulation::MAX_WORKERS)
            || provenance.workspace_cap_bytes != simulation::MAX_WORKSPACE_BYTES
            || provenance.reserved_workspace_bytes > provenance.workspace_cap_bytes
        {
            return Err(error("provenance or archive budget is inconsistent"));
        }
        validate_result(&payloads["model.json"], &payloads["run.json"], &result)?;
        Ok(ArchiveReport::success(&result))
    })();
    action.unwrap_or_else(|e| ArchiveReport::failure(vec![e]))
}

pub fn validate_result(
    model_bytes: &[u8],
    run_bytes: &[u8],
    result: &SimulationResult,
) -> Result<(), Diagnostic> {
    let prepared = simulation::prepare(model_bytes, run_bytes)
        .map_err(|_| error("invalid archived inputs"))?;
    let run = &prepared.run;
    let model = parse_model(model_bytes)?;
    let ensemble = run.ensemble.as_ref().unwrap();
    let n = model.species.len();
    let t = run.time.save_times.len();
    if result.model_sha256 != sha256(model_bytes)
        || result.run_sha256 != sha256(run_bytes)
        || result.sampler != sampling::SAMPLER
        || result.root_seed != ensemble.root_seed
        || result.species != prepared.network.species()
        || result.save_times != run.time.save_times
        || result.requested_replicates != ensemble.replicates
        || result.replicates.len() as u64 != ensemble.replicates
    {
        return Err(error("result identity or dimensions differ from inputs"));
    }
    let mut completed = 0;
    let mut events = 0_u64;
    let mut aggregate = Stop::Completed;
    for (id, record) in result.replicates.iter().enumerate() {
        if record.id != id as u64
            || !record.stop_time.is_finite()
            || record.stop_time < 0.0
            || record.stop_time > run.time.end
            || record.events > run.budget.max_events_per_replicate.unwrap()
            || record.saved_points > t
            || (!record.started
                && (record.events != 0
                    || record.saved_points != 0
                    || !matches!(record.termination, Stop::Cancelled | Stop::WallTimeLimit)))
            || record.termination == Stop::OutputLimit
        {
            return Err(error("invalid replicate accounting"));
        }
        if record.saved_points > 0
            && run.time.save_times[record.saved_points - 1] > record.stop_time
        {
            return Err(error("saved state lies beyond observation time"));
        }
        if record.termination == Stop::Completed {
            if !record.started
                || record.stop_time != run.time.end
                || record.saved_points != t
                || record.diagnostic.is_some()
            {
                return Err(error("incomplete replicate marked complete"));
            }
            completed += 1;
        }
        if matches!(record.termination, Stop::Failed | Stop::PopulationLimit)
            && record.diagnostic.is_none()
        {
            return Err(error("failed/stopped replicate has no diagnostic"));
        }
        if record.termination == Stop::EventLimit
            && record.events != run.budget.max_events_per_replicate.unwrap()
        {
            return Err(error("event-limit accounting mismatch"));
        }
        events = events
            .checked_add(record.events)
            .ok_or_else(|| error("event accounting overflow"))?;
        aggregate = aggregate.max(record.termination);
    }
    if result.completed_replicates != completed || result.events != events {
        return Err(error("aggregate accounting mismatch"));
    }
    let output_limited = result.termination == Stop::OutputLimit;
    let expected_assessment = if output_limited {
        Assessment::Incomplete
    } else if aggregate == Stop::Completed {
        Assessment::ExactEnsemble
    } else if aggregate == Stop::Failed {
        Assessment::FailedChecks
    } else {
        Assessment::Incomplete
    };
    if (!output_limited && result.termination != aggregate)
        || result.assessment != expected_assessment
        || (output_limited
            && (result.statistics.is_some()
                || !result.trajectories.is_empty()
                || !result.diagnostics.iter().any(|d| d.code == "CRN-E064")))
    {
        return Err(error("termination/assessment evidence mismatch"));
    }
    if result.trajectories.len() as u64
        != if output_limited {
            0
        } else {
            run.outputs.keep_trajectories
        }
    {
        return Err(error("retained trajectory count mismatch"));
    }
    for (id, trajectory) in result.trajectories.iter().enumerate() {
        if trajectory.id != id as u64
            || trajectory.counts.len() != result.replicates[id].saved_points
        {
            return Err(error("trajectory identity/length mismatch"));
        }
        for counts in &trajectory.counts {
            prepared
                .network
                .validate_state(counts)
                .map_err(|_| error("invalid retained count state"))?;
        }
        if trajectory
            .counts
            .first()
            .is_some_and(|counts| counts != prepared.network.initial_counts())
        {
            return Err(error("trajectory starts at the wrong initial state"));
        }
    }
    let want_statistics =
        aggregate == Stop::Completed && !output_limited && !run.outputs.statistics.is_empty();
    if result.statistics.is_some() != want_statistics {
        return Err(error(
            "statistics present for incomplete work or absent for completed work",
        ));
    }
    if let Some(stats) = &result.statistics {
        let mean = run
            .outputs
            .statistics
            .contains(&crate::contracts::Statistic::Mean);
        let variance = run
            .outputs
            .statistics
            .contains(&crate::contracts::Statistic::Variance);
        if stats.samples != ensemble.replicates
            || stats.mean.is_some() != mean
            || stats.sample_variance.is_some() != variance
            || stats.mean_standard_error.is_some() != (mean && stats.samples > 1)
        {
            return Err(error("statistic estimator/sample count mismatch"));
        }
        for values in [
            &stats.mean,
            &stats.sample_variance,
            &stats.mean_standard_error,
        ]
        .into_iter()
        .flatten()
        {
            if values.len() != t
                || values
                    .iter()
                    .any(|row| row.len() != n || row.iter().any(|x| !x.is_finite() || *x < 0.0))
            {
                return Err(error("invalid statistic dimensions or values"));
            }
        }
        if result.trajectories.len() as u64 == ensemble.replicates {
            let mut means = vec![0.0; n * t];
            let mut m2 = vec![0.0; n * t];
            for (id, path) in result.trajectories.iter().enumerate() {
                for (i, count) in path.counts.iter().flatten().enumerate() {
                    let delta = *count as f64 - means[i];
                    means[i] += delta / (id + 1) as f64;
                    m2[i] += delta * (*count as f64 - means[i]);
                }
            }
            let variance: Vec<_> = m2
                .iter()
                .map(|x| x / (ensemble.replicates - 1).max(1) as f64)
                .collect();
            if stats
                .mean
                .as_ref()
                .is_some_and(|rows| !rows.iter().flatten().eq(means.iter()))
                || stats
                    .sample_variance
                    .as_ref()
                    .is_some_and(|rows| !rows.iter().flatten().eq(variance.iter()))
                || stats.mean_standard_error.as_ref().is_some_and(|rows| {
                    !rows.iter().flatten().copied().eq(variance
                        .iter()
                        .map(|v| libm::sqrt(v / ensemble.replicates as f64)))
                })
            {
                return Err(error(
                    "statistics disagree with the complete retained ensemble",
                ));
            }
        }
    }
    Ok(())
}
