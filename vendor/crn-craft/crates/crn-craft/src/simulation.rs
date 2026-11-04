//! Direct SSA. Fixed replicate streams and ordered reduction, bounded queues.
use crate::{
    CompiledNetwork, check_bytes,
    contracts::{Method, Run, Statistic},
    diagnostic::Diagnostic,
    parse_model, parse_run, sampling, sha256,
};
use rand_core::RngCore;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{
    path::PathBuf,
    sync::{
        atomic::{AtomicBool, Ordering},
        mpsc,
    },
    time::{Duration, Instant},
};

pub const MAX_REPLICATES: u64 = 1_000_000;
pub const MAX_WORKSPACE_BYTES: u64 = 256 * 1024 * 1024;
pub const MAX_RESULT_BYTES: u64 = 256 * 1024 * 1024;
pub const MAX_WORKERS: usize = 64;
/// Polls between reads of the wall clock and, at most once per
/// `CANCEL_FILE_INTERVAL`, the cancel file. The flag is read at every poll.
pub const CLOCK_POLL_STRIDE: u32 = 64;
pub const CANCEL_FILE_INTERVAL: Duration = Duration::from_millis(1);

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ResultVersion {
    #[serde(rename = "crn-craft.result/0.2")]
    V02,
}
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Stop {
    Completed,
    EventLimit,
    PopulationLimit,
    WallTimeLimit,
    Cancelled,
    Failed,
    OutputLimit,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Assessment {
    ExactEnsemble,
    Incomplete,
    FailedChecks,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Replicate {
    pub id: u64,
    pub started: bool,
    pub termination: Stop,
    pub events: u64,
    pub stop_time: f64,
    pub saved_points: usize,
    pub diagnostic: Option<Diagnostic>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Trajectory {
    pub id: u64,
    pub counts: Vec<Vec<i64>>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Statistics {
    pub samples: u64,
    pub mean: Option<Vec<Vec<f64>>>,
    pub sample_variance: Option<Vec<Vec<f64>>>,
    pub mean_standard_error: Option<Vec<Vec<f64>>>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SimulationResult {
    pub schema: ResultVersion,
    pub model_sha256: String,
    pub run_sha256: String,
    pub sampler: String,
    pub root_seed: u64,
    pub species: Vec<String>,
    pub save_times: Vec<f64>,
    pub termination: Stop,
    pub assessment: Assessment,
    pub requested_replicates: u64,
    pub completed_replicates: u64,
    pub events: u64,
    pub replicates: Vec<Replicate>,
    pub trajectories: Vec<Trajectory>,
    pub statistics: Option<Statistics>,
    pub diagnostics: Vec<Diagnostic>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    pub crate_version: String,
    pub rustc: String,
    pub target_os: String,
    pub target_arch: String,
    pub cargo_lock_sha256: String,
    pub executable_sha256: Option<String>,
    pub workers: usize,
    pub elapsed_seconds: f64,
    pub reserved_workspace_bytes: u64,
    pub workspace_cap_bytes: u64,
}
pub struct Simulation {
    pub result: SimulationResult,
    pub provenance: Provenance,
}

/// Cooperative API cancellation; CLI additionally accepts a caller-owned file.
#[derive(Default)]
pub struct Control {
    cancelled: AtomicBool,
    pub cancel_file: Option<PathBuf>,
    /// Optional integration checkpoint. Invoked at the amortized clock stride.
    pub external_stop: Option<Box<dyn Fn() -> bool + Send + Sync>>,
    #[cfg(test)]
    file_checks: std::sync::atomic::AtomicU64,
    #[cfg(test)]
    polls: std::sync::atomic::AtomicU64,
}
impl Control {
    /// Requests cancellation. Observing the cancel file sets the same flag.
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }
}

/// Per-worker polling state. A stat of the cancel file costs about a
/// microsecond and the clock tens of nanoseconds, against roughly a hundred
/// nanoseconds per event for small models, so both are amortized. The first
/// poll of each worker reads both, and a stop, once observed, is kept.
pub(crate) struct Poller<'a> {
    control: &'a Control,
    start: Instant,
    wall: Option<f64>,
    countdown: u32,
    last_file_check: Option<Instant>,
    stopped: Option<Stop>,
}
impl<'a> Poller<'a> {
    pub(crate) fn new(control: &'a Control, start: Instant, wall: Option<f64>) -> Self {
        Self {
            control,
            start,
            wall,
            countdown: 0,
            last_file_check: None,
            stopped: None,
        }
    }
    pub(crate) fn poll(&mut self) -> Option<Stop> {
        #[cfg(test)]
        self.control.polls.fetch_add(1, Ordering::Relaxed);
        if self.stopped.is_some() {
            return self.stopped;
        }
        if self.control.cancelled.load(Ordering::Relaxed) {
            self.stopped = Some(Stop::Cancelled);
            return self.stopped;
        }
        if self.countdown > 0
            || (self.wall.is_none()
                && self.control.cancel_file.is_none()
                && self.control.external_stop.is_none())
        {
            self.countdown = self.countdown.saturating_sub(1);
            return None;
        }
        self.countdown = CLOCK_POLL_STRIDE - 1;
        if self.control.external_stop.as_ref().is_some_and(|f| f()) {
            self.control.cancel();
            self.stopped = Some(Stop::Cancelled);
            return self.stopped;
        }
        let now = Instant::now();
        if let Some(path) = &self.control.cancel_file {
            if self
                .last_file_check
                .is_none_or(|last| now.duration_since(last) >= CANCEL_FILE_INTERVAL)
            {
                self.last_file_check = Some(now);
                #[cfg(test)]
                self.control.file_checks.fetch_add(1, Ordering::Relaxed);
                if path.exists() {
                    self.control.cancel();
                    self.stopped = Some(Stop::Cancelled);
                    return self.stopped;
                }
            }
        }
        if self
            .wall
            .is_some_and(|seconds| now.duration_since(self.start).as_secs_f64() >= seconds)
        {
            self.stopped = Some(Stop::WallTimeLimit);
        }
        self.stopped
    }
}

pub(crate) fn execution_error(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E060",
        "/simulation",
        message,
        "Use Direct SSA and reduce replicas, saved output or worker count to fit the published execution limits.",
    )
}

#[derive(Debug)]
pub struct Prepared {
    pub(crate) network: CompiledNetwork,
    pub(crate) run: Run,
    pub(crate) model_hash: String,
    pub(crate) run_hash: String,
    pub(crate) workers: usize,
    pub(crate) workspace: u64,
    pub(crate) queue_capacity: usize,
    pub(crate) input_bytes: u64,
}

pub fn prepare(model_bytes: &[u8], run_bytes: &[u8]) -> Result<Prepared, Vec<Diagnostic>> {
    let report = check_bytes(model_bytes, Some(run_bytes));
    if !report.valid {
        return Err(report.diagnostics);
    }
    let model = parse_model(model_bytes).map_err(|e| vec![e])?;
    let run = parse_run(run_bytes).map_err(|e| vec![e])?;
    // Cheap admission first, as before parameter templates: an oversized
    // request is rejected without compiling its dependency graph.
    admission(
        model.species.len(),
        model.reactions.len(),
        &run,
        model_bytes.len(),
        run_bytes.len(),
    )?;
    let network = CompiledNetwork::compile(&model)?;
    prepare_compiled(
        network,
        run,
        sha256(model_bytes),
        sha256(run_bytes),
        model_bytes.len(),
        run_bytes.len(),
    )
}

struct Admission {
    workers: usize,
    workspace: u128,
    queue_capacity: usize,
}

/// Workspace and evidence reservations; depends only on sizes and the run.
fn admission(
    species: usize,
    reactions: usize,
    run: &Run,
    model_len: usize,
    run_len: usize,
) -> Result<Admission, Vec<Diagnostic>> {
    if !matches!(run.method, Method::Ssa { .. }) {
        return Err(vec![execution_error("FSP execution is not implemented")]);
    }
    let count = run
        .ensemble
        .as_ref()
        .expect("validated ensemble")
        .replicates;
    let workers = (run.budget.threads as usize)
        .min(MAX_WORKERS)
        .min(count as usize);
    let cells = (species as u128) * run.time.save_times.len() as u128;
    // Conservative allocation reservation, not a claim to cap allocator/OS RSS.
    let reserve = |buffers_per_worker: u128| {
        (cells * 8 + run.time.save_times.len() as u128 * 24)
            * (buffers_per_worker * workers as u128 + u128::from(run.outputs.keep_trajectories))
            + cells * 24
            + u128::from(count) * 1024
            + workers as u128 * (2 * 1024 * 1024 + reactions as u128 * 8)
            + 32 * 1024 * 1024
    };
    // Preserve legacy admission and archive verification near the memory cap.
    let enlarged = reserve(6);
    let queue_capacity = if workers > 1 && enlarged <= u128::from(MAX_WORKSPACE_BYTES) {
        4
    } else {
        1
    };
    let workspace = if queue_capacity == 4 {
        enlarged
    } else {
        reserve(3)
    };
    let evidence_reserve = model_len as u128
        + run_len as u128
        + 16_384
        + u128::from(count) * 1024
        + species as u128 * 80
        + run.time.save_times.len() as u128 * 32;
    if count > MAX_REPLICATES
        || workspace > u128::from(MAX_WORKSPACE_BYTES)
        || evidence_reserve > u128::from(run.budget.max_output_bytes.min(MAX_RESULT_BYTES))
    {
        return Err(vec![execution_error(format!(
            "execution reservation exceeds limits: {workspace} workspace bytes, {evidence_reserve} evidence bytes"
        ))]);
    }
    // Ensure event accounting cannot wrap, even if every replicate reaches its cap.
    if count
        .checked_mul(run.budget.max_events_per_replicate.unwrap())
        .is_none()
    {
        return Err(vec![execution_error("aggregate event cap exceeds u64")]);
    }
    Ok(Admission {
        workers,
        workspace,
        queue_capacity,
    })
}

/// Admission for already validated structure/run inputs; parameter templates
/// invoke the same reservation checks as ordinary preparation on every bind.
pub(crate) fn prepare_compiled(
    network: CompiledNetwork,
    run: Run,
    model_hash: String,
    run_hash: String,
    model_len: usize,
    run_len: usize,
) -> Result<Prepared, Vec<Diagnostic>> {
    let Admission {
        workers,
        workspace,
        queue_capacity,
    } = admission(
        network.species().len(),
        network.reactions().len(),
        &run,
        model_len,
        run_len,
    )?;
    Ok(Prepared {
        network,
        run,
        model_hash,
        run_hash,
        workers,
        workspace: workspace as u64,
        queue_capacity,
        input_bytes: (model_len + run_len) as u64,
    })
}

pub(crate) struct PathResult {
    pub(crate) record: Replicate,
    /// Saved states, one row of `species` counts per saved point.
    pub(crate) counts: Vec<i64>,
}

/// Per-worker buffers reused across replicates. Each path resets what it reads.
#[derive(Default)]
struct Workspace {
    state: Vec<i64>,
    rates: Vec<f64>,
    prefix: Vec<f64>,
    minima: Vec<f64>,
}

fn trajectory(
    prepared: &Prepared,
    id: u64,
    poller: &mut Poller,
    cached: bool,
    work: &mut Workspace,
) -> PathResult {
    let run = &prepared.run;
    let mut rng = sampling::stream(run.ensemble.as_ref().unwrap().root_seed, id);
    path_with_rng_mode(&prepared.network, run, id, poller, &mut rng, cached, work)
}

#[cfg(test)]
fn path_with_rng(
    network: &CompiledNetwork,
    run: &Run,
    id: u64,
    poller: &mut Poller,
    rng: &mut impl RngCore,
) -> PathResult {
    path_with_rng_mode(
        network,
        run,
        id,
        poller,
        rng,
        false,
        &mut Workspace::default(),
    )
}

fn path_with_rng_mode(
    network: &CompiledNetwork,
    run: &Run,
    id: u64,
    poller: &mut Poller,
    rng: &mut impl RngCore,
    cached: bool,
    work: &mut Workspace,
) -> PathResult {
    let mut record = Replicate {
        id,
        started: false,
        termination: Stop::Completed,
        events: 0,
        stop_time: 0.0,
        saved_points: 0,
        diagnostic: None,
    };
    let mut state = std::mem::take(&mut work.state);
    state.clear();
    state.extend_from_slice(network.initial_counts());
    let mut counts = Vec::with_capacity(run.time.save_times.len() * state.len());
    let mut time = 0.0;
    let mut save = 0;
    // The first event recomputes every rate and running sum before use.
    let mut rates = std::mem::take(&mut work.rates);
    let mut prefix = std::mem::take(&mut work.prefix);
    let mut minima = std::mem::take(&mut work.minima);
    let mut total = 0.0;
    let mut previous: Option<usize> = None;
    let outcome = (|| -> Result<Stop, Diagnostic> {
        if let Some(stop) = poller.poll() {
            return Ok(stop);
        }
        record.started = true;
        loop {
            while save < run.time.save_times.len() && run.time.save_times[save] <= time {
                if let Some(stop) = poller.poll() {
                    return Ok(stop);
                }
                counts.extend_from_slice(&state);
                save += 1;
            }
            if time == run.time.end {
                return Ok(Stop::Completed);
            }
            if let Some(stop) = poller.poll() {
                return Ok(stop);
            }
            if record.events >= run.budget.max_events_per_replicate.unwrap() {
                return Ok(Stop::EventLimit);
            }
            if let (true, Some(changed)) = (cached, previous) {
                let affected = &network.dependencies()[changed];
                for (i, &reaction) in affected.iter().enumerate() {
                    // The event's poll above covers the first 64 channels.
                    if i > 0 && i % 64 == 0 {
                        if let Some(stop) = poller.poll() {
                            return Ok(stop);
                        }
                    }
                    rates[reaction] = network.propensity_in_valid_state(reaction, &state)?;
                }
                // Dependencies are sorted: running sums before the first one are unchanged.
                if let Some(&lowest) = affected.first() {
                    total = sampling::resume_total(&rates, &mut prefix, &mut minima, lowest)?;
                }
            } else {
                rates.clear();
                for reaction in 0..network.reactions().len() {
                    if reaction > 0 && reaction % 64 == 0 {
                        if let Some(stop) = poller.poll() {
                            return Ok(stop);
                        }
                    }
                    rates.push(network.propensity_in_valid_state(reaction, &state)?);
                }
                total = if cached {
                    prefix.resize(rates.len(), 0.0);
                    minima.resize(rates.len(), 0.0);
                    sampling::resume_total(&rates, &mut prefix, &mut minima, 0)?
                } else {
                    sampling::total_with_prefix(&rates, &mut prefix)?
                };
            }
            if total == 0.0 {
                // Absorption consumes no RNG and retains the state to the horizon.
                while save < run.time.save_times.len() {
                    if let Some(stop) = poller.poll() {
                        return Ok(stop);
                    }
                    counts.extend_from_slice(&state);
                    save += 1;
                }
                time = run.time.end;
                return Ok(Stop::Completed);
            }
            let wait = sampling::waiting_time(total, sampling::uniform_open(rng))?;
            let next = sampling::event_time(time, wait)?;
            // Strict inequality: an event exactly at a save time comes first.
            while save < run.time.save_times.len() && run.time.save_times[save] < next {
                if let Some(stop) = poller.poll() {
                    return Ok(stop);
                }
                counts.extend_from_slice(&state);
                save += 1;
            }
            if next > run.time.end {
                time = run.time.end;
                return Ok(Stop::Completed);
            }
            let channel = sampling::select_from_prefix(
                &prefix,
                total,
                sampling::uniform_closed_open(rng.next_u64()),
            )?;
            if let Err(error) = network.fire_enabled_reaction(channel, &mut state) {
                if error.code == "CRN-E044" && error.message.contains("safety limit") {
                    record.diagnostic = Some(error);
                    return Ok(Stop::PopulationLimit);
                }
                return Err(error);
            }
            previous = Some(channel);
            time = next;
            record.events += 1;
        }
    })();
    match outcome {
        Ok(stop) => record.termination = stop,
        Err(error) => {
            record.termination = Stop::Failed;
            record.diagnostic = Some(error);
        }
    }
    record.stop_time = time.max(save.checked_sub(1).map_or(0.0, |i| run.time.save_times[i]));
    record.saved_points = save;
    *work = Workspace {
        state,
        rates,
        prefix,
        minima,
    };
    PathResult { record, counts }
}

pub fn execute(prepared: &Prepared, control: &Control) -> Result<Simulation, Diagnostic> {
    execute_inner(prepared, control, false)
}
/// Exact Direct with dependency invalidation; identical sum/selection order and streams.
pub fn execute_cached(prepared: &Prepared, control: &Control) -> Result<Simulation, Diagnostic> {
    execute_inner(prepared, control, true)
}
fn execute_inner(
    prepared: &Prepared,
    control: &Control,
    cached: bool,
) -> Result<Simulation, Diagnostic> {
    // The sealed archive also carries the model, the request and a fixed
    // envelope, so the legacy numerical payload gets what remains.
    let payload = prepared
        .run
        .budget
        .max_output_bytes
        .min(MAX_RESULT_BYTES)
        .saturating_sub(prepared.input_bytes + 16_384);
    execute_paths::<Workspace, (), _>(
        prepared,
        control,
        sampling::SAMPLER,
        Some(payload),
        |id, poller, work| (trajectory(prepared, id, poller, cached, work), ()),
    )
    .map(|(simulation, _)| simulation)
}

/// Shared bounded scheduling and ordered reduction; each algorithm owns its path state.
/// `payload_limit` bounds the serialized legacy result; callers that encode a
/// different result type pass `None` and apply their own byte budget.
pub(crate) fn execute_paths<W: Default + Send, M: Send, F>(
    prepared: &Prepared,
    control: &Control,
    sampler: &str,
    payload_limit: Option<u64>,
    runner: F,
) -> Result<(Simulation, Vec<M>), Diagnostic>
where
    F: Fn(u64, &mut Poller, &mut W) -> (PathResult, M) + Sync,
{
    let start = Instant::now();
    let run = &prepared.run;
    let ensemble = run.ensemble.as_ref().unwrap();
    let n_species = prepared.network.species().len();
    let cells = n_species * run.time.save_times.len();
    let mut mean = vec![0.0; cells];
    let mut m2 = vec![0.0; cells];
    let mut completed = 0;
    let mut records = Vec::with_capacity(ensemble.replicates as usize);
    let mut metrics = Vec::with_capacity(ensemble.replicates as usize);
    let mut kept = Vec::new();
    {
        let mut consume = |id: u64, (path, metric): (PathResult, M)| {
            metrics.push(metric);
            if path.record.termination == Stop::Completed {
                completed += 1;
                for (i, count) in path.counts.iter().enumerate() {
                    let delta = *count as f64 - mean[i];
                    mean[i] += delta / completed as f64;
                    m2[i] += delta * (*count as f64 - mean[i]);
                }
            }
            if id < run.outputs.keep_trajectories {
                let counts = if n_species == 0 {
                    vec![Vec::new(); path.record.saved_points]
                } else {
                    path.counts.chunks(n_species).map(<[i64]>::to_vec).collect()
                };
                kept.push(Trajectory { id, counts });
            }
            records.push(path.record);
        };
        if prepared.workers == 1 {
            let mut poller = Poller::new(control, start, run.budget.wall_time_s);
            let mut work = W::default();
            for id in 0..ensemble.replicates {
                consume(id, runner(id, &mut poller, &mut work));
            }
        } else {
            std::thread::scope(|scope| -> Result<(), Diagnostic> {
                let mut receivers = Vec::new();
                for worker in 0..prepared.workers {
                    let runner = &runner;
                    let (tx, rx) = mpsc::sync_channel(prepared.queue_capacity);
                    receivers.push(rx);
                    std::thread::Builder::new()
                        .name(format!("crn-{worker}"))
                        .spawn_scoped(scope, move || {
                            let mut poller = Poller::new(control, start, run.budget.wall_time_s);
                            let mut work = W::default();
                            for id in (worker as u64..ensemble.replicates).step_by(prepared.workers)
                            {
                                if tx.send(runner(id, &mut poller, &mut work)).is_err() {
                                    break;
                                }
                            }
                        })
                        .map_err(|e| execution_error(format!("cannot start worker: {e}")))?;
                }
                // Per-worker bounded queues prevent an arbitrarily large reorder buffer.
                for id in 0..ensemble.replicates {
                    let path = receivers[id as usize % prepared.workers]
                        .recv()
                        .map_err(|e| {
                            execution_error(format!("worker did not return a trajectory: {e}"))
                        })?;
                    consume(id, path);
                }
                Ok(())
            })?;
        }
    }
    let termination = records
        .iter()
        .map(|r| r.termination)
        .max()
        .unwrap_or(Stop::Failed);
    let assessment = if termination == Stop::Completed {
        Assessment::ExactEnsemble
    } else if termination == Stop::Failed {
        Assessment::FailedChecks
    } else {
        Assessment::Incomplete
    };
    let statistics = if completed == ensemble.replicates && !run.outputs.statistics.is_empty() {
        let chunks = |values: Vec<f64>| values.chunks(n_species).map(<[f64]>::to_vec).collect();
        let want_mean = run.outputs.statistics.contains(&Statistic::Mean);
        let variance: Vec<_> = m2
            .iter()
            .map(|x| x / (completed - 1).max(1) as f64)
            .collect();
        Some(Statistics {
            samples: completed,
            mean_standard_error: (want_mean && completed > 1).then(|| {
                chunks(
                    variance
                        .iter()
                        .map(|v| libm::sqrt(v / completed as f64))
                        .collect(),
                )
            }),
            mean: want_mean.then(|| chunks(mean)),
            sample_variance: run
                .outputs
                .statistics
                .contains(&Statistic::Variance)
                .then(|| chunks(variance)),
        })
    } else {
        None
    };
    let mut simulation = Simulation {
        result: SimulationResult {
            schema: ResultVersion::V02,
            model_sha256: prepared.model_hash.clone(),
            run_sha256: prepared.run_hash.clone(),
            sampler: sampler.into(),
            root_seed: ensemble.root_seed,
            species: prepared.network.species().to_vec(),
            save_times: run.time.save_times.clone(),
            termination,
            assessment,
            requested_replicates: ensemble.replicates,
            completed_replicates: completed,
            events: records.iter().map(|r| r.events).sum(),
            replicates: records,
            trajectories: kept,
            statistics,
            diagnostics: prepared.network.warnings().to_vec(),
        },
        provenance: Provenance {
            crate_version: env!("CARGO_PKG_VERSION").into(),
            rustc: env!("CRN_RUSTC_VERSION").into(),
            target_os: std::env::consts::OS.into(),
            target_arch: std::env::consts::ARCH.into(),
            cargo_lock_sha256: sha256(include_bytes!("../../../Cargo.lock")),
            executable_sha256: None,
            workers: prepared.workers,
            elapsed_seconds: start.elapsed().as_secs_f64(),
            reserved_workspace_bytes: prepared.workspace,
            workspace_cap_bytes: MAX_WORKSPACE_BYTES,
        },
    };
    if let Some(limit) = payload_limit {
        constrain_output(&mut simulation.result, limit)?;
    }
    simulation.provenance.elapsed_seconds = start.elapsed().as_secs_f64();
    Ok((simulation, metrics))
}

/// Whether `value` serializes as compact JSON plus a trailing newline within
/// `limit` bytes. Counts bytes and stops at the limit without buffering.
pub(crate) fn fits_json(value: &impl serde::Serialize, limit: u64) -> bool {
    struct Counter {
        used: u64,
        limit: u64,
    }
    impl std::io::Write for Counter {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() as u64 > self.limit.saturating_sub(self.used) {
                return Err(std::io::Error::other("output limit"));
            }
            self.used += bytes.len() as u64;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(Counter { used: 1, limit }, value).is_ok()
}

pub(crate) fn constrain_output(
    result: &mut SimulationResult,
    limit: u64,
) -> Result<(), Diagnostic> {
    let fits = |value: &SimulationResult| fits_json(value, limit);
    if fits(result) {
        return Ok(());
    }
    result.termination = Stop::OutputLimit;
    result.assessment = Assessment::Incomplete;
    result.statistics = None;
    result.trajectories.clear();
    let diagnostic = Diagnostic::error(
        "CRN-E064",
        "/outputs",
        "serialized numerical payload exceeded the output byte budget",
        "Increase max_output_bytes; trajectories/statistics were removed, replicate accounting was retained.",
    );
    result.diagnostics.push(diagnostic.clone());
    if fits(result) {
        Ok(())
    } else {
        Err(diagnostic)
    }
}

pub fn simulate_bytes(
    model: &[u8],
    run: &[u8],
    control: &Control,
) -> Result<Simulation, Vec<Diagnostic>> {
    execute(&prepare(model, run)?, control).map_err(|e| vec![e])
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::Error;
    use serde_json::json;

    struct Words {
        values: Vec<u64>,
        index: usize,
    }
    impl RngCore for Words {
        fn next_u32(&mut self) -> u32 {
            panic!("unexpected 32-bit draw")
        }
        fn next_u64(&mut self) -> u64 {
            let value = self.values[self.index];
            self.index += 1;
            value
        }
        fn fill_bytes(&mut self, _: &mut [u8]) {
            panic!("unexpected byte draw")
        }
        fn try_fill_bytes(&mut self, _: &mut [u8]) -> Result<(), Error> {
            panic!("unexpected byte draw")
        }
    }
    #[test]
    fn open_uniform_rejects_zero_and_never_returns_one() {
        let mut rng = Words {
            values: vec![0, u64::MAX],
            index: 0,
        };
        assert_eq!(
            sampling::uniform_open(&mut rng),
            1.0 - sampling::UNIFORM_STEP
        );
        assert_eq!(rng.index, 2);
    }
    #[test]
    fn event_exactly_at_horizon_is_included_and_consumes_two_words() {
        let model=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"boundary","convention":"stochastic_counts","omega":1.0,"parameters":{},
            "species":[{"name":"X","initial_count":0}],"reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":1.0}}]})).unwrap();
        let end = libm::log(2.0);
        let run=serde_json::to_vec(&json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":sha256(&model)},"method":{"kind":"ssa","aggregator":"direct"},
            "ensemble":{"replicates":1,"root_seed":0},"time":{"end":end,"save_times":[0.0,end]},"outputs":{"statistics":["mean"],"keep_trajectories":1},
            "budget":{"max_events_per_replicate":1,"threads":1,"max_output_bytes":100000}})).unwrap();
        let prepared = prepare(&model, &run).unwrap();
        let mut rng = Words {
            values: vec![1_u64 << 63, 0],
            index: 0,
        };
        let control = Control::default();
        let path = path_with_rng(
            &prepared.network,
            &prepared.run,
            0,
            &mut Poller::new(&control, Instant::now(), None),
            &mut rng,
        );
        assert_eq!(path.record.termination, Stop::Completed);
        assert_eq!(path.record.events, 1);
        assert_eq!(path.counts, vec![0, 1]);
        assert_eq!(path.record.saved_points, 2);
        assert_eq!(rng.index, 2);
    }
    #[test]
    fn cancel_file_checks_are_rate_limited() {
        // About 20,000 events in one path: polling must not stat the file per event.
        let model=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"polling","convention":"stochastic_counts","omega":1.0,"parameters":{},
            "species":[{"name":"X","initial_count":0}],"reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":1000.0}}]})).unwrap();
        let run=serde_json::to_vec(&json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":sha256(&model)},"method":{"kind":"ssa","aggregator":"direct"},
            "ensemble":{"replicates":1,"root_seed":0},"time":{"end":20.0,"save_times":[0.0,20.0]},"outputs":{"statistics":["mean"],"keep_trajectories":0},
            "budget":{"max_events_per_replicate":1000000,"threads":1,"max_output_bytes":100000}})).unwrap();
        let prepared = prepare(&model, &run).unwrap();
        let control = Control {
            cancel_file: Some(std::env::temp_dir().join(format!(
                "crn-craft-absent-cancel-file-{}",
                std::process::id()
            ))),
            ..Control::default()
        };
        let start = Instant::now();
        let path = path_with_rng(
            &prepared.network,
            &prepared.run,
            0,
            &mut Poller::new(&control, start, None),
            &mut sampling::stream(1, 0),
        );
        let elapsed = start.elapsed();
        assert_eq!(path.record.termination, Stop::Completed);
        assert!(path.record.events > 15_000);
        let checks = control.file_checks.load(Ordering::Relaxed);
        // At most one check per millisecond after the first.
        assert!(
            checks <= 1 + elapsed.as_millis() as u64,
            "{checks} cancel-file checks in {elapsed:?} for {} events",
            path.record.events
        );
    }
    #[test]
    fn the_event_loop_polls_once_per_event_for_small_networks() {
        let model=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"polls","convention":"stochastic_counts","omega":1.0,"parameters":{},
            "species":[{"name":"X","initial_count":0}],"reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":100.0}},
            {"name":"death","reactants":{"X":1},"products":{},"rate":{"kind":"mass_action","coefficient":1.0}}]})).unwrap();
        let run=serde_json::to_vec(&json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":sha256(&model)},"method":{"kind":"ssa","aggregator":"direct"},
            "ensemble":{"replicates":1,"root_seed":3},"time":{"end":10.0,"save_times":[0.0,5.0,10.0]},"outputs":{"statistics":["mean"],"keep_trajectories":0},
            "budget":{"max_events_per_replicate":1000000,"threads":1,"max_output_bytes":100000}})).unwrap();
        let prepared = prepare(&model, &run).unwrap();
        for cached in [false, true] {
            let control = Control::default();
            let path = path_with_rng_mode(
                &prepared.network,
                &prepared.run,
                0,
                &mut Poller::new(&control, Instant::now(), None),
                &mut sampling::stream(3, 0),
                cached,
                &mut Workspace::default(),
            );
            assert_eq!(path.record.termination, Stop::Completed);
            assert!(path.record.events > 1000);
            let polls = control.polls.load(Ordering::Relaxed);
            // One per event, one per saved point, the first and the final loop.
            assert!(
                polls <= path.record.events + 3 + 2,
                "{polls} polls for {} events",
                path.record.events
            );
        }
    }
    #[test]
    fn reused_workspaces_reproduce_fresh_paths() {
        // Paths end in different states and some absorb early, so a stale
        // state, rate or running sum would change a later replicate.
        let model=serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"reuse","convention":"stochastic_counts","omega":1.0,"parameters":{},
            "species":[{"name":"X","initial_count":5},{"name":"Y","initial_count":0}],"reactions":[
            {"name":"convert","reactants":{"X":1},"products":{"Y":1},"rate":{"kind":"mass_action","coefficient":1.5}},
            {"name":"pair","reactants":{"Y":2},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":0.7}},
            {"name":"lose","reactants":{"Y":1},"products":{},"rate":{"kind":"mass_action","coefficient":0.2}}]})).unwrap();
        let run=serde_json::to_vec(&json!({"schema":"crn-craft.run/0.1","model":{"path":"model.json","sha256":sha256(&model)},"method":{"kind":"ssa","aggregator":"direct"},
            "ensemble":{"replicates":40,"root_seed":11},"time":{"end":6.0,"save_times":[0.0,1.0,3.0,6.0]},"outputs":{"statistics":["mean"],"keep_trajectories":40},
            "budget":{"max_events_per_replicate":100000,"threads":1,"max_output_bytes":1000000}})).unwrap();
        let prepared = prepare(&model, &run).unwrap();
        for cached in [false, true] {
            let ensemble = execute_inner(&prepared, &Control::default(), cached).unwrap();
            let control = Control::default();
            for (id, kept) in ensemble.result.trajectories.iter().enumerate() {
                let fresh = path_with_rng_mode(
                    &prepared.network,
                    &prepared.run,
                    id as u64,
                    &mut Poller::new(&control, Instant::now(), None),
                    &mut sampling::stream(11, id as u64),
                    cached,
                    &mut Workspace::default(),
                );
                let rows: Vec<Vec<i64>> = fresh.counts.chunks(2).map(<[i64]>::to_vec).collect();
                assert_eq!(kept.counts, rows, "replicate {id}");
                assert_eq!(ensemble.result.replicates[id].events, fresh.record.events);
            }
        }
    }
    #[test]
    fn encoded_budget_checks_exact_boundary_without_allocating_payload() {
        let mut r = SimulationResult {
            schema: ResultVersion::V02,
            model_sha256: "a".repeat(64),
            run_sha256: "b".repeat(64),
            sampler: sampling::SAMPLER.into(),
            root_seed: 1,
            species: vec!["X".into()],
            save_times: vec![0.0, 1.0],
            termination: Stop::Completed,
            assessment: Assessment::ExactEnsemble,
            requested_replicates: 0,
            completed_replicates: 0,
            events: 0,
            replicates: vec![],
            trajectories: vec![],
            statistics: None,
            diagnostics: vec![],
        };
        let exact = serde_json::to_vec(&r).unwrap().len() as u64 + 1;
        constrain_output(&mut r, exact).unwrap();
        assert_eq!(r.termination, Stop::Completed);
        assert!(constrain_output(&mut r, exact - 1).is_err());
        assert_eq!(r.termination, Stop::OutputLimit);
    }
}
