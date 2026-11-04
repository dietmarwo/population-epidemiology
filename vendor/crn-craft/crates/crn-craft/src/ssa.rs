//! Versioned alternative exact SSA methods. Legacy Direct contracts are unchanged.
mod bounds;
mod selectors;
mod ziggurat;
pub use ziggurat::exponential_from_words;
mod ziggurat_tables;

use crate::{
    contracts,
    diagnostic::Diagnostic,
    sampling,
    simulation::{
        self, Assessment, Control, PathResult, Poller, Replicate, Statistics, Stop, Trajectory,
    },
};
use rand_core::RngCore;
use schemars::JsonSchema;
use selectors::{Composition, RateTree};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, path::Path, sync::Arc};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum RequestVersion {
    #[serde(rename = "crn-craft.ssa-run/0.1")]
    V01,
    /// Adds the optional `method.fast` switch, default false.
    #[serde(rename = "crn-craft.ssa-run/0.2")]
    V02,
    /// Same fields as 0.2; `method.fast` defaults to true. Set it to false for
    /// the checked exact policy.
    #[serde(rename = "crn-craft.ssa-run/0.3")]
    V03,
    #[serde(rename = "crn-craft.ssa-run/0.4")]
    V04,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum ResultVersion {
    #[serde(rename = "crn-craft.ssa-result/0.1")]
    V01,
    /// Answers an ssa-run/0.2 request; `method` may carry `fast: true`.
    #[serde(rename = "crn-craft.ssa-result/0.2")]
    V02,
    /// Answers an ssa-run/0.3 request; `method.fast` is always present.
    #[serde(rename = "crn-craft.ssa-result/0.3")]
    V03,
    #[serde(rename = "crn-craft.ssa-result/0.4")]
    V04,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
    Direct,
    CachedDirect,
    DirectCr,
    SortingDirect,
    RssaCr,
    NextReaction,
}
impl Algorithm {
    fn name(self) -> &'static str {
        match self {
            Self::Direct => "direct",
            Self::CachedDirect => "cached_direct",
            Self::DirectCr => "direct_cr",
            Self::SortingDirect => "sorting_direct",
            Self::RssaCr => "rssa_cr",
            Self::NextReaction => "next_reaction",
        }
    }
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum WaitingTime {
    Inverse,
    Ziggurat,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    Ssa {
        aggregator: Algorithm,
        waiting_time: WaitingTime,
        /// Not allowed in ssa-run/0.1; default false in 0.2 and true in 0.3.
        /// Emulates JumpProcesses.jl: Xoshiro256++ words, running totals
        /// updated by `+= new - old` without loss or 53-bit resolution
        /// checks, and stop polls per batch of events instead of per random
        /// word. Selection stays exact for the running totals; paths differ
        /// from the checked policy and replay only with the same setting.
        #[serde(
            default,
            deserialize_with = "present_bool",
            skip_serializing_if = "Option::is_none"
        )]
        #[schemars(with = "bool")]
        fast: Option<bool>,
    },
}
/// A present `fast` must be a boolean; `null` is rejected, not read as absent.
fn present_bool<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Option<bool>, D::Error> {
    bool::deserialize(d).map(Some)
}
impl Method {
    fn algorithm(self) -> Algorithm {
        let Self::Ssa { aggregator, .. } = self;
        aggregator
    }
    /// The policy in force: an explicit value, else the version's default.
    fn fast_for(self, version: RequestVersion) -> bool {
        let Self::Ssa { fast, .. } = self;
        fast.unwrap_or(version == RequestVersion::V03)
    }
    /// The resolved method as results echo it: 0.1/0.2 results carry `fast`
    /// only when true (unchanged bytes); 0.3 results always carry it.
    fn resolved(self, version: RequestVersion) -> Self {
        let fast = self.fast_for(version);
        let Self::Ssa {
            aggregator,
            waiting_time,
            ..
        } = self;
        Self::Ssa {
            aggregator,
            waiting_time,
            fast: if matches!(version, RequestVersion::V03 | RequestVersion::V04) {
                Some(fast)
            } else {
                fast.then_some(true)
            },
        }
    }
    fn fast(self) -> bool {
        let Self::Ssa { fast, .. } = self;
        fast.unwrap_or(false)
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    pub max_events_per_replicate: u64,
    /// Includes every random word drawn by every rejection/tail loop.
    pub max_sampling_steps_per_replicate: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wall_time_s: Option<f64>,
    #[schemars(range(min = 1, max = 1024))]
    pub threads: u32,
    #[schemars(range(min = 1))]
    pub max_output_bytes: u64,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub schema: RequestVersion,
    pub model: contracts::ModelRef,
    pub method: Method,
    pub time: contracts::TimeGrid,
    pub ensemble: contracts::Ensemble,
    pub outputs: contracts::Outputs,
    pub budget: Budget,
}
impl Request {
    pub(crate) fn legacy(&self) -> contracts::Run {
        contracts::Run {
            schema: contracts::RunVersion::V01,
            model: self.model.clone(),
            method: contracts::Method::Ssa {
                aggregator: contracts::Aggregator::Direct,
            },
            time: self.time.clone(),
            ensemble: Some(self.ensemble.clone()),
            outputs: self.outputs.clone(),
            budget: contracts::Budget {
                max_events_per_replicate: Some(self.budget.max_events_per_replicate),
                wall_time_s: self.budget.wall_time_s,
                threads: self.budget.threads,
                max_output_bytes: self.budget.max_output_bytes,
            },
        }
    }
}
#[derive(Debug, Default, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Work {
    pub rng_words: u64,
    pub proposals: u64,
    pub propensity_evaluations: u64,
    pub bound_updates: u64,
    pub point_fallbacks: u64,
    pub rejected_candidates: u64,
}
impl Work {
    fn add(&mut self, rhs: &Self) {
        self.rng_words += rhs.rng_words;
        self.proposals += rhs.proposals;
        self.propensity_evaluations += rhs.propensity_evaluations;
        self.bound_updates += rhs.bound_updates;
        self.point_fallbacks += rhs.point_fallbacks;
        self.rejected_candidates += rhs.rejected_candidates;
    }
}
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    Completed,
    EventLimit,
    SamplingLimit,
    PopulationLimit,
    WallTimeLimit,
    Cancelled,
    Failed,
    OutputLimit,
}
fn termination(stop: Stop, diagnostic: Option<&Diagnostic>) -> Termination {
    match stop {
        Stop::Completed => Termination::Completed,
        Stop::EventLimit if diagnostic.is_some_and(|d| d.code == "CRN-E067") => {
            Termination::SamplingLimit
        }
        Stop::EventLimit => Termination::EventLimit,
        Stop::PopulationLimit => Termination::PopulationLimit,
        Stop::WallTimeLimit => Termination::WallTimeLimit,
        Stop::Cancelled => Termination::Cancelled,
        Stop::Failed => Termination::Failed,
        Stop::OutputLimit => Termination::OutputLimit,
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub id: u64,
    pub started: bool,
    pub termination: Termination,
    pub events: u64,
    pub stop_time: f64,
    pub saved_points: usize,
    pub diagnostic: Option<Diagnostic>,
    pub work: Work,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SsaResult {
    pub schema: ResultVersion,
    pub method: Method,
    pub model_sha256: String,
    pub run_sha256: String,
    pub sampler: String,
    pub root_seed: u64,
    pub species: Vec<String>,
    pub save_times: Vec<f64>,
    pub termination: Termination,
    pub assessment: Assessment,
    pub requested_replicates: u64,
    pub completed_replicates: u64,
    pub events: u64,
    pub work: Work,
    pub replicates: Vec<Record>,
    pub trajectories: Vec<Trajectory>,
    pub statistics: Option<Statistics>,
    pub diagnostics: Vec<Diagnostic>,
}

pub struct Prepared {
    pub(crate) core: simulation::Prepared,
    pub(crate) request: Request,
    pub(crate) readers: Arc<Vec<Vec<usize>>>,
    bound_groups: Vec<usize>,
}

pub fn prepare(model: &[u8], request_bytes: &[u8]) -> Result<Prepared, Vec<Diagnostic>> {
    let mut request: Request = crate::strict_json::parse(request_bytes).map_err(|e| vec![e])?;
    if request.method.fast_for(request.schema) && request.schema == RequestVersion::V01 {
        return Err(vec![Diagnostic::error(
            "CRN-E020",
            "/method/fast",
            "fast requires crn-craft.ssa-run/0.2 or 0.3",
            "Use crn-craft.ssa-run/0.3 (fast by default), or remove fast for the checked exact policy.",
        )]);
    }
    // From here on `fast` is explicit wherever it applies.
    request.method = request.method.resolved(request.schema);
    let legacy = serde_json::to_vec(&request.legacy())
        .map_err(|e| vec![simulation::execution_error(e.to_string())])?;
    let core = simulation::prepare(model, &legacy)?;
    finish_prepare(core, request, request_bytes, model.len(), None)
}

pub(crate) fn finish_prepare(
    mut core: simulation::Prepared,
    request: Request,
    request_bytes: &[u8],
    model_len: usize,
    shared_readers: Option<Arc<Vec<Vec<usize>>>>,
) -> Result<Prepared, Vec<Diagnostic>> {
    if request.method.algorithm() == Algorithm::NextReaction
        && (request.schema != RequestVersion::V04
            || request.method.fast()
            || !matches!(
                request.method,
                Method::Ssa {
                    waiting_time: WaitingTime::Inverse,
                    ..
                }
            ))
    {
        return Err(vec![simulation::execution_error(
            "next_reaction requires ssa-run/0.4, inverse waiting time, and checked policy",
        )]);
    }
    // Additional selector and bounded-range buffers, conservatively reserved.
    let reactions = core.network.reactions().len() as u64;
    let species = core.network.species().len() as u64;
    let extra =
        (core.workers as u128) * (reactions as u128 * 384 + species as u128 * 128 + 200_000);
    let reads: usize = core.network.reactions().iter().map(|r| r.reads.len()).sum();
    let workspace = core.workspace as u128 + extra + reads as u128 * 8 + species as u128 * 24;
    let max_work = (u128::from(request.budget.max_events_per_replicate) + 1)
        * (u128::from(reactions) * 260 + u128::from(species) + 1)
        * u128::from(request.ensemble.replicates);
    let sampling_work = (u128::from(request.budget.max_sampling_steps_per_replicate) + 1)
        * u128::from(request.ensemble.replicates);
    if workspace > u128::from(simulation::MAX_WORKSPACE_BYTES)
        || max_work + sampling_work > u128::from(u64::MAX)
    {
        return Err(vec![simulation::execution_error(
            "alternative SSA workspace or work accounting exceeds limits",
        )]);
    }
    core.workspace = workspace as u64;
    core.run_hash = crate::sha256(request_bytes);
    core.input_bytes = (model_len + request_bytes.len()) as u64;
    let readers = shared_readers.unwrap_or_else(|| {
        let mut readers = vec![Vec::new(); species as usize];
        for (rx, structure) in core.network.reactions().iter().enumerate() {
            for &s in &structure.reads {
                readers[s].push(rx);
            }
        }
        Arc::new(readers)
    });
    let mut group_keys = BTreeMap::new();
    let mut bound_groups = vec![usize::MAX; reactions as usize];
    for (rx, slot) in bound_groups.iter_mut().enumerate() {
        if let Some(key) = core.network.compiled_rate(rx).bounds_key() {
            let next = group_keys.len();
            *slot = *group_keys.entry(key).or_insert(next);
        }
    }
    Ok(Prepared {
        core,
        request,
        readers,
        bound_groups,
    })
}

fn numeric(message: impl Into<String>) -> Diagnostic {
    Diagnostic::error(
        "CRN-E061",
        "/simulation",
        message,
        "Rescale the rates; the exact sampler does not silently discard probability.",
    )
}
enum Bail {
    Stopped(Stop),
    SamplingLimit,
    Numerical(Diagnostic),
}
impl From<Diagnostic> for Bail {
    fn from(value: Diagnostic) -> Self {
        Self::Numerical(value)
    }
}
struct Draws<'a, 'b> {
    limit: u64,
    work: &'a mut Work,
    poller: &'a mut Poller<'b>,
    /// `fast` polls for stops per batch of event-loop iterations and on
    /// rejection retries, not before every random word.
    fast: bool,
}
impl Draws<'_, '_> {
    fn poll(&mut self) -> Result<(), Bail> {
        if let Some(s) = self.poller.poll() {
            Err(Bail::Stopped(s))
        } else {
            Ok(())
        }
    }
    /// Before redrawing in a rejection loop (exact runs already polled per word).
    fn retry(&mut self) -> Result<(), Bail> {
        if self.fast { self.poll() } else { Ok(()) }
    }
    fn word(&mut self, rng: &mut impl RngCore) -> Result<u64, Bail> {
        if !self.fast {
            self.poll()?;
        }
        if self.work.rng_words >= self.limit {
            return Err(Bail::SamplingLimit);
        }
        self.work.rng_words += 1;
        Ok(rng.next_u64())
    }
    fn closed(&mut self, rng: &mut impl RngCore) -> Result<f64, Bail> {
        Ok(sampling::uniform_closed_open(self.word(rng)?))
    }
    fn open(&mut self, rng: &mut impl RngCore) -> Result<f64, Bail> {
        loop {
            let u = self.closed(rng)?;
            if u > 0.0 {
                return Ok(u);
            }
            self.retry()?;
        }
    }
    fn index(&mut self, rng: &mut impl RngCore, n: usize) -> Result<usize, Bail> {
        let n = n as u64;
        let cutoff = n.wrapping_neg() % n;
        loop {
            let product = u128::from(self.word(rng)?) * u128::from(n);
            if product as u64 >= cutoff {
                return Ok((product >> 64) as usize);
            }
            self.retry()?;
        }
    }
}

/// Model-specific propensity evaluators must match the checked reference in
/// every valid state, including gates, rounded arithmetic and diagnostics.
/// Generated kernels are guarded by the fixed-topology fingerprint on execute.
pub trait PropensityKernel: Sync {
    fn topology_sha256(&self) -> &str;
    fn evaluate(
        &self,
        network: &crate::CompiledNetwork,
        reaction: usize,
        state: &[i64],
    ) -> Result<f64, Diagnostic>;
}
struct DynamicKernel;
impl PropensityKernel for DynamicKernel {
    fn topology_sha256(&self) -> &str {
        ""
    }
    #[inline]
    fn evaluate(
        &self,
        network: &crate::CompiledNetwork,
        reaction: usize,
        state: &[i64],
    ) -> Result<f64, Diagnostic> {
        network.propensity_in_valid_state(reaction, state)
    }
}

#[derive(Default)]
struct Workspace {
    state: Vec<i64>,
    rates: Vec<f64>,
    prefix: Vec<f64>,
    minima: Vec<f64>,
    tree: RateTree,
    composition: Composition,
    order: Vec<usize>,
    bounds: bounds::Bounds,
}

fn path<const MODE: u8, const ZIG: bool, const FAST: bool, K: PropensityKernel>(
    prepared: &Prepared,
    id: u64,
    poller: &mut Poller,
    space: &mut Workspace,
    kernel: &K,
) -> (PathResult, Work) {
    let seed = prepared.core.run.ensemble.as_ref().unwrap().root_seed;
    if FAST {
        run_path::<MODE, ZIG, FAST, _, _>(
            prepared,
            id,
            poller,
            space,
            sampling::fast_stream(seed, id),
            kernel,
        )
    } else {
        run_path::<MODE, ZIG, FAST, _, _>(
            prepared,
            id,
            poller,
            space,
            sampling::stream(seed, id),
            kernel,
        )
    }
}

/// Stop polls in `fast` runs: one per this many event-loop iterations.
const FAST_POLL_STRIDE: u64 = 256;

fn run_path<const MODE: u8, const ZIG: bool, const FAST: bool, R: RngCore, K: PropensityKernel>(
    prepared: &Prepared,
    id: u64,
    poller: &mut Poller,
    space: &mut Workspace,
    mut rng: R,
    kernel: &K,
) -> (PathResult, Work) {
    let network = &prepared.core.network;
    let run = &prepared.core.run;
    let mut record = Replicate {
        id,
        started: false,
        termination: Stop::Completed,
        events: 0,
        stop_time: 0.0,
        saved_points: 0,
        diagnostic: None,
    };
    let mut work = Work::default();
    let mut draws = Draws {
        limit: prepared.request.budget.max_sampling_steps_per_replicate,
        work: &mut work,
        poller,
        fast: FAST,
    };
    let state = &mut space.state;
    state.clear();
    state.extend_from_slice(network.initial_counts());
    let n = network.reactions().len();
    space.rates.resize(n, 0.0);
    space.prefix.resize(n, 0.0);
    space.minima.resize(n, 0.0);
    if MODE >= 2 && MODE != 5 {
        space.tree.reset_fast(n, FAST);
        if MODE == 2 || MODE == 4 {
            space.composition.reset(n);
        }
        if MODE == 3 {
            space.order.clear();
            space.order.extend(0..n);
        }
    }
    if MODE == 4 {
        space.bounds.reset(network, state);
    }
    let mut counts = Vec::with_capacity(run.time.save_times.len() * state.len());
    let mut time = 0.0;
    let mut save = 0;
    let mut previous: Option<usize> = None;
    let mut total = 0.0;
    let mut bounds_initialized = false;
    let mut absorption_checked = false;
    let mut absorbed = false;
    let mut iterations: u64 = 0;
    let mut clocks = crate::next_reaction::Clocks::default();
    if MODE == 5 {
        clocks.reset(n);
    }
    let outcome = (|| -> Result<Stop, Bail> {
        draws.poll()?;
        record.started = true;
        loop {
            while save < run.time.save_times.len() && run.time.save_times[save] <= time {
                draws.poll()?;
                counts.extend_from_slice(state);
                save += 1;
            }
            if time == run.time.end {
                return Ok(Stop::Completed);
            }
            if !FAST || iterations % FAST_POLL_STRIDE == 0 {
                draws.poll()?;
            }
            iterations += 1;
            if record.events >= run.budget.max_events_per_replicate.unwrap() {
                return Ok(Stop::EventLimit);
            }
            let (next, channel) = if MODE == 5 {
                let mut update = |rx, fired| -> Result<(), Bail> {
                    draws.poll()?;
                    let a = kernel.evaluate(network, rx, state)?;
                    draws.work.propensity_evaluations += 1;
                    let fresh = if clocks.needs_draw(rx, a, fired) {
                        Some(-libm::log(draws.open(&mut rng)?))
                    } else {
                        None
                    };
                    if !clocks.update(rx, a, time, fresh) {
                        return Err(numeric("next reaction clock failed to advance").into());
                    }
                    Ok(())
                };
                if let Some(rx) = previous {
                    update(rx, true)?;
                    for &other in &network.dependencies()[rx] {
                        if other != rx {
                            update(other, false)?;
                        }
                    }
                } else {
                    for rx in 0..n {
                        update(rx, false)?;
                    }
                }
                let candidate = clocks.next();
                let until = candidate.map_or(f64::INFINITY, |(t, _)| t);
                while save < run.time.save_times.len() && run.time.save_times[save] < until {
                    draws.poll()?;
                    counts.extend_from_slice(state);
                    save += 1;
                }
                if until > run.time.end {
                    time = run.time.end;
                    return Ok(Stop::Completed);
                }
                candidate.expect("finite clock")
            } else {
                if MODE == 4 {
                    if !bounds_initialized || previous.is_some() {
                        space.bounds.refresh(
                            network,
                            &prepared.readers,
                            &prepared.bound_groups,
                            state,
                            previous.take(),
                            &mut space.tree,
                            &mut space.composition,
                            &mut draws,
                        )?;
                        bounds_initialized = true;
                    }
                    if !absorption_checked && space.bounds.positive_lower == 0 {
                        absorbed = true;
                        for rx in 0..n {
                            draws.poll()?;
                            draws.work.propensity_evaluations += 1;
                            if kernel.evaluate(network, rx, state)? > 0.0 {
                                absorbed = false;
                                break;
                            }
                        }
                    }
                    absorption_checked = true;
                } else if MODE > 0 && previous.is_some() {
                    let affected = &network.dependencies()[previous.unwrap()];
                    for (k, &rx) in affected.iter().enumerate() {
                        if k % 64 == 0 && (!FAST || k > 0) {
                            draws.poll()?;
                        }
                        let a = kernel.evaluate(network, rx, state)?;
                        draws.work.propensity_evaluations += 1;
                        space.rates[rx] = a;
                        if MODE >= 2 {
                            space
                                .tree
                                .set(rx, a, if a > 0.0 { a } else { f64::INFINITY });
                            if MODE == 2 {
                                space.composition.update(rx, a)?;
                            }
                        }
                    }
                    if MODE == 1 {
                        if let Some(&first) = affected.first() {
                            total = if FAST {
                                sampling::fast_prefix(&space.rates, &mut space.prefix, first)?
                            } else {
                                sampling::resume_total(
                                    &space.rates,
                                    &mut space.prefix,
                                    &mut space.minima,
                                    first,
                                )?
                            };
                        }
                    }
                } else {
                    for rx in 0..n {
                        if rx % 64 == 0 {
                            draws.poll()?;
                        }
                        let a = kernel.evaluate(network, rx, state)?;
                        draws.work.propensity_evaluations += 1;
                        space.rates[rx] = a;
                        if MODE >= 2 {
                            space
                                .tree
                                .set(rx, a, if a > 0.0 { a } else { f64::INFINITY });
                            if MODE == 2 {
                                space.composition.update(rx, a)?;
                            }
                        }
                    }
                    if FAST && MODE <= 1 {
                        total = sampling::fast_prefix(&space.rates, &mut space.prefix, 0)?;
                    } else if MODE == 0 {
                        total = sampling::total_with_prefix(&space.rates, &mut space.prefix)?;
                    } else if MODE == 1 {
                        total = sampling::resume_total(
                            &space.rates,
                            &mut space.prefix,
                            &mut space.minima,
                            0,
                        )?;
                    }
                }
                if MODE >= 2 {
                    total = if MODE == 4 && absorbed {
                        0.0
                    } else {
                        space.tree.total()?
                    };
                }
                if total == 0.0 {
                    while save < run.time.save_times.len() {
                        draws.poll()?;
                        counts.extend_from_slice(state);
                        save += 1;
                    }
                    time = run.time.end;
                    return Ok(Stop::Completed);
                }
                let wait = if ZIG {
                    ziggurat::wait(total, &mut rng, &mut draws)?
                } else {
                    sampling::waiting_time(total, draws.open(&mut rng)?)?
                };
                let next = sampling::event_time(time, wait)?;
                while save < run.time.save_times.len() && run.time.save_times[save] < next {
                    draws.poll()?;
                    counts.extend_from_slice(state);
                    save += 1;
                }
                if next > run.time.end {
                    time = run.time.end;
                    return Ok(Stop::Completed);
                }
                let channel = if MODE == 2 || MODE == 4 {
                    space
                        .composition
                        .sample(&space.tree, &mut rng, &mut draws)?
                } else if MODE == 3 {
                    draws.work.proposals += 1;
                    selectors::sorting_sample(
                        &mut space.order,
                        &space.rates,
                        total,
                        &mut rng,
                        &mut draws,
                    )?
                } else {
                    draws.work.proposals += 1;
                    sampling::select_from_prefix(&space.prefix, total, draws.closed(&mut rng)?)?
                };
                (next, channel)
            };
            if MODE == 4 {
                let upper = space.tree.get(channel);
                let threshold = draws.closed(&mut rng)? * upper;
                if threshold >= upper {
                    return Err(numeric("RSSA acceptance rounded to its bound").into());
                }
                if threshold >= space.bounds.lower[channel] {
                    let a = kernel.evaluate(network, channel, state)?;
                    draws.work.propensity_evaluations += 1;
                    if a < space.bounds.lower[channel] || a > upper {
                        return Err(numeric(
                            "RSSA propensity escaped its certified integer-box bounds",
                        )
                        .into());
                    }
                    if threshold >= a {
                        time = next;
                        draws.work.rejected_candidates += 1;
                        continue;
                    }
                }
            }
            if let Err(error) = network.fire_enabled_reaction(channel, state) {
                if error.code == "CRN-E044" && error.message.contains("safety limit") {
                    record.diagnostic = Some(error);
                    return Ok(Stop::PopulationLimit);
                }
                return Err(error.into());
            }
            previous = Some(channel);
            absorption_checked = false;
            time = next;
            record.events += 1;
        }
    })();
    match outcome {
        Ok(stop) | Err(Bail::Stopped(stop)) => record.termination = stop,
        Err(Bail::SamplingLimit) => {
            record.termination = Stop::EventLimit;
            record.diagnostic = Some(Diagnostic::error(
                "CRN-E067",
                "/budget/max_sampling_steps_per_replicate",
                "sampling work limit reached",
                "Increase the explicit sampling-step budget; incomplete ensembles withhold statistics.",
            ));
        }
        Err(Bail::Numerical(error)) => {
            record.termination = Stop::Failed;
            record.diagnostic = Some(error);
        }
    }
    record.stop_time = time.max(save.checked_sub(1).map_or(0.0, |i| run.time.save_times[i]));
    record.saved_points = save;
    (PathResult { record, counts }, work)
}

pub fn execute(prepared: &Prepared, control: &Control) -> Result<SsaResult, Diagnostic> {
    execute_inner(prepared, control, &DynamicKernel)
}

/// Execute an optional trusted Rust kernel against its exact topology. As with
/// custom Rust simulation code, callers must establish evaluator correctness;
/// this fingerprint prevents accidental use with a different network.
pub fn execute_with_kernel<K: PropensityKernel>(
    prepared: &Prepared,
    control: &Control,
    kernel: &K,
) -> Result<SsaResult, Diagnostic> {
    if kernel.topology_sha256() != prepared.core.network.topology_sha256() {
        return Err(Diagnostic::error(
            "CRN-E020",
            "/kernel",
            "generated kernel topology mismatch",
            "Regenerate the kernel from this parameter template's model.",
        ));
    }
    execute_inner(prepared, control, kernel)
}

fn execute_inner<K: PropensityKernel>(
    prepared: &Prepared,
    control: &Control,
    kernel: &K,
) -> Result<SsaResult, Diagnostic> {
    let algorithm = prepared.request.method.algorithm();
    let Method::Ssa { waiting_time, .. } = prepared.request.method;
    // Resolved at preparation: explicit for 0.3, present only when true before.
    let fast = prepared.request.method.fast();
    let wait = match waiting_time {
        WaitingTime::Inverse => "inverse",
        WaitingTime::Ziggurat => "ziggurat256-julia",
    };
    let sampler = if fast {
        format!("xoshiro256pp-u53/{}/{wait}/fast-1", algorithm.name())
    } else {
        format!(
            "chacha8-u53/{}/{wait}/{}",
            algorithm.name(),
            // Version 2: composition draws member and acceptance from one word.
            if matches!(algorithm, Algorithm::DirectCr | Algorithm::RssaCr) {
                2
            } else {
                1
            }
        )
    };
    macro_rules! run {
        ($mode:literal, $zig:literal, $id:expr, $poller:expr, $work:expr) => {
            if fast {
                path::<$mode, $zig, true, _>(prepared, $id, $poller, $work, kernel)
            } else {
                path::<$mode, $zig, false, _>(prepared, $id, $poller, $work, kernel)
            }
        };
    }
    let (simulation, metrics) = simulation::execute_paths::<Workspace, Work, _>(
        &prepared.core,
        control,
        &sampler,
        None,
        |id, poller, work| match (algorithm, waiting_time) {
            (Algorithm::Direct, WaitingTime::Inverse) => run!(0, false, id, poller, work),
            (Algorithm::CachedDirect, WaitingTime::Inverse) => run!(1, false, id, poller, work),
            (Algorithm::DirectCr, WaitingTime::Inverse) => run!(2, false, id, poller, work),
            (Algorithm::SortingDirect, WaitingTime::Inverse) => run!(3, false, id, poller, work),
            (Algorithm::RssaCr, WaitingTime::Inverse) => run!(4, false, id, poller, work),
            (Algorithm::Direct, WaitingTime::Ziggurat) => run!(0, true, id, poller, work),
            (Algorithm::CachedDirect, WaitingTime::Ziggurat) => run!(1, true, id, poller, work),
            (Algorithm::DirectCr, WaitingTime::Ziggurat) => run!(2, true, id, poller, work),
            (Algorithm::SortingDirect, WaitingTime::Ziggurat) => run!(3, true, id, poller, work),
            (Algorithm::RssaCr, WaitingTime::Ziggurat) => run!(4, true, id, poller, work),
            (Algorithm::NextReaction, _) => {
                path::<5, false, false, _>(prepared, id, poller, work, kernel)
            }
        },
    )?;
    let r = simulation.result;
    let mut work = Work::default();
    let replicates: Vec<_> = r
        .replicates
        .into_iter()
        .zip(metrics)
        .map(|(r, w)| {
            work.add(&w);
            Record {
                id: r.id,
                started: r.started,
                termination: termination(r.termination, r.diagnostic.as_ref()),
                events: r.events,
                stop_time: r.stop_time,
                saved_points: r.saved_points,
                diagnostic: r.diagnostic,
                work: w,
            }
        })
        .collect();
    let mut stop = replicates
        .iter()
        .map(|r| r.termination)
        .max()
        .unwrap_or(Termination::Failed);
    if r.termination == Stop::OutputLimit {
        stop = Termination::OutputLimit;
    }
    let mut result = SsaResult {
        schema: match prepared.request.schema {
            RequestVersion::V01 => ResultVersion::V01,
            RequestVersion::V02 => ResultVersion::V02,
            RequestVersion::V03 => ResultVersion::V03,
            RequestVersion::V04 => ResultVersion::V04,
        },
        method: prepared.request.method,
        model_sha256: r.model_sha256,
        run_sha256: r.run_sha256,
        sampler: r.sampler,
        root_seed: r.root_seed,
        species: r.species,
        save_times: r.save_times,
        termination: stop,
        assessment: r.assessment,
        requested_replicates: r.requested_replicates,
        completed_replicates: r.completed_replicates,
        events: r.events,
        work,
        replicates,
        trajectories: r.trajectories,
        statistics: r.statistics,
        diagnostics: r.diagnostics,
    };
    let limit = prepared
        .request
        .budget
        .max_output_bytes
        .min(simulation::MAX_RESULT_BYTES);
    // The budget covers this compact result JSON and its trailing newline.
    if !simulation::fits_json(&result, limit) {
        result.termination = Termination::OutputLimit;
        result.assessment = Assessment::Incomplete;
        result.statistics = None;
        result.trajectories.clear();
        result.diagnostics.push(Diagnostic::error(
            "CRN-E064",
            "/outputs",
            "alternative SSA output exceeds its byte budget",
            "Increase max_output_bytes; statistics and trajectories have been removed.",
        ));
        if !simulation::fits_json(&result, limit) {
            return Err(simulation::execution_error(
                "SSA evidence exceeds output budget",
            ));
        }
    }
    Ok(result)
}
pub fn solve(
    model: &[u8],
    request: &[u8],
    control: &Control,
) -> Result<SsaResult, Vec<Diagnostic>> {
    execute(&prepare(model, request)?, control).map_err(|e| vec![e])
}
pub fn solve_files(
    model_path: &Path,
    request_path: &Path,
    control: &Control,
) -> Result<SsaResult, Vec<Diagnostic>> {
    let model = crate::read_document(model_path).map_err(|e| vec![e])?;
    let bytes = crate::read_document(request_path).map_err(|e| vec![e])?;
    let request: Request = crate::strict_json::parse(&bytes).map_err(|e| vec![e])?;
    let parent = request_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let resolved = (
        parent.canonicalize(),
        parent.join(&request.model.path).canonicalize(),
        model_path.canonicalize(),
    );
    if !crate::validate::valid_relative_path(&request.model.path)
        || !matches!(resolved,(Ok(ref base),Ok(ref reference),Ok(ref supplied)) if reference.starts_with(base) && reference==supplied)
    {
        return Err(vec![Diagnostic::error(
            "CRN-E022",
            "/model/path",
            "SSA request does not identify the supplied model below its directory",
            "Use the referenced model within the request directory.",
        )]);
    }
    solve(&model, &bytes, control)
}
