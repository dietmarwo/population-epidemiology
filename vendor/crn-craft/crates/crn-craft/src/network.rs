//! Immutable, validated indexed representation and checked reaction operations.
use crate::contracts::{Model, Rate, Scalar};
use crate::diagnostic::{Diagnostic, Severity, pointer};
use crate::rates::{CompiledRate, MassKernel};
use crate::validate::validate_model;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

pub const MAX_DEPENDENCY_EDGES: usize = 1_000_000;

/// Sparse entry. Indices refer to the model's species declaration order.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Term {
    pub species: usize,
    pub coefficient: i64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReactionStructure {
    pub name: String,
    pub reactants: Vec<Term>,
    pub products: Vec<Term>,
    /// One stoichiometric matrix column: products minus reactants, zeros omitted.
    pub net_changes: Vec<Term>,
    /// Conservative reads, including reactant gates and custom-rate regulators.
    pub reads: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct CompiledNetwork {
    species: Arc<Vec<String>>,
    initial: Arc<Vec<i64>>,
    limits: Arc<Vec<Option<i64>>>,
    reactions: Arc<Vec<ReactionStructure>>,
    rates: Vec<CompiledRate>,
    dependencies: Arc<Vec<Vec<usize>>>,
    warnings: Vec<Diagnostic>,
    topology_hash: Arc<str>,
}

impl CompiledNetwork {
    pub(crate) fn compiled_rate(&self, reaction: usize) -> &CompiledRate {
        &self.rates[reaction]
    }
    pub fn compile(model: &Model) -> Result<Self, Vec<Diagnostic>> {
        let warnings = validate_model(model);
        if warnings.iter().any(|d| d.severity == Severity::Error) {
            return Err(warnings);
        }
        let indices: BTreeMap<&str, usize> = model
            .species
            .iter()
            .enumerate()
            .map(|(i, s)| (s.name.as_str(), i))
            .collect();
        let resolve = |value: &Scalar| match value {
            Scalar::Literal(x) => *x,
            Scalar::Parameter(name) => model.parameters[name],
        };
        let sparse = |side: &BTreeMap<String, i64>| {
            let mut entries: Vec<_> = side
                .iter()
                .map(|(name, c)| Term {
                    species: indices[name.as_str()],
                    coefficient: *c,
                })
                .collect();
            entries.sort_by_key(|term| term.species);
            entries
        };
        let mut reactions = Vec::with_capacity(model.reactions.len());
        let mut rates = Vec::with_capacity(model.reactions.len());
        let mut readers = vec![Vec::new(); model.species.len()];
        let mut hill_laws = 0;
        for (id, reaction) in model.reactions.iter().enumerate() {
            let reactants = sparse(&reaction.reactants);
            let rate = compile_rate(
                &reaction.rate,
                &reactants,
                &indices,
                &resolve,
                &mut hill_laws,
            );
            let products = sparse(&reaction.products);
            let mut changes = BTreeMap::new();
            for term in &reactants {
                changes.insert(term.species, -term.coefficient);
            }
            for term in &products {
                *changes.entry(term.species).or_insert(0) += term.coefficient;
            }
            let net_changes = changes
                .into_iter()
                .filter(|(_, c)| *c != 0)
                .map(|(species, coefficient)| Term {
                    species,
                    coefficient,
                })
                .collect();
            let mut reads: BTreeSet<usize> = reactants.iter().map(|t| t.species).collect();
            if let Some(species) = rate.extra_read() {
                reads.insert(species);
            }
            for species in &reads {
                readers[*species].push(id);
            }
            reactions.push(ReactionStructure {
                name: reaction.name.clone(),
                reactants,
                products,
                net_changes,
                reads: reads.into_iter().collect(),
            });
            rates.push(rate);
        }
        let mut edge_count = 0;
        let mut dependencies = Vec::with_capacity(reactions.len());
        for (id, reaction) in reactions.iter().enumerate() {
            let mut affected = BTreeSet::new();
            for change in &reaction.net_changes {
                for &reader in &readers[change.species] {
                    if affected.insert(reader) {
                        edge_count += 1;
                        if edge_count > MAX_DEPENDENCY_EDGES {
                            return Err(vec![Diagnostic::error(
                                "CRN-E040",
                                &format!("/reactions/{id}"),
                                "dependency graph exceeds 1000000 edges",
                                "Reduce network coupling/size; analysis does not silently drop edges.",
                            )]);
                        }
                    }
                }
            }
            dependencies.push(affected.into_iter().collect());
        }
        Ok(Self {
            species: Arc::new(model.species.iter().map(|s| s.name.clone()).collect()),
            initial: Arc::new(model.species.iter().map(|s| s.initial_count).collect()),
            limits: Arc::new(model.species.iter().map(|s| s.population_limit).collect()),
            reactions: Arc::new(reactions),
            rates,
            dependencies: Arc::new(dependencies),
            warnings,
            topology_hash: topology_sha256(model).into(),
        })
    }

    /// Rebind values to the template's unchanged topology. Callers in this crate
    /// obtain `model` by replacing only the parameter map of a validated model.
    pub(crate) fn rebind(&self, model: &Model) -> Result<Self, Vec<Diagnostic>> {
        let warnings = validate_model(model);
        if warnings.iter().any(|d| d.severity == Severity::Error) {
            return Err(warnings);
        }
        let indices = self
            .species
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        let resolve = |value: &Scalar| match value {
            Scalar::Literal(x) => *x,
            Scalar::Parameter(name) => model.parameters[name],
        };
        let mut hill_laws = 0;
        let rates = model
            .reactions
            .iter()
            .zip(self.reactions.iter())
            .map(|(r, structure)| {
                compile_rate(
                    &r.rate,
                    &structure.reactants,
                    &indices,
                    &resolve,
                    &mut hill_laws,
                )
            })
            .collect();
        Ok(Self {
            species: Arc::clone(&self.species),
            initial: Arc::clone(&self.initial),
            limits: Arc::clone(&self.limits),
            reactions: Arc::clone(&self.reactions),
            dependencies: Arc::clone(&self.dependencies),
            rates,
            warnings,
            topology_hash: Arc::clone(&self.topology_hash),
        })
    }

    /// Stable identity of structure, literals and fixed populations, with named
    /// parameter values removed. It guards optional generated propensity code.
    pub fn topology_sha256(&self) -> &str {
        &self.topology_hash
    }

    /// Runtime coefficient access for model-specific mass-action kernels.
    #[inline]
    pub fn mass_action_coefficient(&self, reaction: usize) -> Option<f64> {
        match self.rates.get(reaction)? {
            CompiledRate::MassAction(c, _) => Some(*c),
            _ => None,
        }
    }

    #[cfg(test)]
    pub(crate) fn shares_structure(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.reactions, &other.reactions)
            && Arc::ptr_eq(&self.dependencies, &other.dependencies)
            && Arc::ptr_eq(&self.species, &other.species)
    }

    pub fn species(&self) -> &[String] {
        &self.species
    }
    pub fn initial_counts(&self) -> &[i64] {
        &self.initial
    }
    pub fn reactions(&self) -> &[ReactionStructure] {
        &self.reactions
    }
    /// Propensity invalidation only. A scheduler must separately reschedule a fired channel.
    pub fn dependencies(&self) -> &[Vec<usize>] {
        &self.dependencies
    }
    pub fn warnings(&self) -> &[Diagnostic] {
        &self.warnings
    }

    pub fn validate_state(&self, state: &[i64]) -> Result<(), Diagnostic> {
        if state.len() != self.species.len() {
            return Err(Diagnostic::error(
                "CRN-E041",
                "/counts",
                "state length does not match species count",
                "Provide one count per species in declaration order.",
            ));
        }
        for (i, &count) in state.iter().enumerate() {
            if count < 0 || self.limits[i].is_some_and(|limit| count > limit) {
                return Err(Diagnostic::error(
                    "CRN-E041",
                    &pointer("/counts", &self.species[i]),
                    "invalid count or population safety limit exceeded",
                    "Use non-negative i64 counts within each declared population_limit.",
                ));
            }
        }
        Ok(())
    }

    pub fn named_state(&self, counts: &BTreeMap<String, i64>) -> Result<Vec<i64>, Diagnostic> {
        let mut result = Vec::with_capacity(self.species.len());
        for name in self.species.iter() {
            let count = counts.get(name).ok_or_else(|| {
                Diagnostic::error(
                    "CRN-E041",
                    &pointer("/counts", name),
                    "missing species count",
                    "Provide every species exactly once.",
                )
            })?;
            result.push(*count);
        }
        for name in counts.keys() {
            if !self.species.contains(name) {
                return Err(Diagnostic::error(
                    "CRN-E041",
                    &pointer("/counts", name),
                    "unknown species in state",
                    "Use the model's species names.",
                ));
            }
        }
        self.validate_state(&result)?;
        Ok(result)
    }

    pub fn propensities(&self, state: &[i64]) -> Result<Vec<f64>, Diagnostic> {
        self.validate_state(state)?;
        self.rates
            .iter()
            .zip(self.reactions.iter())
            .enumerate()
            .map(|(i, (rate, structure))| rate.evaluate(state, &structure.reactants, i))
            .collect()
    }

    pub fn propensity(&self, reaction: usize, state: &[i64]) -> Result<f64, Diagnostic> {
        self.validate_state(state)?;
        self.propensity_in_valid_state(reaction, state)
    }

    /// For [`crate::ssa::PropensityKernel`] implementations: the engine only
    /// passes states built from `initial_counts` and checked `fire_reaction`
    /// updates, so the O(species) state scan of [`Self::propensity`] is
    /// skipped. The length is still checked; callers must not pass negative
    /// counts.
    pub fn propensity_in_scheduler_state(
        &self,
        reaction: usize,
        state: &[i64],
    ) -> Result<f64, Diagnostic> {
        if state.len() != self.species.len() {
            return self.propensity(reaction, state);
        }
        self.propensity_in_valid_state(reaction, state)
    }

    /// Internal scheduler fast path: state comes from initial_counts or a
    /// successful checked fire_reaction. Avoid O(species) validation per channel.
    pub(crate) fn propensity_in_valid_state(
        &self,
        reaction: usize,
        state: &[i64],
    ) -> Result<f64, Diagnostic> {
        let rate = self
            .rates
            .get(reaction)
            .ok_or_else(|| index_error(reaction))?;
        rate.evaluate(state, &self.reactions[reaction].reactants, reaction)
    }

    /// Apply one enabled reaction atomically. No RNG, elapsed time or scheduling.
    pub fn fire_reaction(&self, reaction: usize, state: &mut [i64]) -> Result<(), Diagnostic> {
        if self.propensity(reaction, state)? == 0.0 {
            return Err(Diagnostic::error(
                "CRN-E043",
                &format!("/reactions/{reaction}"),
                "reaction is disabled at this state",
                "Select a channel with a positive propensity and sufficient reactants.",
            ));
        }
        self.fire_enabled_reaction(reaction, state)
    }

    /// Internal primitive: caller established a positive channel in a valid state.
    /// Updates remain checked and transactional, including all population limits.
    pub(crate) fn fire_enabled_reaction(
        &self,
        reaction: usize,
        state: &mut [i64],
    ) -> Result<(), Diagnostic> {
        let changes = &self
            .reactions
            .get(reaction)
            .ok_or_else(|| index_error(reaction))?
            .net_changes;
        // Validate all updates first: an overflow/limit error cannot partially mutate state.
        for term in changes {
            let next = state[term.species]
                .checked_add(term.coefficient)
                .ok_or_else(|| {
                    Diagnostic::error(
                        "CRN-E044",
                        &pointer("/counts", &self.species[term.species]),
                        "count update overflows i64",
                        "Stop before exceeding the supported count range.",
                    )
                })?;
            // A negative count means the channel was not enabled: a scheduler
            // defect, reported as such and never as a population-limit stop.
            if next < 0 {
                return Err(Diagnostic::error(
                    "CRN-E043",
                    &format!("/reactions/{reaction}"),
                    "reaction would make a count negative; it is disabled at this state",
                    "Select a channel with a positive propensity and sufficient reactants.",
                ));
            }
            if self.limits[term.species].is_some_and(|limit| next > limit) {
                return Err(Diagnostic::error(
                    "CRN-E044",
                    &pointer("/counts", &self.species[term.species]),
                    "reaction would exceed a population safety limit",
                    "Stop this execution or explicitly revise the model's safety limit.",
                ));
            }
        }
        for term in changes {
            state[term.species] += term.coefficient;
        }
        Ok(())
    }
}

/// Fingerprint of fixed topology, literals, initial populations and limits.
/// It does not cover the CRN-craft version: regenerate kernels after upgrading.
pub fn topology_sha256(model: &Model) -> String {
    let mut shape = model.clone();
    for value in shape.parameters.values_mut() {
        *value = 0.0;
    }
    crate::sha256(&serde_json::to_vec(&shape).expect("validated finite model"))
}

fn compile_rate(
    law: &Rate,
    reactants: &[Term],
    indices: &BTreeMap<&str, usize>,
    resolve: &impl Fn(&Scalar) -> f64,
    hill_laws: &mut usize,
) -> CompiledRate {
    match law {
        Rate::MassAction { coefficient } => {
            let c = resolve(coefficient);
            CompiledRate::MassAction(c, MassKernel::compile(c, reactants))
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
        } => CompiledRate::Hill {
            regulator: indices[regulator.as_str()],
            vmax: resolve(vmax),
            basal: resolve(basal),
            k: resolve(k),
            n: resolve(n),
            repression: matches!(law, Rate::HillRepression { .. }),
            memo: {
                *hill_laws += 1;
                (*hill_laws <= crate::rates::MEMO_LAWS).then(Default::default)
            },
        },
        Rate::MichaelisMenten {
            substrate,
            vmax,
            km,
        } => CompiledRate::MichaelisMenten {
            substrate: indices[substrate.as_str()],
            vmax: resolve(vmax),
            km: resolve(km),
        },
    }
}

fn index_error(reaction: usize) -> Diagnostic {
    Diagnostic::error(
        "CRN-E043",
        &format!("/reactions/{reaction}"),
        "reaction index out of range",
        "Use an index from the model's reaction declaration order.",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_unenabled_channel_is_a_disabled_reaction_not_a_population_stop() {
        let model: Model = serde_json::from_value(serde_json::json!({"schema":"crn-craft.model/0.1","name":"gate","convention":"stochastic_counts","omega":1.0,"parameters":{},
            "species":[{"name":"X","initial_count":1,"population_limit":5}],
            "reactions":[{"name":"pair","reactants":{"X":2},"products":{},"rate":{"kind":"mass_action","coefficient":1.0}},
                {"name":"grow","reactants":{},"products":{"X":1},"rate":{"kind":"mass_action","coefficient":1.0}}]})).unwrap();
        let network = CompiledNetwork::compile(&model).unwrap();
        let mut state = vec![1];
        let error = network.fire_enabled_reaction(0, &mut state).unwrap_err();
        assert_eq!(error.code, "CRN-E043");
        assert!(!error.message.contains("safety limit"));
        assert_eq!(state, vec![1]);
        let mut full = vec![5];
        let error = network.fire_enabled_reaction(1, &mut full).unwrap_err();
        assert_eq!(error.code, "CRN-E044");
        assert!(error.message.contains("safety limit"));
        assert_eq!(full, vec![5]);
    }
}
