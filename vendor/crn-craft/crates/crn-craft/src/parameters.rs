//! Reusable fixed-topology parameter bindings. No candidate mutates its template.
use crate::{contracts::Model, diagnostic::Diagnostic, simulation, ssa};
use std::{collections::BTreeMap, sync::Arc};

/// Validated topology and run settings. Values are ordered by `parameter_names`.
/// Binding still validates the candidate model and serializes auditable inputs;
/// the compiled stoichiometry, dependencies and reverse readers are shared.
pub struct ParameterTemplate {
    model: Model,
    prepared: ssa::Prepared,
    names: Vec<String>,
}

/// Independent execution inputs with full snapshots for replay and attribution.
pub struct BoundParameters {
    pub prepared: ssa::Prepared,
    pub model_bytes: Vec<u8>,
    pub request_bytes: Vec<u8>,
}

fn invalid(message: impl Into<String>) -> Vec<Diagnostic> {
    vec![Diagnostic::error(
        "CRN-E015",
        "/parameters",
        message,
        "Supply exactly one finite value per declared parameter in parameter_names order; satisfy each rate law's domain.",
    )]
}

/// Wrong number of values or a key set different from the declared names.
fn mismatch(message: impl Into<String>) -> Vec<Diagnostic> {
    vec![Diagnostic::error(
        "CRN-E014",
        "/parameters",
        message,
        "Supply exactly one value per declared parameter, in parameter_names order (sorted by name) or by name.",
    )]
}

impl ParameterTemplate {
    pub fn prepare(model: &[u8], request: &[u8]) -> Result<Self, Vec<Diagnostic>> {
        let prepared = ssa::prepare(model, request)?;
        let model = crate::parse_model(model).map_err(|e| vec![e])?;
        let names = model.parameters.keys().cloned().collect();
        Ok(Self {
            model,
            prepared,
            names,
        })
    }

    pub fn parameter_names(&self) -> &[String] {
        &self.names
    }
    pub fn default_values(&self) -> Vec<f64> {
        self.model.parameters.values().copied().collect()
    }
    pub fn topology_sha256(&self) -> &str {
        self.prepared.core.network.topology_sha256()
    }

    /// Values include every declared parameter, even unused ones. Seed overrides
    /// the run root seed only; it does not alter the topology or rate bindings.
    pub fn bind(
        &self,
        values: &[f64],
        seed: Option<u64>,
    ) -> Result<BoundParameters, Vec<Diagnostic>> {
        if values.len() != self.names.len() {
            return Err(mismatch(format!(
                "expected {} values, got {}",
                self.names.len(),
                values.len()
            )));
        }
        let mut model = self.model.clone();
        for (name, value) in self.names.iter().zip(values) {
            model.parameters.insert(name.clone(), *value);
        }
        // Rebuild every resolved law and its memo, never the topology/graph.
        let network = self.prepared.core.network.rebind(&model)?;
        let model_bytes = serde_json::to_vec(&model).map_err(|e| invalid(e.to_string()))?;
        let mut request = self.prepared.request.clone();
        request.model.sha256 = crate::sha256(&model_bytes);
        if let Some(seed) = seed {
            request.ensemble.root_seed = seed;
        }
        let request_bytes = serde_json::to_vec(&request).map_err(|e| invalid(e.to_string()))?;
        if model_bytes.len() > crate::strict_json::MAX_DOCUMENT_BYTES
            || request_bytes.len() > crate::strict_json::MAX_DOCUMENT_BYTES
        {
            return Err(vec![Diagnostic::error(
                "CRN-E004",
                "",
                "document exceeds 16 MiB",
                "Split or reduce the model.",
            )]);
        }
        let legacy = request.legacy();
        let legacy_bytes = serde_json::to_vec(&legacy).map_err(|e| invalid(e.to_string()))?;
        let core = simulation::prepare_compiled(
            network,
            legacy,
            request.model.sha256.clone(),
            crate::sha256(&legacy_bytes),
            model_bytes.len(),
            legacy_bytes.len(),
        )?;
        let prepared = ssa::finish_prepare(
            core,
            request,
            &request_bytes,
            model_bytes.len(),
            Some(Arc::clone(&self.prepared.readers)),
        )?;
        Ok(BoundParameters {
            prepared,
            model_bytes,
            request_bytes,
        })
    }

    /// Named binding requires exactly the declared key set; no ignored values.
    pub fn bind_named(
        &self,
        values: &BTreeMap<String, f64>,
        seed: Option<u64>,
    ) -> Result<BoundParameters, Vec<Diagnostic>> {
        if values.len() != self.names.len() || !values.keys().eq(self.names.iter()) {
            let missing: Vec<_> = self
                .names
                .iter()
                .filter(|n| !values.contains_key(*n))
                .collect();
            let extra: Vec<_> = values.keys().filter(|n| !self.names.contains(n)).collect();
            return Err(mismatch(format!(
                "parameter names do not match template: missing {missing:?}, unknown {extra:?}"
            )));
        }
        self.bind(&values.values().copied().collect::<Vec<_>>(), seed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn model(custom: bool) -> Vec<u8> {
        let laws = if custom {
            vec![
                json!({"kind":"hill_repression","regulator":"X","vmax":"a","basal":0.1,"k":"k","n":2.0}),
                json!({"kind":"hill_repression","regulator":"X","vmax":"b","basal":0.1,"k":"k","n":2.0}),
            ]
        } else {
            vec![
                json!({"kind":"mass_action","coefficient":"a"}),
                json!({"kind":"mass_action","coefficient":"b"}),
            ]
        };
        serde_json::to_vec(&json!({"schema":"crn-craft.model/0.1","name":"parameters","convention":"stochastic_counts","omega":1.0,
            "parameters":{"a":2.0,"b":2.0,"k":3.0},"species":[{"name":"X","initial_count":4}],
            "reactions":[{"name":"birth","reactants":{},"products":{"X":1},"rate":laws[0]},
                {"name":"death","reactants":{"X":1},"products":{},"rate":laws[1]}]})).unwrap()
    }
    fn request(model: &[u8], method: &str, wait: &str, fast: bool, threads: u32) -> Vec<u8> {
        serde_json::to_vec(&json!({"schema":"crn-craft.ssa-run/0.2","model":{"path":"model.json","sha256":crate::sha256(model)},
            "method":{"kind":"ssa","aggregator":method,"waiting_time":wait,"fast":fast},
            "ensemble":{"replicates":4,"root_seed":4242},"time":{"end":2.0,"save_times":[0.0,1.0,2.0]},
            "outputs":{"statistics":["mean","variance"],"keep_trajectories":4},
            "budget":{"max_events_per_replicate":10000,"max_sampling_steps_per_replicate":1000000,"threads":threads,"max_output_bytes":1000000}})).unwrap()
    }
    #[test]
    fn shares_structure_and_matches_fresh_preparation_all_modes() {
        for custom in [false, true] {
            let model = model(custom);
            for method in [
                "direct",
                "cached_direct",
                "direct_cr",
                "sorting_direct",
                "rssa_cr",
            ] {
                for wait in ["inverse", "ziggurat"] {
                    for fast in [false, true] {
                        let template = ParameterTemplate::prepare(
                            &model,
                            &request(&model, method, wait, fast, 1),
                        )
                        .unwrap();
                        // Populate a template's Hill memo; the next binding must not inherit it.
                        template.prepared.core.network.propensities(&[4]).unwrap();
                        for values in [[0.0, 2.0, 3.0], [6.0, 0.7, 0.5], [2.0, 2.0, 3.0]] {
                            let bound = template.bind(&values, Some(7301)).unwrap();
                            assert!(
                                template
                                    .prepared
                                    .core
                                    .network
                                    .shares_structure(&bound.prepared.core.network)
                            );
                            assert!(Arc::ptr_eq(
                                &template.prepared.readers,
                                &bound.prepared.readers
                            ));
                            let fresh =
                                ssa::prepare(&bound.model_bytes, &bound.request_bytes).unwrap();
                            let control = simulation::Control::default();
                            let a = ssa::execute(&bound.prepared, &control).unwrap();
                            let b = ssa::execute(&fresh, &control).unwrap();
                            assert_eq!(
                                serde_json::to_vec(&a).unwrap(),
                                serde_json::to_vec(&b).unwrap()
                            );
                            assert_eq!(
                                bound.prepared.core.model_hash,
                                crate::sha256(&bound.model_bytes)
                            );
                            assert_eq!(
                                bound.prepared.core.run_hash,
                                crate::sha256(&bound.request_bytes)
                            );
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn domains_named_keys_and_failed_bindings_do_not_mutate_template() {
        let model = model(true);
        let template =
            ParameterTemplate::prepare(&model, &request(&model, "rssa_cr", "ziggurat", true, 1))
                .unwrap();
        for values in [
            vec![],
            vec![1.0],
            vec![1.0, 2.0, 0.0],
            vec![-1.0, 2.0, 1.0],
            vec![f64::NAN, 2.0, 1.0],
            vec![1.0, f64::INFINITY, 1.0],
        ] {
            assert!(template.bind(&values, None).is_err());
        }
        assert!(
            template
                .bind_named(&BTreeMap::from([("missing".into(), 1.0)]), None)
                .is_err()
        );
        assert_eq!(template.default_values(), vec![2.0, 2.0, 3.0]);
        let a = template
            .bind_named(
                &BTreeMap::from([("a".into(), 4.0), ("b".into(), 1.0), ("k".into(), 2.0)]),
                None,
            )
            .unwrap();
        let b = template.bind(&[4.0, 1.0, 2.0], None).unwrap();
        assert_eq!(a.model_bytes, b.model_bytes);
    }
    #[test]
    fn bindings_replay_across_threads_and_preserve_cancellation() {
        let model = model(false);
        let one =
            ParameterTemplate::prepare(&model, &request(&model, "direct_cr", "ziggurat", true, 1))
                .unwrap();
        let four =
            ParameterTemplate::prepare(&model, &request(&model, "direct_cr", "ziggurat", true, 4))
                .unwrap();
        let a = one.bind(&[3.0, 0.5, 2.0], None).unwrap();
        let b = four.bind(&[3.0, 0.5, 2.0], None).unwrap();
        let control = simulation::Control::default();
        let mut x = serde_json::to_value(ssa::execute(&a.prepared, &control).unwrap()).unwrap();
        let mut y = serde_json::to_value(ssa::execute(&b.prepared, &control).unwrap()).unwrap();
        x["run_sha256"] = json!(null);
        y["run_sha256"] = json!(null);
        assert_eq!(x, y);
        control.cancel();
        let stopped = ssa::execute(&a.prepared, &control).unwrap();
        assert_eq!(stopped.termination, ssa::Termination::Cancelled);
        assert!(stopped.statistics.is_none());
    }
    #[test]
    fn topology_identity_ignores_parameters_but_keeps_literals_and_populations() {
        let raw = model(false);
        let template =
            ParameterTemplate::prepare(&raw, &request(&raw, "direct", "inverse", false, 1))
                .unwrap();
        let bound = template.bind(&[4.0, 1.0, 2.0], None).unwrap();
        assert_eq!(
            template.topology_sha256(),
            bound.prepared.core.network.topology_sha256()
        );
        let mut model = crate::parse_model(&raw).unwrap();
        model.species[0].initial_count += 1;
        assert_ne!(
            template.topology_sha256(),
            crate::network::topology_sha256(&model)
        );
    }
    #[test]
    fn binding_errors_name_the_mismatch() {
        let raw = model(false);
        let template =
            ParameterTemplate::prepare(&raw, &request(&raw, "direct", "inverse", false, 1))
                .unwrap();
        let count = template.bind(&[1.0, 2.0], None).err().unwrap();
        assert_eq!(
            (count[0].code.as_str(), count[0].path.as_str()),
            ("CRN-E014", "/parameters")
        );
        let named = BTreeMap::from([
            ("a".to_string(), 1.0),
            ("b".to_string(), 1.0),
            ("z".to_string(), 1.0),
        ]);
        let names = template.bind_named(&named, None).err().unwrap();
        assert_eq!(names[0].code, "CRN-E014");
        assert!(
            names[0].message.contains("missing [\"k\"]")
                && names[0].message.contains("unknown [\"z\"]")
        );
        let value = template.bind(&[f64::NAN, 1.0, 1.0], None).err().unwrap();
        assert_eq!(value[0].code, "CRN-E015");
    }
    #[test]
    fn scheduler_state_evaluator_matches_checked_propensity() {
        let raw = model(true);
        let template =
            ParameterTemplate::prepare(&raw, &request(&raw, "direct", "inverse", false, 1))
                .unwrap();
        let network = &template.prepared.core.network;
        for x in [0, 1, 2, 4, 1000] {
            for rx in 0..2 {
                assert_eq!(
                    network
                        .propensity_in_scheduler_state(rx, &[x])
                        .unwrap()
                        .to_bits(),
                    network.propensity(rx, &[x]).unwrap().to_bits()
                );
            }
        }
        // A wrong length still gets the checked path and its diagnostic.
        assert_eq!(
            network
                .propensity_in_scheduler_state(0, &[1, 2])
                .unwrap_err()
                .code,
            network.propensity(0, &[1, 2]).unwrap_err().code
        );
    }
}
