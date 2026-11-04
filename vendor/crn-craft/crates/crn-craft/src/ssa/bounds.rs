//! Integer-box bounds for the actual gated evaluator, not idealized rate formulas.
//! RSSA family: Thanh et al., J. Chemical Physics 141, 134116 (2014).
use super::{
    Bail, Draws,
    selectors::{Composition, RateTree},
};
use crate::{network::CompiledNetwork, rates::CompiledRate};

#[derive(Clone, Copy)]
struct CachedRange {
    start: i64,
    end: i64,
    minimum: f64,
    maximum: f64,
    positive: f64,
}
#[derive(Default)]
pub(super) struct Bounds {
    low: Vec<i64>,
    high: Vec<i64>,
    scratch: Vec<i64>,
    pinned: Vec<bool>,
    pub(super) lower: Vec<f64>,
    pub(super) positive_lower: usize,
    marks: Vec<u32>,
    generation: u32,
    dirty: Vec<usize>,
    ranges: Vec<Option<CachedRange>>,
}

fn interval(count: i64, pinned: bool) -> (i64, i64) {
    if pinned || count == 0 {
        return (count, count);
    }
    let width = (count / 10).clamp(4, 64);
    (
        count.saturating_sub(width).max(0),
        count.saturating_add(width),
    )
}

fn elementary(rate: &CompiledRate, rx: &crate::network::ReactionStructure) -> bool {
    match rate {
        CompiledRate::MassAction(c, _) => {
            if *c == 0.0 || rx.reactants.is_empty() {
                return true;
            }
            let order: i64 = rx.reactants.iter().map(|t| t.coefficient).sum();
            // Positive normal coefficient and at most two factors >= 1. Every
            // active product stays normal and finite, hence takes the monotone
            // fast branch in rates.rs (including the n=2/3 homodimer cases).
            order <= 2 && c.is_normal() && *c <= f64::MAX / (i64::MAX as f64).powi(2)
        }
        _ => true,
    }
}

impl Bounds {
    pub(super) fn reset(&mut self, network: &CompiledNetwork, state: &[i64]) {
        self.low.resize(state.len(), 0);
        self.high.resize(state.len(), 0);
        self.pinned.resize(state.len(), false);
        self.pinned.fill(false);
        for (i, rx) in network.reactions().iter().enumerate() {
            if !elementary(network.compiled_rate(i), rx) {
                for &s in &rx.reads {
                    self.pinned[s] = true;
                }
            }
        }
        for (s, &x) in state.iter().enumerate() {
            (self.low[s], self.high[s]) = interval(x, self.pinned[s]);
        }
        self.scratch.clone_from(&self.high);
        self.lower.resize(network.reactions().len(), 0.0);
        self.lower.fill(0.0);
        self.positive_lower = 0;
        self.marks.resize(network.reactions().len(), 0);
        self.marks.fill(0);
        self.generation = 0;
        self.dirty.clear();
        self.ranges.fill(None);
    }

    fn evaluate(
        network: &CompiledNetwork,
        rx: usize,
        state: &[i64],
        draws: &mut Draws,
    ) -> Result<f64, Bail> {
        draws.poll()?;
        draws.work.propensity_evaluations += 1;
        Ok(network.propensity_in_valid_state(rx, state)?)
    }

    fn bounded(
        &mut self,
        network: &CompiledNetwork,
        rx: usize,
        group: usize,
        draws: &mut Draws,
    ) -> Result<(f64, f64, f64), Bail> {
        let structure = &network.reactions()[rx];
        if structure
            .reactants
            .iter()
            .any(|t| self.high[t.species] < t.coefficient)
        {
            return Ok((0.0, 0.0, f64::INFINITY));
        }
        let gated = structure
            .reactants
            .iter()
            .any(|t| self.low[t.species] < t.coefficient);
        let rate = network.compiled_rate(rx);
        if let Some(s) = rate.extra_read() {
            // Other coordinates influence only availability. Enumerate all
            // feasible integer regulator values; this also handles a consumed
            // repressor whose gated propensity is nonmonotone over the box.
            let gate = structure
                .reactants
                .iter()
                .find(|t| t.species == s)
                .map_or(0, |t| t.coefficient);
            let start = self.low[s].max(gate);
            if group != usize::MAX {
                if let Some(Some(cached)) = self.ranges.get(group) {
                    if cached.start == start && cached.end == self.high[s] {
                        draws.poll()?;
                        return Ok((
                            if gated { 0.0 } else { cached.minimum },
                            cached.maximum,
                            cached.positive,
                        ));
                    }
                }
            }
            let mut minimum = f64::INFINITY;
            let mut positive = f64::INFINITY;
            let mut maximum: f64 = 0.0;
            for value in start..=self.high[s] {
                self.scratch[s] = value;
                let result = Self::evaluate(network, rx, &self.scratch, draws);
                self.scratch[s] = self.high[s];
                let a = result?;
                minimum = minimum.min(a);
                maximum = maximum.max(a);
                if a > 0.0 {
                    positive = positive.min(a);
                }
            }
            if group != usize::MAX {
                if self.ranges.len() <= group {
                    self.ranges.resize(group + 1, None);
                }
                self.ranges[group] = Some(CachedRange {
                    start,
                    end: self.high[s],
                    minimum,
                    maximum,
                    positive,
                });
            }
            Ok((if gated { 0.0 } else { minimum }, maximum, positive))
        } else {
            let high = Self::evaluate(network, rx, &self.high, draws)?;
            let low = if gated {
                0.0
            } else {
                Self::evaluate(network, rx, &self.low, draws)?
            };
            if !gated {
                return Ok((low, high, if low > 0.0 { low } else { f64::INFINITY }));
            }
            // Lowest *enabled* state supplies the minimum positive mass-action
            // rate even when a lower corner is disabled by an availability gate.
            for t in &structure.reactants {
                self.scratch[t.species] = self.low[t.species].max(t.coefficient);
            }
            let result = Self::evaluate(network, rx, &self.scratch, draws);
            for t in &structure.reactants {
                self.scratch[t.species] = self.high[t.species];
            }
            let a = result?;
            Ok((low, high, if a > 0.0 { a } else { f64::INFINITY }))
        }
    }

    #[cfg(test)]
    fn bound(
        &mut self,
        network: &CompiledNetwork,
        rx: usize,
        state: &[i64],
        draws: &mut Draws,
    ) -> Result<(f64, f64, f64), Bail> {
        self.bound_group(network, rx, state, usize::MAX, draws)
    }
    fn bound_group(
        &mut self,
        network: &CompiledNetwork,
        rx: usize,
        state: &[i64],
        group: usize,
        draws: &mut Draws,
    ) -> Result<(f64, f64, f64), Bail> {
        if network.reactions()[rx]
            .reads
            .iter()
            .all(|&s| self.low[s] == self.high[s])
        {
            let a = Self::evaluate(network, rx, state, draws)?;
            if !elementary(network.compiled_rate(rx), &network.reactions()[rx]) {
                draws.work.point_fallbacks += 1;
            }
            return Ok((a, a, if a > 0.0 { a } else { f64::INFINITY }));
        }
        match self.bounded(network, rx, group, draws) {
            Ok(value) => Ok(value),
            Err(Bail::Numerical(_)) => {
                // A neighboring state may be unrepresentable even if the
                // current state is valid. Shrink every read coordinate to a
                // point and keep it pinned for the remainder of this path.
                let a = Self::evaluate(network, rx, state, draws)?;
                for &s in &network.reactions()[rx].reads {
                    self.pinned[s] = true;
                    self.low[s] = state[s];
                    self.high[s] = state[s];
                    self.scratch[s] = state[s];
                }
                draws.work.point_fallbacks += 1;
                Ok((a, a, if a > 0.0 { a } else { f64::INFINITY }))
            }
            Err(other) => Err(other),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn refresh(
        &mut self,
        network: &CompiledNetwork,
        readers: &[Vec<usize>],
        groups: &[usize],
        state: &[i64],
        previous: Option<usize>,
        tree: &mut RateTree,
        composition: &mut Composition,
        draws: &mut Draws,
    ) -> Result<(), Bail> {
        self.dirty.clear();
        if let Some(fired) = previous {
            self.generation = self.generation.wrapping_add(1);
            if self.generation == 0 {
                self.marks.fill(0);
                self.generation = 1;
            }
            for term in &network.reactions()[fired].net_changes {
                let s = term.species;
                if state[s] >= self.low[s] && state[s] <= self.high[s] {
                    continue;
                }
                (self.low[s], self.high[s]) = interval(state[s], self.pinned[s]);
                self.scratch[s] = self.high[s];
                for (k, &rx) in readers[s].iter().enumerate() {
                    if k % 64 == 0 {
                        draws.poll()?;
                    }
                    if self.marks[rx] != self.generation {
                        self.marks[rx] = self.generation;
                        self.dirty.push(rx);
                    }
                }
            }
            // Fixed index order makes fallback shrinking deterministic.
            self.dirty.sort_unstable();
        } else {
            self.dirty.extend(0..network.reactions().len());
        }
        for k in 0..self.dirty.len() {
            draws.poll()?;
            let rx = self.dirty[k];
            let (low, high, minimum) = self.bound_group(network, rx, state, groups[rx], draws)?;
            self.positive_lower -= usize::from(self.lower[rx] > 0.0);
            self.positive_lower += usize::from(low > 0.0);
            self.lower[rx] = low;
            tree.set(rx, high, minimum);
            composition.update(rx, high)?;
            draws.work.bound_updates += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        simulation::{Control, Poller},
        ssa::Work,
    };
    use serde_json::json;
    use std::time::Instant;

    fn network(rate: serde_json::Value, reactants: serde_json::Value) -> CompiledNetwork {
        let value = json!({"schema":"crn-craft.model/0.1","name":"Bounds","convention":"stochastic_counts","omega":1.,"parameters":{},
            "species":[{"name":"X","initial_count":3},{"name":"G","initial_count":3}],
            "reactions":[{"name":"r","reactants":reactants,"products":{},"rate":rate}]});
        CompiledNetwork::compile(&crate::parse_model(&serde_json::to_vec(&value).unwrap()).unwrap())
            .unwrap()
    }

    #[test]
    fn shared_ranges_preserve_distinct_gates_and_feasible_regulator_intervals() {
        let value = json!({"schema":"crn-craft.model/0.1","name":"Shared","convention":"stochastic_counts","omega":1.,"parameters":{},
          "species":[{"name":"X","initial_count":3},{"name":"G","initial_count":3}],
          "reactions":[
          {"name":"first","reactants":{"G":1},"products":{},"rate":{"kind":"michaelis_menten","substrate":"X","vmax":3.,"km":2.}},
          {"name":"second","reactants":{"G":2,"X":1},"products":{},"rate":{"kind":"michaelis_menten","substrate":"X","vmax":3.,"km":2.}}]});
        let net = CompiledNetwork::compile(
            &crate::parse_model(&serde_json::to_vec(&value).unwrap()).unwrap(),
        )
        .unwrap();
        assert_eq!(
            net.compiled_rate(0).bounds_key(),
            net.compiled_rate(1).bounds_key()
        );
        let control = Control::default();
        let mut poll = Poller::new(&control, Instant::now(), None);
        let mut work = Work::default();
        let mut draws = Draws {
            limit: 0,
            work: &mut work,
            poller: &mut poll,
            fast: false,
        };
        let mut b = Bounds::default();
        b.reset(&net, &[3, 3]);
        b.low = vec![2, 2];
        b.high = vec![6, 6];
        b.scratch.clone_from(&b.high);
        let first = b.bound_group(&net, 0, &[3, 3], 0, &mut draws).ok().unwrap();
        let count = draws.work.propensity_evaluations;
        assert_eq!(
            first,
            b.bound_group(&net, 1, &[3, 3], 0, &mut draws).ok().unwrap()
        );
        assert_eq!(count, draws.work.propensity_evaluations);
        b.low = vec![0, 0];
        for rx in 0..2 {
            let (low, high, _) = b
                .bound_group(&net, rx, &[3, 3], 0, &mut draws)
                .ok()
                .unwrap();
            for x in 0..=6 {
                for g in 0..=6 {
                    let a = net.propensity(rx, &[x, g]).unwrap();
                    assert!(a >= low && a <= high);
                }
            }
        }
        assert!(draws.work.propensity_evaluations > count);
        b.reset(&net, &[3, 3]);
        assert!(b.ranges.iter().all(Option::is_none));
    }
    #[test]
    fn enabled_mass_action_corner_is_reused_as_the_positive_minimum() {
        let net = network(
            json!({"kind":"mass_action","coefficient":0.1}),
            json!({"X":2}),
        );
        let control = Control::default();
        let mut poll = Poller::new(&control, Instant::now(), None);
        let mut work = Work::default();
        let mut draws = Draws {
            limit: 0,
            work: &mut work,
            poller: &mut poll,
            fast: false,
        };
        let mut b = Bounds::default();
        b.reset(&net, &[4, 3]);
        b.low = vec![4, 0];
        b.high = vec![8, 6];
        b.scratch.clone_from(&b.high);
        let (low, high, minimum) = b.bound(&net, 0, &[4, 3], &mut draws).ok().unwrap();
        assert_eq!(low, minimum);
        assert_eq!(low, net.propensity(0, &[4, 0]).unwrap());
        assert_eq!(high, net.propensity(0, &[8, 6]).unwrap());
        assert_eq!(draws.work.propensity_evaluations, 2);
    }
    #[test]
    fn every_integer_state_is_enclosed_in_gated_custom_and_mass_action_boxes() {
        let mut cases = vec![
            (
                json!({"kind":"mass_action","coefficient":0.1}),
                json!({"X":2}),
            ),
            (
                json!({"kind":"mass_action","coefficient":0.03}),
                json!({"X":1,"G":1}),
            ),
            (
                json!({"kind":"michaelis_menten","substrate":"X","vmax":3.,"km":2.}),
                json!({"G":2}),
            ),
            (
                json!({"kind":"michaelis_menten","substrate":"X","vmax":3.,"km":2.}),
                json!({"X":1,"G":2}),
            ),
        ];
        for kind in ["hill_activation", "hill_repression"] {
            for gate in [json!({"G":2}), json!({"X":2,"G":1})] {
                cases.push((
                    json!({"kind":kind,"regulator":"X","vmax":3.,"basal":0.,"k":2.,"n":2.5}),
                    gate,
                ));
            }
        }
        for (rate, gate) in cases {
            let net = network(rate, gate);
            let control = Control::default();
            let mut poll = Poller::new(&control, Instant::now(), None);
            let mut work = Work::default();
            let mut draws = Draws {
                limit: 0,
                work: &mut work,
                poller: &mut poll,
                fast: false,
            };
            let mut b = Bounds::default();
            b.reset(&net, &[3, 3]);
            for start in 0..=3 {
                b.low = vec![start, 0];
                b.high = vec![6, 6];
                b.scratch.clone_from(&b.high);
                let (low, high, minimum) = b.bound(&net, 0, &[3, 3], &mut draws).ok().unwrap();
                let mut actual_min = f64::INFINITY;
                for x in start..=6 {
                    for g in 0..=6 {
                        // Independent enumeration through the public checked API.
                        let a = net.propensity(0, &[x, g]).unwrap();
                        assert!(a >= low && a <= high, "{x}/{g}: {low} <= {a} <= {high}");
                        if a > 0. {
                            actual_min = actual_min.min(a);
                        }
                    }
                }
                assert_eq!(minimum, actual_min);
            }
        }
    }

    #[test]
    fn unsupported_mass_action_uses_exact_points_and_refreshes_on_change() {
        let net = network(
            json!({"kind":"mass_action","coefficient":0.1}),
            json!({"X":3}),
        );
        let mut b = Bounds::default();
        b.reset(&net, &[3, 3]);
        assert_eq!(b.low[0], 3);
        assert_eq!(b.high[0], 3);
        assert!(b.pinned[0]);
        let control = Control::default();
        let mut poll = Poller::new(&control, Instant::now(), None);
        let mut work = Work::default();
        let mut draws = Draws {
            limit: 0,
            work: &mut work,
            poller: &mut poll,
            fast: false,
        };
        let mut tree = RateTree::default();
        tree.reset(1);
        let mut comp = Composition::default();
        comp.reset(1);
        let readers = vec![vec![0], vec![]];
        b.refresh(
            &net,
            &readers,
            &[usize::MAX],
            &[3, 3],
            None,
            &mut tree,
            &mut comp,
            &mut draws,
        )
        .ok()
        .unwrap();
        assert_eq!(tree.total().unwrap(), 0.1);
        b.refresh(
            &net,
            &readers,
            &[usize::MAX],
            &[0, 3],
            Some(0),
            &mut tree,
            &mut comp,
            &mut draws,
        )
        .ok()
        .unwrap();
        assert_eq!(tree.total().unwrap(), 0.0);
        assert_eq!(draws.work.point_fallbacks, 2);
    }

    #[test]
    fn neighboring_underflow_pins_reads_without_masking_current_state_errors() {
        let net = network(
            json!({"kind":"hill_repression","regulator":"X","vmax":1.,"basal":0.,"k":1.,"n":1000.}),
            json!({}),
        );
        let control = Control::default();
        let mut poll = Poller::new(&control, Instant::now(), None);
        let mut work = Work::default();
        let mut draws = Draws {
            limit: 0,
            work: &mut work,
            poller: &mut poll,
            fast: false,
        };
        let mut b = Bounds::default();
        b.reset(&net, &[1, 3]);
        let (low, high, _) = b.bound(&net, 0, &[1, 3], &mut draws).ok().unwrap();
        assert_eq!((low, high), (0.5, 0.5));
        assert!(b.pinned[0]);
        assert_eq!(draws.work.point_fallbacks, 1);
        b.reset(&net, &[5, 3]);
        assert!(matches!(
            b.bound(&net, 0, &[5, 3], &mut draws),
            Err(Bail::Numerical(_))
        ));
    }
}
