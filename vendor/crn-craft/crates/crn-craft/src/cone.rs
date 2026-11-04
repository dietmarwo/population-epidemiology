//! Exact extreme rays of {l >= 0 : N^T l = 0}, with bounded double description.
//! Algorithm/adjacency reference: Burton, arXiv:0808.4050, Sections 3.1 and 4.1.
//! Independent implementation; no third-party enumeration code is incorporated.
use crate::{conservation::ConservationLaw, diagnostic::Diagnostic, network::CompiledNetwork};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{cmp::Ordering, collections::BTreeMap};

pub const MAX_CONE_SPECIES: usize = 64;
pub const MAX_CONE_REACTIONS: usize = 256;
pub const MAX_CONE_RAYS: usize = 4096;
pub const MAX_CONE_WORK: u64 = 2_000_000;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ConeLimits {
    pub max_rays: usize,
    pub max_work: u64,
}
impl Default for ConeLimits {
    fn default() -> Self {
        Self {
            max_rays: MAX_CONE_RAYS,
            max_work: MAX_CONE_WORK,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ConeStatus {
    Exact,
    SizeLimit,
    RayLimit,
    WorkLimit,
    ArithmeticOverflow,
    VerificationFailed,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ConeReport {
    pub status: ConeStatus,
    pub limits: ConeLimits,
    pub work_units: u64,
    pub extreme_rays: Vec<ConservationLaw>,
    pub initial_class_upper_bounds: Vec<Option<String>>,
    pub diagnostics: Vec<Diagnostic>,
}
struct Budget {
    used: u64,
    limit: u64,
}
impl Budget {
    fn charge(&mut self) -> Result<(), ConeStatus> {
        if self.used == self.limit {
            return Err(ConeStatus::WorkLimit);
        }
        self.used += 1;
        Ok(())
    }
}
#[derive(Clone)]
struct Ray {
    coordinates: Vec<i128>,
    support: u64,
}
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
fn insert_ray(next: &mut BTreeMap<Vec<i128>, u64>, ray: Ray, max: usize) -> Result<(), ConeStatus> {
    if next.len() == max && !next.contains_key(&ray.coordinates) {
        return Err(ConeStatus::RayLimit);
    }
    next.insert(ray.coordinates, ray.support);
    Ok(())
}

pub fn analyze_cone(network: &CompiledNetwork) -> ConeReport {
    analyze_cone_with_limits(network, ConeLimits::default())
}
/// Callers may lower the deterministic caps; values above hard caps are clamped.
pub fn analyze_cone_with_limits(network: &CompiledNetwork, limits: ConeLimits) -> ConeReport {
    let limits = ConeLimits {
        max_rays: limits.max_rays.min(MAX_CONE_RAYS),
        max_work: limits.max_work.min(MAX_CONE_WORK),
    };
    let mut budget = Budget {
        used: 0,
        limit: limits.max_work,
    };
    let mut report = ConeReport {
        status: ConeStatus::Exact,
        limits,
        work_units: 0,
        extreme_rays: Vec::new(),
        initial_class_upper_bounds: vec![None; network.species().len()],
        diagnostics: Vec::new(),
    };
    match enumerate(network, limits, &mut budget)
        .and_then(|rays| summarize(network, rays, &mut budget))
    {
        Ok((rays, bounds)) => {
            report.extreme_rays = rays;
            report.initial_class_upper_bounds = bounds;
        }
        Err(status) => {
            report.status = status;
            let (code, message) = match status {
                ConeStatus::SizeLimit => (
                    "CRN-W012",
                    "cone exceeds the 64-species/256-reaction size cap",
                ),
                ConeStatus::RayLimit => (
                    "CRN-W013",
                    "cone enumeration exhausted its intermediate-ray cap",
                ),
                ConeStatus::WorkLimit => (
                    "CRN-W013",
                    "cone enumeration exhausted its deterministic work cap",
                ),
                ConeStatus::ArithmeticOverflow => {
                    ("CRN-W014", "checked i128 cone arithmetic overflowed")
                }
                _ => ("CRN-W015", "internal exact cone verification failed"),
            };
            report.diagnostics.push(Diagnostic::warning(code, "/nonnegative_conservation", message,
                "No partial cone is returned. Use an independent analyzer; absence of laws or bounds is not implied."));
        }
    }
    report.work_units = budget.used;
    report
}

fn enumerate(
    network: &CompiledNetwork,
    limits: ConeLimits,
    budget: &mut Budget,
) -> Result<Vec<Ray>, ConeStatus> {
    let n = network.species().len();
    if n > MAX_CONE_SPECIES || network.reactions().len() > MAX_CONE_REACTIONS {
        return Err(ConeStatus::SizeLimit);
    }
    if n > limits.max_rays {
        return Err(ConeStatus::RayLimit);
    }
    let mut rays: Vec<_> = (0..n)
        .map(|i| {
            let mut coordinates = vec![0; n];
            coordinates[i] = 1;
            Ray {
                coordinates,
                support: 1_u64 << i,
            }
        })
        .collect();
    for reaction in network.reactions() {
        let mut dots = Vec::with_capacity(rays.len());
        let mut positive = Vec::new();
        let mut negative = Vec::new();
        let mut next = BTreeMap::new();
        for (i, ray) in rays.iter().enumerate() {
            budget.charge()?;
            let mut value = 0_i128;
            for term in &reaction.net_changes {
                budget.charge()?;
                value = value
                    .checked_add(
                        ray.coordinates[term.species]
                            .checked_mul(i128::from(term.coefficient))
                            .ok_or(ConeStatus::ArithmeticOverflow)?,
                    )
                    .ok_or(ConeStatus::ArithmeticOverflow)?;
            }
            dots.push(value);
            match value.cmp(&0) {
                Ordering::Greater => positive.push(i),
                Ordering::Less => negative.push(i),
                Ordering::Equal => insert_ray(&mut next, ray.clone(), limits.max_rays)?,
            }
        }
        for &p in &positive {
            for &q in &negative {
                budget.charge()?;
                let support = rays[p].support | rays[q].support;
                // In an orthant/equality cone, a pair is adjacent iff no third
                // extreme ray has support contained in their union.
                let mut adjacent = true;
                for (k, other) in rays.iter().enumerate() {
                    budget.charge()?;
                    if k != p && k != q && other.support & !support == 0 {
                        adjacent = false;
                        break;
                    }
                }
                if !adjacent {
                    continue;
                }
                let minus_q = dots[q]
                    .checked_neg()
                    .ok_or(ConeStatus::ArithmeticOverflow)?;
                let divisor = gcd(dots[p], minus_q);
                let a = minus_q / divisor;
                let b = dots[p] / divisor;
                let mut coordinates = Vec::with_capacity(n);
                let mut common = 0;
                for (&x, &y) in rays[p].coordinates.iter().zip(&rays[q].coordinates) {
                    budget.charge()?;
                    let value = x
                        .checked_mul(a)
                        .ok_or(ConeStatus::ArithmeticOverflow)?
                        .checked_add(y.checked_mul(b).ok_or(ConeStatus::ArithmeticOverflow)?)
                        .ok_or(ConeStatus::ArithmeticOverflow)?;
                    common = gcd(common, value);
                    coordinates.push(value);
                }
                for value in &mut coordinates {
                    *value /= common;
                }
                insert_ray(
                    &mut next,
                    Ray {
                        coordinates,
                        support,
                    },
                    limits.max_rays,
                )?;
            }
        }
        rays = next
            .into_iter()
            .map(|(coordinates, support)| Ray {
                coordinates,
                support,
            })
            .collect();
        if rays.is_empty() {
            break;
        }
    }
    rays.sort_by(|a, b| a.coordinates.cmp(&b.coordinates));
    Ok(rays)
}

type LawsAndBounds = (Vec<ConservationLaw>, Vec<Option<String>>);
fn summarize(
    network: &CompiledNetwork,
    rays: Vec<Ray>,
    budget: &mut Budget,
) -> Result<LawsAndBounds, ConeStatus> {
    let mut laws = Vec::new();
    let mut bounds = vec![None::<i128>; network.species().len()];
    for ray in rays {
        for reaction in network.reactions() {
            let mut residual = 0_i128;
            for term in &reaction.net_changes {
                budget.charge()?;
                residual = residual
                    .checked_add(
                        ray.coordinates[term.species]
                            .checked_mul(i128::from(term.coefficient))
                            .ok_or(ConeStatus::ArithmeticOverflow)?,
                    )
                    .ok_or(ConeStatus::ArithmeticOverflow)?;
            }
            if residual != 0 {
                return Err(ConeStatus::VerificationFailed);
            }
        }
        let mut total = 0_i128;
        for (&weight, &count) in ray.coordinates.iter().zip(network.initial_counts()) {
            budget.charge()?;
            total = total
                .checked_add(
                    weight
                        .checked_mul(i128::from(count))
                        .ok_or(ConeStatus::ArithmeticOverflow)?,
                )
                .ok_or(ConeStatus::ArithmeticOverflow)?;
        }
        for (i, &weight) in ray.coordinates.iter().enumerate() {
            if weight > 0 {
                let bound = total / weight;
                bounds[i] = Some(bounds[i].map_or(bound, |old| old.min(bound)));
            }
        }
        laws.push(ConservationLaw {
            coefficients: ray.coordinates.iter().map(ToString::to_string).collect(),
            initial_total: total.to_string(),
        });
    }
    Ok((
        laws,
        bounds
            .into_iter()
            .map(|v| v.map(|x| x.to_string()))
            .collect(),
    ))
}
