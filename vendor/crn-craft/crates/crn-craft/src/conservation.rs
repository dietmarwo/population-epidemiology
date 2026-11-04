//! Small-network exact left-nullspace analysis. No floating-point rank decisions.
use crate::{diagnostic::Diagnostic, network::CompiledNetwork};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

pub const MAX_CONSERVATION_SPECIES: usize = 64;
pub const MAX_CONSERVATION_REACTIONS: usize = 256;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ConservationStatus {
    Exact,
    SizeLimit,
    ArithmeticOverflow,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ConservationLaw {
    /// Primitive integer coefficients, encoded as strings to preserve i128 values.
    pub coefficients: Vec<String>,
    pub initial_total: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ConservationReport {
    pub status: ConservationStatus,
    pub rank: Option<usize>,
    pub basis: Vec<ConservationLaw>,
    /// Proven using non-negative basis rows only; None is not an unboundedness proof.
    pub initial_class_upper_bounds: Vec<Option<String>>,
    pub diagnostics: Vec<Diagnostic>,
}

fn gcd(mut a: i128, mut b: i128) -> Result<i128, ()> {
    a = a.checked_abs().ok_or(())?;
    b = b.checked_abs().ok_or(())?;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    Ok(a)
}

fn primitive(row: &mut [i128]) -> Result<(), ()> {
    let mut divisor = 0;
    for &value in row.iter() {
        divisor = gcd(divisor, value)?;
    }
    if divisor == 0 {
        return Ok(());
    }
    let negative = row.iter().find(|&&x| x != 0).is_some_and(|&x| x < 0);
    for value in row {
        *value /= divisor;
        if negative {
            *value = value.checked_neg().ok_or(())?;
        }
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct Rational {
    n: i128,
    d: i128,
}
impl Rational {
    const ZERO: Self = Self { n: 0, d: 1 };
    fn new(n: i128, d: i128) -> Result<Self, ()> {
        if d == 0 {
            return Err(());
        }
        let divisor = gcd(n, d)?;
        let (mut n, mut d) = (n / divisor, d / divisor);
        if d < 0 {
            n = n.checked_neg().ok_or(())?;
            d = d.checked_neg().ok_or(())?;
        }
        Ok(Self { n, d })
    }
    fn add(self, other: Self) -> Result<Self, ()> {
        let divisor = gcd(self.d, other.d)?;
        let left = other.d / divisor;
        let right = self.d / divisor;
        let n = self
            .n
            .checked_mul(left)
            .ok_or(())?
            .checked_add(other.n.checked_mul(right).ok_or(())?)
            .ok_or(())?;
        Self::new(n, self.d.checked_mul(left).ok_or(())?)
    }
    fn times(self, value: i128) -> Result<Self, ()> {
        let divisor = gcd(value, self.d)?;
        Self::new(
            self.n.checked_mul(value / divisor).ok_or(())?,
            self.d / divisor,
        )
    }
    fn divided_by(self, value: i128) -> Result<Self, ()> {
        if value == 0 {
            return Err(());
        }
        let divisor = gcd(self.n, value)?;
        Self::new(
            self.n / divisor,
            self.d.checked_mul(value / divisor).ok_or(())?,
        )
    }
}

/// Integer row elimination with gcd reduction; checked rational back substitution.
fn nullspace(mut rows: Vec<Vec<i128>>, width: usize) -> Result<(usize, Vec<Vec<i128>>), ()> {
    let mut pivots = Vec::new();
    for column in 0..width {
        let rank = pivots.len();
        let Some(pivot) = (rank..rows.len()).find(|&i| rows[i][column] != 0) else {
            continue;
        };
        rows.swap(rank, pivot);
        primitive(&mut rows[rank])?;
        for i in rank + 1..rows.len() {
            if rows[i][column] == 0 {
                continue;
            }
            let divisor = gcd(rows[rank][column], rows[i][column])?;
            let a = rows[rank][column] / divisor;
            let b = rows[i][column] / divisor;
            let (previous, remaining) = rows.split_at_mut(i);
            for (value, &pivot_value) in remaining[0].iter_mut().zip(&previous[rank]).skip(column) {
                *value = value
                    .checked_mul(a)
                    .ok_or(())?
                    .checked_sub(pivot_value.checked_mul(b).ok_or(())?)
                    .ok_or(())?;
            }
            primitive(&mut remaining[0])?;
        }
        pivots.push(column);
        if pivots.len() == rows.len() {
            break;
        }
    }
    let mut basis = Vec::new();
    for free in (0..width).filter(|i| !pivots.contains(i)) {
        let mut vector = vec![Rational::ZERO; width];
        vector[free] = Rational { n: 1, d: 1 };
        for (row, &pivot) in pivots.iter().enumerate().rev() {
            let mut sum = Rational::ZERO;
            for (j, &value) in vector.iter().enumerate().skip(pivot + 1) {
                sum = sum.add(value.times(rows[row][j])?)?;
            }
            vector[pivot] = sum.times(-1)?.divided_by(rows[row][pivot])?;
        }
        let mut common = 1_i128;
        for value in &vector {
            common = (common / gcd(common, value.d)?)
                .checked_mul(value.d)
                .ok_or(())?;
        }
        let mut integers: Vec<i128> = vector
            .iter()
            .map(|v| v.n.checked_mul(common / v.d).ok_or(()))
            .collect::<Result<_, _>>()?;
        primitive(&mut integers)?;
        basis.push(integers);
    }
    expose_nonnegative_laws(&mut basis)?;
    Ok((pivots.len(), basis))
}

// Elementary row combinations preserve the basis span/rank. This limited
// heuristic exposes familiar totals (e.g. substrate+complex+product); it does
// not enumerate the non-negative conservation cone.
fn expose_nonnegative_laws(basis: &mut [Vec<i128>]) -> Result<(), ()> {
    for _ in 0..basis.len() {
        let mut changed = false;
        for i in 0..basis.len() {
            if basis[i].iter().all(|&x| x >= 0) {
                continue;
            }
            let mut candidate = basis[i].clone();
            for column in 0..candidate.len() {
                if candidate[column] >= 0 {
                    continue;
                }
                if let Some(other) = basis
                    .iter()
                    .enumerate()
                    .find(|(j, row)| *j != i && row[column] > 0 && row.iter().all(|&x| x >= 0))
                    .map(|(_, row)| row)
                {
                    let needed = candidate[column].checked_neg().ok_or(())?;
                    let multiple = (needed - 1) / other[column] + 1;
                    for (value, &coefficient) in candidate.iter_mut().zip(other) {
                        *value = value
                            .checked_add(multiple.checked_mul(coefficient).ok_or(())?)
                            .ok_or(())?;
                    }
                }
            }
            if candidate.iter().all(|&x| x >= 0) {
                primitive(&mut candidate)?;
                basis[i] = candidate;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    Ok(())
}

pub fn analyze_conservation(network: &CompiledNetwork) -> ConservationReport {
    if network.species().len() > MAX_CONSERVATION_SPECIES
        || network.reactions().len() > MAX_CONSERVATION_REACTIONS
    {
        return unavailable(
            network,
            ConservationStatus::SizeLimit,
            "CRN-W010",
            "exact conservation analysis exceeds its 64-species/256-reaction cap",
        );
    }
    match exact_report(network) {
        Ok(report) => report,
        Err(()) => unavailable(
            network,
            ConservationStatus::ArithmeticOverflow,
            "CRN-W011",
            "checked i128 conservation arithmetic overflowed; no partial basis is reported",
        ),
    }
}

fn unavailable(
    network: &CompiledNetwork,
    status: ConservationStatus,
    code: &str,
    message: &str,
) -> ConservationReport {
    ConservationReport {
        status,
        rank: None,
        basis: Vec::new(),
        initial_class_upper_bounds: vec![None; network.species().len()],
        diagnostics: vec![Diagnostic::warning(
            code,
            "/conservation",
            message,
            "Use a smaller network or an independent arbitrary-precision analyzer; no absence of conservation laws is implied.",
        )],
    }
}

fn exact_report(network: &CompiledNetwork) -> Result<ConservationReport, ()> {
    let n = network.species().len();
    // Rows of N^T, so its nullspace is the species conservation space.
    let rows: Vec<Vec<i128>> = network
        .reactions()
        .iter()
        .map(|reaction| {
            let mut row = vec![0; n];
            for term in &reaction.net_changes {
                row[term.species] = i128::from(term.coefficient);
            }
            row
        })
        .collect();
    let (rank, basis) = nullspace(rows, n)?;
    let mut laws = Vec::with_capacity(basis.len());
    let mut bounds = vec![None::<i128>; n];
    for vector in basis {
        // Verify every returned law against the original sparse columns.
        for reaction in network.reactions() {
            let mut residual = 0_i128;
            for term in &reaction.net_changes {
                residual = residual
                    .checked_add(
                        vector[term.species]
                            .checked_mul(i128::from(term.coefficient))
                            .ok_or(())?,
                    )
                    .ok_or(())?;
            }
            if residual != 0 {
                return Err(());
            }
        }
        let mut total = 0_i128;
        for (&weight, &count) in vector.iter().zip(network.initial_counts()) {
            total = total
                .checked_add(weight.checked_mul(i128::from(count)).ok_or(())?)
                .ok_or(())?;
        }
        if vector.iter().all(|&x| x >= 0) {
            for (i, &weight) in vector.iter().enumerate() {
                if weight > 0 {
                    let bound = total / weight;
                    bounds[i] = Some(bounds[i].map_or(bound, |old| old.min(bound)));
                }
            }
        }
        laws.push(ConservationLaw {
            coefficients: vector.iter().map(ToString::to_string).collect(),
            initial_total: total.to_string(),
        });
    }
    Ok(ConservationReport {
        status: ConservationStatus::Exact,
        rank: Some(rank),
        basis: laws,
        initial_class_upper_bounds: bounds
            .into_iter()
            .map(|x| x.map(|v| v.to_string()))
            .collect(),
        diagnostics: Vec::new(),
    })
}
