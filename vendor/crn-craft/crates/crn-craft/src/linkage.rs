//! Reaction complexes and graph components. No kinetic theorem is inferred.
use crate::network::{CompiledNetwork, Term};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ComplexEdge {
    pub source: usize,
    pub target: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LinkageReport {
    /// First occurrence order while visiting each reaction's reactants/products.
    pub complexes: Vec<Vec<Term>>,
    pub reaction_edges: Vec<ComplexEdge>,
    pub linkage_classes: Vec<Vec<usize>>,
    pub strong_linkage_classes: Vec<Vec<usize>>,
    pub terminal_strong_linkage_classes: Vec<Vec<usize>>,
    pub weakly_reversible: bool,
    /// Absent when the exact stoichiometric rank is unavailable.
    pub deficiency: Option<usize>,
}

pub fn analyze_linkage(network: &CompiledNetwork, exact_rank: Option<usize>) -> LinkageReport {
    let mut complexes = Vec::new();
    let mut index = BTreeMap::new();
    let mut identify = |side: &[Term]| {
        let key: Vec<_> = side.iter().map(|t| (t.species, t.coefficient)).collect();
        *index.entry(key).or_insert_with(|| {
            let id = complexes.len();
            complexes.push(side.to_vec());
            id
        })
    };
    let edges: Vec<_> = network
        .reactions()
        .iter()
        .map(|r| ComplexEdge {
            source: identify(&r.reactants),
            target: identify(&r.products),
        })
        .collect();
    let n = complexes.len();
    let mut out = vec![Vec::new(); n];
    let mut reverse = vec![Vec::new(); n];
    let mut undirected = vec![Vec::new(); n];
    for edge in &edges {
        out[edge.source].push(edge.target);
        reverse[edge.target].push(edge.source);
        undirected[edge.source].push(edge.target);
        undirected[edge.target].push(edge.source);
    }
    let linkage_classes = components(&undirected, 0..n);
    // Iterative Kosaraju: no call-stack dependence on network depth.
    let mut visited = vec![false; n];
    let mut finish = Vec::with_capacity(n);
    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut stack = vec![(start, 0)];
        while let Some((node, next)) = stack.last_mut() {
            if *next == out[*node].len() {
                finish.push(*node);
                stack.pop();
            } else {
                let neighbor = out[*node][*next];
                *next += 1;
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push((neighbor, 0));
                }
            }
        }
    }
    let strong = components(&reverse, finish.into_iter().rev());
    let mut membership = vec![0; n];
    for (id, class) in strong.iter().enumerate() {
        for &node in class {
            membership[node] = id;
        }
    }
    let mut terminal = vec![true; strong.len()];
    for edge in &edges {
        if membership[edge.source] != membership[edge.target] {
            terminal[membership[edge.source]] = false;
        }
    }
    let terminal_strong_linkage_classes = strong
        .iter()
        .zip(terminal)
        .filter(|(_, t)| *t)
        .map(|(c, _)| c.clone())
        .collect();
    let weakly_reversible = linkage_classes == strong;
    let deficiency =
        exact_rank.and_then(|rank| n.checked_sub(linkage_classes.len())?.checked_sub(rank));
    LinkageReport {
        complexes,
        reaction_edges: edges,
        linkage_classes,
        strong_linkage_classes: strong,
        terminal_strong_linkage_classes,
        weakly_reversible,
        deficiency,
    }
}
fn components(graph: &[Vec<usize>], order: impl Iterator<Item = usize>) -> Vec<Vec<usize>> {
    let mut visited = vec![false; graph.len()];
    let mut classes = Vec::new();
    for start in order {
        if visited[start] {
            continue;
        }
        let mut class = Vec::new();
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(node) = stack.pop() {
            class.push(node);
            for &next in &graph[node] {
                if !visited[next] {
                    visited[next] = true;
                    stack.push(next);
                }
            }
        }
        class.sort_unstable();
        classes.push(class);
    }
    classes.sort();
    classes
}
