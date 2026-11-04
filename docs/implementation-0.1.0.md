# population-epidemiology implementation log

3 November 2026. Version 0.1.0. This delivery uses the bounded research
scope recorded below. The independent references and optional comparisons
are identified separately for each increment.

## Increment results

| Increment | Status | Implementation and evidence |
| --- | --- | --- |
| 0 | Public fallback | The original Numba research model was unavailable. Audited CRN-craft FSP and checked SSA and chose transparent fixed-population SIR and neutral two-mutant Moran fixtures. |
| 1 | Implemented | SIR model construction uses supported CRN count mass-action laws. The native SIR path calls existing CRN-craft FSP. Moran uses a separate tiny discrete recurrence. |
| 2 | Implemented | Independent SciPy matrix exponential and an exhaustive labeled-population oracle verify generator convention, distributions, absorbing behavior and fixation. Checked SSA agrees within a predeclared sampling bound. |
| 3 | Bounded fallback | Returns complete two-mutant marginals for the declared neutral model and optionally its joint law. SIR retains CRN FSP provenance and numerical diagnostics. This is a separate workload from the unavailable original. |
| 4 | Implemented with reuse limit | Prepared companion models and bounded parallel candidate batches replay across worker counts. The existing FSP API still constructs its projection for each rate binding. |
| 5 | Implemented on bounded subset | Exact piecewise-rate two-state calendar and analytical boundary check. The SIR workflow retains constant rates. |
| 6 | Implemented on synthetic inputs | Twelve grid points fit exact synthetic infected-count means, with held-out validation and tighter FSP reevaluation. No empirical inference or intervention recommendation is claimed. |
| 7 | Deferred by plan gate | A structured extension requires a concrete researcher question. The existing finite-law study supplies the bounded release. |
| 8 | Local fallback | SQLite manifest, replay and cancellation/resume evidence. External campaign contracts and independent researcher review remain open. |
| 9 | Source and wheel release prepared | The sdist includes a pinned MIT CRN core, and a clean wheel built from it. A standalone registry crate failed because CRN-craft is a local source dependency. The ZIP still builds the Rust companion offline when its locked dependencies are cached. |

## Numerical results

The maximum FSP state-probability difference versus dense propagation was **8.91e-13**. The neutral fixation marginal converged to 0.333333333. The checked SSA infected-mean deviation was 0.0358 against a predeclared 0.4024 bound. The held-out SIR fitting MSE was 4.40e-24.

The optional in-process SciPy differential-evolution demonstration made
48 evaluations and replayed its selected candidate.
It used the packaged example and its own local bounds. The domain-specific
grid and holdout evidence are in `evidence/validation.json`.

## Installed performance and replay

The following are medians of three 96-candidate batch observations,
measured in a clean CPython 3.12 environment with the installed wheel.
The scenarios and candidate inputs are identical within a row pair.
Times include native crossing and compact JSON result conversion.

| Boundary | Time |
| --- | ---: |
| Cold preparation and batch, one worker | 17.850 ms |
| Prepared batch, one worker | 17.826 ms |
| Prepared batch, four workers | 5.968 ms |

The replay test reversed candidate order and changed worker count. It
recovered 12 completed SQLite records across cancellation and retained
21 attempts. This campaign uses a local database and package-build
identity. Cross-machine determinism and distributed worker behavior
were outside the test scope.

## Release evidence

- `cargo fmt`, Clippy with warnings denied and locked Rust release tests passed.
- The native sweep and numerical validation ran against a release build.
- A wheel built from the extracted sdist and installed offline into a fresh
  environment. The wheel target is Linux x86-64, CPython 3.10+ ABI. Actual
  verification used CPython 3.12 on x86_64.
- The standalone registry crate check failed at the vendored CRN path dependency. Use the bundled Rust source, Python sdist or wheel.
- Model and source manifests are in `docs/`. Raw timing, numerical and
  clean-install results are in `evidence/`.

The source was staged in successive local commits. Each increment above
records its own outcome, while related implementation changes share
commits.

- `0946062` Define versioned evaluation and packaging boundaries
- `9d0e091` Implement bounded domain evaluator and resumable local execution
- `20003f1` Validate model laws and add bounded search and installation examples
- `724eee4` Prepare source distributions and document validated research scope
- `6451511` Canonicalize campaign replay output
