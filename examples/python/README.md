# Python examples

Each script runs against an installed `population_epidemiology` 0.3.0
wheel. The optimization scripts also need NumPy, SciPy and fcmaes 2.1.0.
Run a script from any directory, for example
`python examples/python/sir_rate_landscape.py`. The output of the run
recorded for the 0.3.0 release is saved under [`outputs/`](outputs/).
Timings in the outputs come from a shared 2-vCPU host and vary.

## Application examples

| Script | Question | Model and APIs |
| --- | --- | --- |
| [`sir_rate_landscape.py`](sir_rate_landscape.py) | Which infection and recovery rates explain the packaged outbreak, and how sharply are they identified? | SIR; one `score_matrix` call over a 41 x 41 grid, full replay with `evaluate` |
| [`sir_design_comparison.py`](sir_design_comparison.py) | How do early, spread and late sampling designs constrain beta, gamma and R0? | SIR with synthetic observations; `score_matrix` perturbations, `batch` finite-difference sensitivities, tolerance diagnostics |
| [`moran_fixation_study.py`](moran_fixation_study.py) | How quickly does neutral drift fix one of two mutants, and does it reach the initial-frequency limit? | Moran populations 4 to 32; traced joint law, `max_steps` budgets, `work_limit` |
| [`calendar_schedule_comparison.py`](calendar_schedule_comparison.py) | Which on-rate schedule with the same exposure leaves a two-state system most likely on? | Calendar; boundary traces, closed-form check, per-segment sensitivity |

## Optimization examples

| Script | Optimizer | Evaluation path |
| --- | --- | --- |
| [`optimize_scipy_de.py`](optimize_scipy_de.py) | SciPy `differential_evolution` (`vectorized=True`), then Nelder-Mead | one `score_matrix` call per generation, penalty for non-complete rows, traced replay |
| [`optimize_fcmaes_cmaes.py`](optimize_fcmaes_cmaes.py) | fcmaes `cmaes.Cmaes` ask/tell | one `score_matrix` call per population on two native workers, penalty for non-complete rows |
| [`optimize_fcmaes_minimize.py`](optimize_fcmaes_minimize.py) | fcmaes `bitecpp.minimize` and `cmaescpp.minimize` | scalar objective through `score_matrix([x])`, graded failure penalty outside the supported rate range |

The Rust example [`../matrix_scan.rs`](../matrix_scan.rs) scans SIR rates
with `Prepared::score_matrix` and replays the best row as a JSON record:
`cargo run --release --locked --example matrix_scan -- models/example.json`.
