![population-epidemiology — stochastic population models](docs/images/banner.png)

# population-epidemiology 0.3.0

A bounded research evaluator for small stochastic population models: a
fixed-population SIR model solved by finite state projection, neutral
two-mutant Moran drift, and an exact two-state rate calendar. See the
[user manual](docs/user-manual.md) for every model field, API and status,
the [tutorial](docs/tutorial.md) for a short walk-through, the
[model card](docs/model.md) for units and equations, and the
[architecture guide](docs/architecture.md) for the implementation and
measured performance. The [0.3.0 performance log](population-epidemiology-performance-log-0.3.0.md)
records the ten performance rounds. The [0.1.0 implementation record](docs/implementation-0.1.0.md)
preserves the original scope.

## Rust CLI and native examples

```sh
cargo run --release --locked -- describe
cargo run --release --locked -- run models/request.json
cargo run --release --locked --example sweep -- models/example.json
cargo run --release --locked --example matrix_scan -- models/example.json
```

The CLI uses versioned JSON with exact model-byte SHA-256, candidate IDs and
explicit work and worker limits. `run -` accepts stdin. Exit 1 indicates a
rejected request. Exit 2 retains classified candidate records when at least
one candidate is incomplete or infeasible. Output traces are opt-in.

## Python

Install the matching wheel from `packages/` in a clean environment.
The tested wheel target is Linux x86-64 with a CPython 3.10+ ABI.
Actual clean-environment verification used CPython 3.12.

```python
from population_epidemiology import Prepared, load_example

case = load_example()
prepared = Prepared(case["model"])
row = prepared.evaluate(case["parameters"], trace=True)
assert row["status"] == "complete"
scores, statuses, steps = prepared.score_matrix([case["parameters"], [0.9, 0.5]])
print(row["score"], scores, statuses)
```

`Prepared.batch` and `Prepared.scores` cross into native code once per
batch with JSON candidates and records. `Prepared.score_matrix` takes a
list of parameter rows or a 2-D NumPy array and returns scores, statuses
and steps without JSON; use it inside optimizers. Workers evaluate
candidates in isolation and results come back in input order. The local
Python campaign commits results to SQLite, preserves interrupted attempts
and reuses completed records by candidate, scenario, settings and build
identity. Run `validation/run_campaign.py MODULE models/campaign.json campaign.sqlite`
with the installed package to evaluate or resume the included manifest.

[`examples/python/`](examples/python/README.md) has application studies
(rate landscapes, sampling designs, Moran fixation, calendar schedules) and
optimization examples with SciPy and fcmaes, with their recorded outputs.

## Performance

On a 2-vCPU host, 0.3.0 evaluates the 96-candidate SIR batch 5.6x faster
than 0.2.0 with full records and 11.1x faster through the matrix API, and
the installed SciPy search runs 2.2x to 2.9x faster end to end. See the
architecture guide's performance section and `evidence/performance-0.3.0/`.

## Verify and rebuild

```sh
cargo test --release --locked
cargo build --release --locked
python validation/validate.py ./target/release/population-epidemiology
maturin sdist --out dist
maturin build --release --locked --out dist
```

Independent numerical scripts use Python 3.12. The oscillator and population
oracles use NumPy 2.3.5 and SciPy 1.17.0. The optional
`validation/optimize.py population_epidemiology` demonstrates an in-process SciPy DE objective.
The native sweep needs only the Rust crate. See `evidence/` for actual runs.
The CRN-craft core is vendored by path because it is not on crates.io, so
the crate cannot be packaged for a registry; the source distribution and
wheel carry it. The ordered change series is in `patches/0.3.0/` (from the
0.2.0 import) and `patches/0.2.0/` (the 0.2.0 series from 0.1.0).
