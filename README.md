# population-epidemiology 0.2.0

A bounded research evaluator. See the [model card](docs/model.md) for units,
equations and supported behavior. The accompanying implementation log
records which milestones ran and which use a bounded fallback.
See the [tutorial](docs/tutorial.md) and [architecture](docs/architecture.md).
The [0.1.0 implementation record](docs/implementation-0.1.0.md) preserves the original scope.

## Rust CLI and native example

```sh
cargo run --release --locked -- describe
cargo run --release --locked -- run models/request.json
cargo run --release --locked --example sweep -- models/example.json
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
print(row["score"])
```

`Prepared.batch` crosses into native code once per batch. A set of scoped
workers evaluates candidates in isolation and returns records in input order.
The local Python campaign commits results to SQLite, preserves interrupted
attempts and reuses completed records by candidate, scenario, settings and
build identity. Fixed-work replay is scoped to the tested build.
Run `validation/run_campaign.py MODULE models/campaign.json campaign.sqlite`
with the installed package to evaluate or resume the included manifest.

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
