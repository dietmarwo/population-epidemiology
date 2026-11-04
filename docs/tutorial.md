# Population epidemiology tutorial

Evaluate a small, fixed-population SIR model, neutral two-mutant Moran drift, and an exact two-state rate calendar. This guide uses the examples shipped with version 0.2.0.

## 1. Install and run the example

Open the package root after extracting the ZIP. The bundled wheel was
built for Linux x86-64 with a CPython 3.10+ stable ABI and tested with
CPython 3.12. Create a virtual environment with a compatible Python.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --no-deps packages/*.whl
```

If the wheel tag is incompatible, build the supplied Rust/Python source
on the target system with Rust 1.88+ and Maturin 1.9.6. Cargo needs the
locked dependencies available from its registry or cache.

```sh
python -m pip install 'maturin==1.9.6'
maturin build --release --locked --out dist
python -m pip install --no-deps dist/*.whl
```

```python
from population_epidemiology import Prepared, load_example

case = load_example()
prepared = Prepared(case["model"])
row = prepared.evaluate(case["parameters"], trace=True)
assert row["status"] == "complete", row
print(row["id"], row["steps"], row["score"])
print(list(row["metrics"]))
```

The model is validated when you create `Prepared`. A malformed model
raises `ValueError`. A well-formed batch can still contain an invalid,
infeasible, or interrupted candidate. Its row carries a status and
detail. The optimizer convention is to minimize `score`.

## 2. Run the domain workflow

The following examples continue in the same Python session.

### Fit infection and recovery

The packaged fixture starts with 9 susceptible, 2 infected and 1 recovered.
Time is measured in the model's time unit. The two candidate coordinates
are infection rate `beta` and recovery rate `gamma`, both per time unit.
The infection propensity is `beta*S*I/N` and recovery is `gamma*I`.
The objective is mean squared error on the first `train_count` sample
positions that have observations. JSON `null` masks a measurement.

```python
for beta, gamma in [(0.5, 0.4), (0.7, 0.4), (0.9, 0.5)]:
    row = prepared.evaluate([beta, gamma], trace=True)
    if row["status"] == "complete":
        metrics = row["metrics"]
        print(beta, gamma, row["score"], metrics["holdout_mse"])
        print(metrics["infected_means"], metrics["error"])
```

Keep the observation grid and training boundary fixed while searching.
After selecting rates, evaluate again with a tighter supported
`tolerance` and compare the time series, holdout loss and numerical
diagnostics. FSP state bounds in this companion follow the fixed
population. The supplied reference compares probabilities against a
dense independent generator on a small case.

### Replay neutral drift and a rate calendar

```python
neutral = Prepared({"kind": "moran", "population": 8,
                    "mutant_a": 2, "mutant_b": 1, "updates": 40})
m = neutral.evaluate([], trace=True)["metrics"]
print(m["marginal_a"], m["marginal_b"], m["mass"])
calendar = Prepared({"kind": "calendar", "initial_on_probability": 0,
                     "segments": [{"duration": 2, "on_rate": 1, "off_rate": 0},
                                  {"duration": 1, "on_rate": 0, "off_rate": 2}]})
print(calendar.evaluate([])["metrics"]["boundary_time_and_on_probability"])
```

Moran updates sample parent and death with replacement. The joint trace
lists its triangular state order. Marginals include extinction at index
0 and fixation at index `population`. Calendar rates are fixed by the
model. Change the model to change a segment rate. A zero total rate holds
the probability fixed. The calendar path can be checked against
`(1 - exp(-2)) * exp(-2)` at the final boundary.

## 3. Batch candidates and replay the winner

```python
candidates = [
    {"id": "first", "parameters": case["parameters"]},
    {"id": "second", "parameters": case["parameters"]},
]
scores = prepared.scores(candidates, workers=2)["records"]
assert [row["id"] for row in scores] == ["first", "second"]
winner = min((r for r in scores if r["status"] == "complete"),
             key=lambda r: r["score"])
candidate = next(c for c in candidates if c["id"] == winner["id"])
replay = prepared.evaluate(candidate["parameters"], trace=True)
assert replay["status"] == "complete", replay
print(replay["score"], replay["metrics"])
```

Replace the duplicate vectors with distinct choices for a real search.
Use `batch` if you need the full metrics for every row. `scores` returns
ordered records with score, status, steps and detail, and `metrics: null`.
It rejects any candidate with `trace: true`. A one-worker batch or a
singleton runs inline. Larger batches can try two or four workers.
Measure before setting a worker count for a real workload.

The optional differential-evolution driver times preparation, seeded
search and traced replay as one run. Install SciPy for this example.

```sh
python -m pip install scipy
python validation/optimize_performance.py population_epidemiology --scores
python validation/optimize_performance.py population_epidemiology
```

The driver uses small local bounds and a short fixed search budget.
For a scientific result, use problem-specific bounds, an explicit
feasibility penalty or constraint policy, and a held-out evaluation.

## 4. Run and resume a local campaign

```sh
python validation/run_campaign.py population_epidemiology models/campaign.json campaign.sqlite > first.json
python validation/run_campaign.py population_epidemiology models/campaign.json campaign.sqlite > replay.json
cmp first.json replay.json
```

The manifest holds scenario models and candidate vectors. Keep IDs
unique inside each list. Reuse the database to resume identical work.
The native build, model, scenario, candidate and work settings are
part of cache identity. A new wheel re-evaluates cached work under
its own build identity. A cancel file passed as `cancel_file` or a
`wall_ms` budget can stop new evaluations cooperatively. Inspect
statuses in the attempt journal before rerunning.

## 5. Use the Rust CLI and benchmark

```sh
cargo build --release --locked --bins --examples
target/release/population-epidemiology describe
target/release/population-epidemiology run models/request.json
target/release/examples/sweep models/example.json
target/release/examples/performance models/example.json serial
CRAFT_SCORES=1 target/release/examples/performance models/example.json serial
```

The CLI reads a request file or standard input with `run -`. Exit code
0 indicates complete rows, 1 indicates a rejected request and 2
indicates one or more unsuccessful rows in a well-formed request.
The benchmark prints elapsed seconds per prepared batch after warm-up.
Use `parallel` for four workers, `trace` for traced evaluations,
or `domain` and `domain2` for the package's supplemental cases.
With a previous benchmark executable built from the baseline and this
one, `validation/performance_pair.py BEFORE AFTER models/example.json serial`
checks output parity and times both. See the included performance log
and raw evidence for the exact workloads and comparison boundaries.

## 6. Verify results

```sh
cargo test --release --locked
cargo clippy --locked --all-targets --features python -- -D warnings
python validation/validate.py target/release/population-epidemiology
```

Consult the [model card](model.md) for supported inputs, the
[architecture guide](architecture.md) for data flow and numerical
assumptions, and the [source notes](sources.md) for references.
