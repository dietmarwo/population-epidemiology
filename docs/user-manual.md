# population-epidemiology user manual

Version 0.3.0. This manual is for researchers and engineers who evaluate,
fit or screen small stochastic population models with this package. It
covers installation, the model and request formats, every status and
error, budgets, the Python, Rust and command-line interfaces,
reproducibility and performance. The [model card](model.md) states the
equations and research scope, the [tutorial](tutorial.md) is a short
walk-through, and the [architecture guide](architecture.md) describes the
implementation and its measured performance.

Every Python block in this manual runs as written, in order, in one
Python session; `validation/release_checks.py` executes them for each
release. Shell transcripts were copied from real runs of the 0.3.0 CLI.

## Contents

1. [Install](#1-install)
2. [Concepts](#2-concepts)
3. [Model reference](#3-model-reference)
4. [Requests and results](#4-requests-and-results)
5. [Statuses and errors](#5-statuses-and-errors)
6. [Budgets](#6-budgets)
7. [Python API](#7-python-api)
8. [Campaigns](#8-campaigns)
9. [Rust API](#9-rust-api)
10. [Command line](#10-command-line)
11. [Reproducibility and provenance](#11-reproducibility-and-provenance)
12. [Performance](#12-performance)
13. [Troubleshooting](#13-troubleshooting)
14. [Example programs](#14-example-programs)

## 1. Install

### Wheel

The release ZIP ships a wheel for Linux x86-64 with the CPython 3.10+
stable ABI. It was tested with CPython 3.12. The package itself has no
Python dependencies. NumPy is optional: `score_matrix` accepts NumPy
arrays. The optimization examples also use SciPy and fcmaes.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --no-deps packages/population_epidemiology-0.3.0-cp310-abi3-manylinux_2_34_x86_64.whl
python -m pip install numpy==2.3.5 scipy==1.17.0 fcmaes   # optional, for the examples
```

### From source

Building needs Rust 1.88 or newer and maturin 1.9.6. The 0.3.0 source
distribution contains the pinned CRN-craft core under `vendor/`, so no
other Rust source is fetched. Cargo still needs the locked crates.io
dependencies from its registry or cache. A checkout of `main` instead takes
`crn-craft` 0.5.0 from crates.io with the other locked dependencies.

```sh
python -m pip install maturin==1.9.6
tar xzf packages/population_epidemiology-0.3.0.tar.gz
cd population_epidemiology-0.3.0
maturin build --release --locked --out dist
python -m pip install --no-deps dist/*.whl
```

### Command-line tool

```sh
cargo build --release --locked --bins --examples
target/release/population-epidemiology describe
```

Up to 0.3.0 the crate could not be packaged for crates.io, because it
depended on a vendored CRN-craft core by path. A checkout of `main` depends
on `crn-craft` 0.5.0 from crates.io instead, so `cargo package` works there;
the crate is not published.

## 2. Concepts

- **Model.** A JSON object with a `kind` of `sir`, `moran` or `calendar`
  and that kind's fields. `Prepared(model)` validates it once and builds
  everything that does not depend on a candidate: the SIR state space and
  generator structure, the Moran transition stencil, or the calendar's
  exact boundary probabilities.
- **Candidate.** An identifier, a parameter vector and an optional
  `trace` flag. SIR candidates have two parameters, `[beta, gamma]`.
  Moran and calendar candidates have none, `[]`, because their model fixes
  every rate.
- **Batch.** Up to 4,096 candidates evaluated in one native call. Records
  come back in input order, whatever the worker count.
- **Record.** `id`, `status`, `steps`, `score`, `metrics` and `detail`.
  Optimizers minimize `score`. It is set only when the status is
  `complete`.
- **Logical steps.** A deterministic work counter defined by the domain.
  `max_steps` limits it per candidate. The count does not depend on
  workers, timing or which internal kernel ran, so fixed-work results
  replay exactly.
- **Three result paths.** Full records (`batch`, `evaluate`) carry
  metrics and optional traces. Lean records (`scores`) carry the same
  keys with `metrics: null`. The matrix path (`score_matrix`) returns only
  scores, statuses and steps, without JSON in either direction. All three
  give identical statuses, steps and scores, bit for bit.

```python
from population_epidemiology import Prepared, load_example

case = load_example()
prepared = Prepared(case["model"])
row = prepared.evaluate(case["parameters"])
print(row["id"], row["status"], row["steps"], row["score"])
lean = prepared.scores([{"id": "x", "parameters": case["parameters"]}])["records"][0]
scores, statuses, steps = prepared.score_matrix([case["parameters"]])
assert (lean["score"], lean["steps"]) == (row["score"], row["steps"]) == (scores[0], steps[0])
assert lean["metrics"] is None and statuses == ["complete"]
```

## 3. Model reference

All fields are required and unknown fields are rejected, so there are no
defaults. Times are in the model's own time unit, and rates are per time
unit. A model that violates a bound raises `ValueError` in Python
(exit code 1 from the CLI) with the message shown for its kind.

### 3.1 SIR (`"kind": "sir"`)

A closed population of N = S + I + R individuals with stochastic counts.
Infection S + I → 2I has propensity `beta * S * I / N`, and recovery
I → R has propensity `gamma * I`. The CRN-craft finite state
projection (FSP) gives the exact count distribution at each sample time
up to its documented series tail.

| Field | Type | Bounds | Meaning |
| --- | --- | --- | --- |
| `susceptible` | integer | ≥ 0 | initial S |
| `infected` | integer | ≥ 0 | initial I |
| `recovered` | integer | ≥ 0 | initial R; N = S + I + R must be 1–30 |
| `sample_times` | array of numbers | 2–64 entries, each in [0, 120], strictly increasing | observation times; CRN-craft requires the first to be 0 |
| `observations` | array of numbers or `null` | same length as `sample_times`, each in [0, N] | observed infected counts; `null` masks a point |
| `train_count` | integer | 1 ≤ `train_count` < number of sample times | the first `train_count` points train the fit, the rest are held out |
| `tolerance` | number | [1e-13, 1e-5] | FSP Poisson-series tail tolerance |

At least one of the first `train_count` observations must be non-null.
Violations raise `invalid SIR population, time grid, observations or
tolerance`. A model whose first sample time is not 0 passes this check,
but CRN-craft rejects its documents, so every candidate is `infeasible`
(section 5).

**Candidate parameters:** `[beta, gamma]`, each in [0, 5] per time unit.

**Score:** mean squared error between the expected infected count and the
non-null training observations.

**Metrics (full records):**

| Key | Meaning |
| --- | --- |
| `infected_means` | expected infected count at each sample time |
| `holdout_mse` | mean squared error over non-null held-out points, or `null` |
| `error` | per sample time: `physical_leakage`, `series_tail_upper_bound`, `mass_defect`, `roundoff_mass_residual` |
| `fsp_work` | CRN-craft's FSP work count |
| `numerical_error_certified` | always `false`; the series tail is a bound, roundoff is not certified |
| `crn_model_sha256`, `crn_run_sha256` | SHA-256 of the CRN-craft model and run documents for this rate binding |
| `fsp` | traced candidates only: the full CRN-craft FSP result (states, per-time probabilities, moments and diagnostics); otherwise `null` |

**Steps:** 2 per candidate, one before and one after the solve. With
`max_steps: 1` an SIR candidate stops with `work_limit` before solving.

```python
sir = {"kind": "sir", "susceptible": 9, "infected": 2, "recovered": 1,
       "sample_times": [0.0, 1.0, 2.0, 3.0, 5.0],
       "observations": [2.0, 2.16, 2.14, None, 1.44],
       "train_count": 3, "tolerance": 1e-11}
fit = Prepared(sir)
row = fit.evaluate([0.7, 0.4])
m = row["metrics"]
print(row["score"], m["holdout_mse"])
print([round(x, 4) for x in m["infected_means"]])
print(sorted(m["error"][1]), m["fsp_work"], m["fsp"])
assert m["fsp"] is None and len(m["infected_means"]) == 5

# A traced candidate adds CRN-craft's full FSP result.
fsp = fit.evaluate([0.7, 0.4], trace=True)["metrics"]["fsp"]
print(fsp["schema"], len(fsp["states"]), fsp["states"][:3], fsp["uniformization_rate"])
snapshot = fsp["snapshots"][-1]
print(snapshot["time"], round(sum(snapshot["probabilities"]), 12), snapshot["projected_first_moment"])
```

Masking every held-out point gives `holdout_mse: null`. Tightening the
tolerance shrinks `series_tail_upper_bound` at the cost of more FSP work:

```python
no_holdout = Prepared(sir | {"observations": [2.0, 2.16, 2.14, None, None]})
assert no_holdout.evaluate([0.7, 0.4])["metrics"]["holdout_mse"] is None
for tolerance in (1e-5, 1e-13):
    m = Prepared(sir | {"tolerance": tolerance}).evaluate([0.7, 0.4])["metrics"]
    print(tolerance, max(e["series_tail_upper_bound"] for e in m["error"]), m["fsp_work"])
```

Rates of exactly 0 are allowed. They change which states are reachable,
so the package hands them to the CRN-craft solver instead of its
prepared projection; the result is the same either way.

```python
for rates in ([0.0, 0.4], [0.7, 0.0], [5.0, 5.0]):
    print(rates, fit.evaluate(rates)["status"], fit.evaluate(rates)["score"])
```

Models outside the bounds are rejected when you prepare them:

```python
for bad in (sir | {"susceptible": 40},              # N above 30
            sir | {"tolerance": 1e-3},              # tolerance above 1e-5
            sir | {"train_count": 5},               # no held-out position left
            sir | {"sample_times": [0, 2, 1, 3, 5]},  # not increasing
            sir | {"observations": [None, None, None, 1, 1]},  # no training point
            sir | {"extra": 1}):                    # unknown field
    try:
        Prepared(bad)
    except ValueError as error:
        print(type(error).__name__, str(error)[:70])
    else:
        raise AssertionError(bad)
```

### 3.2 Neutral two-mutant Moran (`"kind": "moran"`)

A population of N individuals of three types: wild type, mutant A and
mutant B. Each discrete update draws a uniform parent and a uniform
individual to replace, with replacement, so self-replacement and null
updates are possible. No type has a fitness advantage. The package
propagates the exact joint law over the (N+1)(N+2)/2 states (a, b).

| Field | Type | Bounds | Meaning |
| --- | --- | --- | --- |
| `population` | integer | 2–32 | N |
| `mutant_a` | integer | ≥ 0 | initial A count |
| `mutant_b` | integer | ≥ 0, `mutant_a + mutant_b` ≤ N | initial B count |
| `updates` | integer | 0–10,000 | number of update attempts |

Violations raise `invalid Moran state or update count`.

**Candidate parameters:** none, `[]`.

**Score:** P(a = N) + P(b = N), the probability that one of the mutants
has fixed.

**Metrics:** `marginal_a` and `marginal_b` (length N+1, from extinction at
index 0 to fixation at index N), `mass` (total probability, 1 up to
roundoff), `updates`, and `joint`, which is `null` unless the candidate is
traced and then holds `states` (the triangular (a, b) order) and
`probabilities`.

**Steps:** (N+1)(N+2)/2 × (updates + 1). Long runs with large populations
need a larger `max_steps` than the default 1,000,000.

```python
moran = {"kind": "moran", "population": 8, "mutant_a": 2, "mutant_b": 1, "updates": 40}
drift = Prepared(moran)
row = drift.evaluate([], trace=True)
m = row["metrics"]
print(row["steps"], round(row["score"], 6), round(m["mass"], 12))
print([round(x, 4) for x in m["marginal_a"]])
print(len(m["joint"]["states"]), m["joint"]["states"][:4])
assert row["steps"] == 45 * 41

long = Prepared(moran | {"population": 32, "mutant_a": 8, "mutant_b": 8, "updates": 4000})
print(long.evaluate([])["status"], long.evaluate([], max_steps=3_000_000)["status"])
print(drift.evaluate([1.0])["status"], drift.evaluate([1.0])["detail"])
```

The neutral fixation probability of each mutant tends to its initial
frequency:

```python
far = Prepared({"kind": "moran", "population": 6, "mutant_a": 2, "mutant_b": 1,
                "updates": 2000}).evaluate([])["metrics"]
print(round(far["marginal_a"][6], 9), round(far["marginal_b"][6], 9))
assert abs(far["marginal_a"][6] - 2 / 6) < 1e-9 and abs(far["marginal_b"][6] - 1 / 6) < 1e-9
```

### 3.3 Two-state rate calendar (`"kind": "calendar"`)

A two-state (off/on) continuous-time chain whose rates are constant within
each segment and change at segment boundaries. The probability of being
on is propagated exactly with `p += (p_eq - p) * (1 - exp(-(on + off) * duration))`,
where `p_eq = on / (on + off)`. A segment whose two rates are both zero
holds the probability.

| Field | Type | Bounds | Meaning |
| --- | --- | --- | --- |
| `initial_on_probability` | number | [0, 1] | probability of the on state at time 0 |
| `segments` | array of objects | 1–1,024 segments | in time order |
| `segments[].duration` | number | [1e-10, 120] | segment length |
| `segments[].on_rate` | number | [0, 1000] | off → on rate in this segment |
| `segments[].off_rate` | number | [0, 1000] | on → off rate in this segment |

Violations raise `invalid calendar probability, duration or rates`.

**Candidate parameters:** none, `[]`. To change a rate, change the model.

**Score:** the on-probability at the final boundary. Negate it in your
objective if you want to maximize it.

**Metrics:** `boundary_time_and_on_probability`, one `[time, probability]`
pair per segment end.

**Steps:** one per segment.

```python
import math

calendar = Prepared({"kind": "calendar", "initial_on_probability": 0.0,
                     "segments": [{"duration": 2, "on_rate": 1, "off_rate": 0},
                                  {"duration": 1, "on_rate": 0, "off_rate": 2}]})
row = calendar.evaluate([])
print(row["steps"], row["metrics"]["boundary_time_and_on_probability"])
assert abs(row["score"] - (1 - math.exp(-2)) * math.exp(-2)) < 1e-15
held = Prepared({"kind": "calendar", "initial_on_probability": 0.25,
                 "segments": [{"duration": 5, "on_rate": 0, "off_rate": 0}]})
assert held.evaluate([])["score"] == 0.25
print(calendar.evaluate([0.5])["detail"])
```

## 4. Requests and results

### Request (CLI and `population_epidemiology::run`)

```json
{
  "schema": "population-epidemiology.request/0.1",
  "model": {"kind": "calendar", "initial_on_probability": 0,
            "segments": [{"duration": 1, "on_rate": 0.3, "off_rate": 0}]},
  "candidates": [{"id": "a", "parameters": [], "trace": false}],
  "budget": {"workers": 1, "max_steps": 1000000}
}
```

| Field | Type | Rules |
| --- | --- | --- |
| `schema` | string | exactly `population-epidemiology.request/0.1` |
| `model` | object | section 3; its exact bytes are hashed into `model_sha256` |
| `candidates` | array | 1–4,096 entries |
| `candidates[].id` | string | 1–128 bytes, unique within the batch |
| `candidates[].parameters` | array of numbers | at most 64 values |
| `candidates[].trace` | boolean | optional, default `false`; at most 16 traced per batch |
| `budget` | object | optional; section 6 |

Python's `batch` and `scores` take the `candidates` list and the budget as
keyword arguments; `Prepared` holds the model.

### Result

| Field | Meaning |
| --- | --- |
| `model_sha256` | SHA-256 of the model bytes the evaluator was prepared with |
| `records` | one record per candidate, in input order |
| `schema` | `population-epidemiology.result/0.1` |
| `version` | package version, `0.3.0` |
| `records[].id` | the candidate's identifier |
| `records[].status` | section 5 |
| `records[].steps` | logical steps charged |
| `records[].score` | number when complete, otherwise `null` |
| `records[].metrics` | kind-specific object in full records when complete; `null` in lean records and for failures |
| `records[].detail` | `null` when complete; otherwise a reason |

The request and result schema tags stay at 0.1 in 0.3.0; no field
changed. The matrix API is not part of the JSON schema.

```python
import hashlib
import json

result = prepared.batch([{"id": "fit", "parameters": [0.7, 0.4]},
                         {"id": "slow", "parameters": [0.2, 0.1]}], workers=2)
print(sorted(result), result["schema"], result["version"])
print([(r["id"], r["status"], r["steps"]) for r in result["records"]])
print(sorted(result["records"][0]))
```

## 5. Statuses and errors

### Record statuses

| Status | Meaning | `steps` | Typical `detail` |
| --- | --- | --- | --- |
| `complete` | evaluated; `score` is finite | logical steps used | `null` |
| `invalid` | the candidate does not fit the model | steps used before the check (usually 0) | `expected beta and gamma per time unit`, `rates must be in [0,5]`, `neutral Moran has no candidate parameters`, `calendar rates come from the model`, `nonfinite objective` |
| `infeasible` | CRN-craft rejected or could not finish the FSP | 2 for SIR | the CRN-craft diagnostics (warnings and the error) as a JSON array string |
| `work_limit` | `max_steps` reached | equal to `max_steps` | `candidate step budget reached` |
| `wall_limit` | the batch deadline (`wall_ms`) passed | steps used so far | `batch deadline reached` |
| `cancelled` | the cancel file existed at a step check | steps used so far | `cancel file exists` |

Only `complete` rows have a score. A batch is still well formed when some
rows fail. The CLI then exits with code 2 and prints every record.

```python
import tempfile
from pathlib import Path

statuses = {}
statuses["complete"] = prepared.evaluate([0.7, 0.4])
statuses["invalid"] = prepared.evaluate([7.0, 0.4])
statuses["work_limit"] = prepared.evaluate([0.7, 0.4], max_steps=1)
late_start = Prepared(sir | {"sample_times": [1.0, 2.0, 3.0, 4.0, 5.0]})
statuses["infeasible"] = late_start.evaluate([0.7, 0.4])
slow = Prepared({"kind": "moran", "population": 32, "mutant_a": 8, "mutant_b": 8,
                 "updates": 10_000})
statuses["wall_limit"] = slow.evaluate([], wall_ms=1, max_steps=10_000_000)
with tempfile.TemporaryDirectory() as folder:
    flag = Path(folder) / "stop"
    flag.touch()
    statuses["cancelled"] = prepared.evaluate([0.7, 0.4], cancel_file=str(flag))
for name, row in statuses.items():
    assert row["status"] == name, row
    print(f"{name:11s} steps={row['steps']:<6} detail={str(row['detail'])[:60]}")
diagnostics = json.loads(statuses["infeasible"]["detail"])
error = next(d for d in diagnostics if d["severity"] == "error")
print(error["code"], error["path"], error["message"])
```

### Batch errors

These reject the whole call: `ValueError` in Python, exit code 1 and
`{"status": "rejected", "detail": ...}` on standard error from the CLI.

| Message | Cause |
| --- | --- |
| `invalid budget or batch size` | workers outside 1–64, `max_steps` outside 1–100,000,000, `wall_ms` of 0 or above 86,400,000, or 0 or more than 4,096 candidates |
| `invalid or duplicate candidate id, or parameter vector exceeds limit` | empty, too long or repeated ID, or more than 64 parameters |
| `at most 16 traced candidates per batch` | more than 16 `trace: true` |
| `scores mode accepts untraced candidates only` | `scores` with a traced candidate |
| `unsupported request schema` | CLI request with another `schema` value |
| `request exceeds 16 MiB` | CLI input above 16 MiB |
| `matrix rows must be in 1..=4096` | `score_matrix` with no rows or too many |
| `parameter count is not a multiple of the row count` | Rust `score_matrix` with a ragged flat buffer |
| `parameter vector exceeds limit` | `score_matrix` rows longer than 64 |
| `ragged parameter rows` | Python `score_matrix` rows of different lengths |
| `score_matrix expects a 2-D array of parameter rows` | a NumPy array that is not 2-D |
| serde messages such as `unknown field ...` or `missing field ...` | malformed model, candidate or budget JSON |
| `invalid SIR ...`, `invalid Moran ...`, `invalid calendar ...` | model bounds, section 3 |

```python
def rejected(call):
    try:
        call()
    except ValueError as error:
        return str(error)
    raise AssertionError("accepted")

print(rejected(lambda: prepared.batch([], workers=1)))
print(rejected(lambda: prepared.batch([{"id": "a", "parameters": []}] * 2)))
print(rejected(lambda: prepared.batch([{"id": "a", "parameters": []}], workers=65)))
print(rejected(lambda: prepared.scores([{"id": "a", "parameters": [0.7, 0.4], "trace": True}])))
print(rejected(lambda: prepared.score_matrix([[0.7, 0.4], [0.7]])))
print(rejected(lambda: prepared.batch([{"id": "a", "parameters": [], "extra": 1}])))
print(rejected(lambda: prepared.batch([{"id": "a", "parameters": []}], unknown=1)))
```

## 6. Budgets

Budgets apply per batch in `batch`, `scores`, `score_matrix`, `evaluate`
and the CLI's `budget` object. Every field is optional.

| Field | Type | Default | Bounds | Effect |
| --- | --- | --- | --- | --- |
| `workers` | integer | 1 | 1–64 | number of scoped worker threads; records stay in input order and results do not depend on it |
| `max_steps` | integer | 1,000,000 | 1–100,000,000 | logical steps per candidate; reaching it gives `work_limit` |
| `wall_ms` | integer | none | 1–86,400,000 | batch deadline in milliseconds; checked at each logical step, giving `wall_limit` |
| `cancel_file` | string path | none | any path | if the file exists at a step check, the candidate stops with `cancelled` |

Deadlines and cancellation are cooperative: they are checked at logical
steps, so a running step finishes first. SIR's FSP solve also polls them.
A batch with `wall_ms` or `cancel_file` is *monitored*. Monitored SIR
candidates run through the CRN-craft solver, which polls the
deadline and cancel file inside the solve; this is about ten times slower
per candidate than the prepared projection (section 12) but gives the
same results. Moran and calendar candidates check every logical step when
monitored. Unmonitored work can only stop at `max_steps`, so its step
count is charged in bulk.

```python
cases = [{"id": f"c{i}", "parameters": [0.6 + 0.02 * i, 0.4]} for i in range(8)]
one = prepared.scores(cases, workers=1)["records"]
four = prepared.scores(cases, workers=4)["records"]
assert one == four
with tempfile.TemporaryDirectory() as folder:
    absent = str(Path(folder) / "not-created")
    monitored = prepared.scores(cases, cancel_file=absent, wall_ms=60_000)["records"]
assert monitored == one
limited = Prepared(moran).scores([{"id": "m", "parameters": []}], max_steps=100)["records"][0]
print(limited["status"], limited["steps"])
```

## 7. Python API

```python
from population_epidemiology import Prepared, load_example
```

| Call | Returns |
| --- | --- |
| `Prepared(model)` | a prepared evaluator; raises `ValueError` for an invalid model |
| `Prepared.batch(candidates, **budget)` | the result dictionary with full records |
| `Prepared.scores(candidates, **budget)` | the result dictionary with `metrics: null`; rejects traced candidates |
| `Prepared.score_matrix(rows, **budget)` | `(scores, statuses, steps)`, three lists in row order |
| `Prepared.evaluate(parameters, *, trace=False, **budget)` | one full record with ID `"one"` |
| `load_example()` | `{"model": ..., "parameters": ...}`, the packaged SIR fixture and its generating rates |

The model and candidates are encoded to JSON with sorted keys and finite
numbers (`allow_nan=False`), so NaN or infinity in `batch` or `scores`
raises `ValueError` before native code runs. Each call crosses into native
code once and releases the GIL during evaluation, so other Python threads
keep running.

### `score_matrix`

`rows` is a list of equal-length lists or a 2-D NumPy array, one row per
candidate. Rows are numbered by position, so there are no IDs and no
duplicate check. A score is `nan` unless the status is `complete`.
Statuses, steps and scores equal what `scores` reports for the same
vectors, bit for bit. Values JSON cannot carry, such as NaN, are allowed
and classified per row.

```python
import numpy as np

rows = np.array([[0.7, 0.4], [0.9, 0.5], [6.0, 0.4]])
scores, statuses, steps = prepared.score_matrix(rows, workers=2)
print([round(s, 6) for s in scores[:2]], statuses, steps)
assert np.isnan(scores[2])
listed = prepared.scores([{"id": str(i), "parameters": r} for i, r in enumerate(rows.tolist())])
assert [r["status"] for r in listed["records"]] == statuses
assert listed["records"][0]["score"] == scores[0]
print(prepared.score_matrix([[float("nan"), 0.4]])[1])
# Dimension-0 rows for models without candidate parameters.
print(Prepared(moran).score_matrix(np.zeros((3, 0)))[0])
```

### Choosing an API

```python
from time import perf_counter

grid = [[b, g] for b in np.linspace(0.4, 1.0, 20).tolist() for g in np.linspace(0.2, 0.6, 20).tolist()]
t0 = perf_counter()
matrix = prepared.score_matrix(grid)
t1 = perf_counter()
lean = prepared.scores([{"id": str(i), "parameters": r} for i, r in enumerate(grid)])
t2 = perf_counter()
full = prepared.batch([{"id": str(i), "parameters": r} for i, r in enumerate(grid)])
t3 = perf_counter()
print(f"400 candidates: matrix {1e3 * (t1 - t0):.1f} ms, lean {1e3 * (t2 - t1):.1f} ms, "
      f"full {1e3 * (t3 - t2):.1f} ms")
best = int(np.nanargmin(matrix[0]))
print("best", grid[best], prepared.evaluate(grid[best])["metrics"]["holdout_mse"])
```

Use `score_matrix` inside optimizers, `scores` when you need IDs or
`detail` strings for screening, and `batch`/`evaluate` with `trace=True`
to inspect a chosen candidate.

### An optimizer objective

The usual pattern maps non-complete rows to a penalty. With SciPy's
vectorized differential evolution, one `score_matrix` call scores a whole
generation:

```python
from scipy.optimize import differential_evolution

def objective(x):
    scores, statuses, _ = prepared.score_matrix(np.atleast_2d(x.T))
    return np.where(np.array(statuses) == "complete", scores, 1e6)

found = differential_evolution(objective, [(0.0, 5.0), (0.0, 5.0)], seed=3, popsize=8,
                               maxiter=40, polish=False, vectorized=True,
                               updating="deferred", tol=0)
print(np.round(found.x, 3), f"{found.fun:.2e}")
assert np.allclose(found.x, [0.7, 0.4], atol=0.05)
```

## 8. Campaigns

`population_epidemiology.campaign.run(database, scenarios, candidates, *,
workers=1, chunk_size=16, **budget)` evaluates every candidate against
every scenario model and returns rows in manifest order. Each row holds
`scenario`, `key`, `build_sha256`, `model_sha256` and the full `record`.
Results are committed to a SQLite database chunk by chunk:

- `results` keeps complete, invalid and infeasible records. They depend
  only on the build, scenario, model, candidate and work settings, which
  form the SHA-256 cache key.
- `attempts` is an append-only journal of every evaluated record,
  including cancelled and time-limited ones, which are retried on the next
  run.
- `build_sha256` hashes the installed native extension, so a new build
  re-evaluates instead of reusing old records.
- `cancel_file` is a runtime control and not part of the key. `chunk_size`
  must be 1–256, and scenario and candidate IDs must be unique.

```python
import sqlite3

from population_epidemiology import campaign

scenarios = [{"id": "baseline", "model": case["model"]},
             {"id": "longer", "model": case["model"] | {"tolerance": 1e-13}}]
trial = [{"id": f"r{i}", "parameters": [0.5 + 0.1 * i, 0.4]} for i in range(5)]
with tempfile.TemporaryDirectory() as folder:
    db = Path(folder) / "campaign.sqlite"
    flag = Path(folder) / "stop"
    flag.touch()
    stopped = campaign.run(db, scenarios, trial, workers=2, cancel_file=str(flag))
    print({r["record"]["status"] for r in stopped})
    flag.unlink()
    first = campaign.run(db, scenarios, trial, workers=2)
    again = campaign.run(db, scenarios, trial, workers=1)
    assert first == again
    with sqlite3.connect(db) as conn:
        print(conn.execute("SELECT count(*) FROM results").fetchone()[0],
              conn.execute("SELECT count(*) FROM attempts").fetchone()[0])
print(first[0]["scenario"], first[0]["record"]["id"], first[0]["record"]["status"])
```

The included manifest runner evaluates `models/campaign.json`:

```sh
python validation/run_campaign.py population_epidemiology models/campaign.json campaign.sqlite > first.json
python validation/run_campaign.py population_epidemiology models/campaign.json campaign.sqlite > replay.json
cmp first.json replay.json
```

## 9. Rust API

The library crate `population_epidemiology` exposes:

| Item | Purpose |
| --- | --- |
| `Prepared::new(model_json: &str) -> Result<Prepared, String>` | validate and prepare a model |
| `Prepared::evaluate_json(candidates, budget) -> Result<String, String>` | full records as the result JSON |
| `Prepared::scores_json(candidates, budget) -> Result<String, String>` | lean records, untraced candidates only |
| `Prepared::score_matrix(rows, parameters, budget) -> Result<Vec<MatrixScore>, String>` | typed scores for row-major parameters; `parameters.len() / rows` values per row, dimension 0 allowed |
| `Prepared::model_sha256() -> &str` | SHA-256 of the prepared model bytes |
| `MatrixScore { status: &'static str, steps: u64, score: Option<f64> }` | one matrix row; `score` is `Some` only for `"complete"` |
| `run(request_json) -> Result<serde_json::Value, String>` | the CLI's request path |
| `capabilities() -> serde_json::Value` | the `describe` output |
| `Candidate` | the deserializable candidate type (`id`, `parameters`, `trace`) |

Budgets are passed as JSON strings, for example `r#"{"workers":2}"#`.

```rust
use population_epidemiology::Prepared;

let prepared = Prepared::new(&model_json)?;
// Two rows of (beta, gamma).
let rows = prepared.score_matrix(2, &[0.7, 0.4, 0.9, 0.5], r#"{"workers":2}"#)?;
for row in &rows {
    println!("{} {} {:?}", row.status, row.steps, row.score);
}
```

[`examples/matrix_scan.rs`](../examples/matrix_scan.rs) is a complete
program: `cargo run --release --locked --example matrix_scan -- models/example.json`.
`examples/sweep.rs` evaluates five scaled vectors with two workers.
`examples/performance.rs` and `examples/breakdown.rs` are the benchmark
programs described in section 12.

## 10. Command line

```text
population-epidemiology describe
population-epidemiology run FILE
population-epidemiology run -        # read the request from standard input
```

| Exit code | Meaning |
| --- | --- |
| 0 | every record is `complete` |
| 1 | the request was rejected; one JSON line on standard error |
| 2 | the request was valid but at least one record is not `complete`; all records are printed |

```sh
$ target/release/population-epidemiology describe
{"batch":true,"models":["sir_via_crn_fsp","finite_neutral_two_mutant_moran","exact_two_state_rate_calendar"],"schema":"population-epidemiology.capabilities/0.1","version":"0.3.0","worker_limit":64}
```

```sh
$ echo '{"schema":"population-epidemiology.request/0.1","model":{"kind":"calendar","initial_on_probability":0,"segments":[{"duration":1,"on_rate":0.3,"off_rate":0},{"duration":2,"on_rate":0,"off_rate":0.4}]},"candidates":[{"id":"a","parameters":[]}]}' | target/release/population-epidemiology run -
{"model_sha256":"201430af95cef942c385efdb17403f88a39a9392635e819ddcfce4bfb40cbcab","records":[{"detail":null,"id":"a","metrics":{"boundary_time_and_on_probability":[[1.0,0.2591817793182821],[3.0,0.11645788041914201]]},"score":0.11645788041914201,"status":"complete","steps":2}],"schema":"population-epidemiology.result/0.1","version":"0.3.0"}
$ echo $?
0
```

A valid request with an invalid row exits with 2:

```sh
$ echo '{"schema":"population-epidemiology.request/0.1","model":{"kind":"moran","population":3,"mutant_a":1,"mutant_b":1,"updates":2},"candidates":[{"id":"ok","parameters":[]},{"id":"bad","parameters":[0.5]}],"budget":{"max_steps":100}}' | target/release/population-epidemiology run -
{"model_sha256":"6583f5c626c0b76d99c9c525b689ddc2ad546216a3879b42f340268ab6016b67","records":[{"detail":null,"id":"ok","metrics":{"joint":null,"marginal_a":[0.345679012345679,0.35802469135802467,0.24691358024691357,0.04938271604938271],"marginal_b":[0.345679012345679,0.35802469135802467,0.24691358024691357,0.04938271604938271],"mass":0.9999999999999999,"updates":2},"score":0.09876543209876543,"status":"complete","steps":30},{"detail":"neutral Moran has no candidate parameters","id":"bad","metrics":null,"score":null,"status":"invalid","steps":0}],"schema":"population-epidemiology.result/0.1","version":"0.3.0"}
$ echo $?
2
```

A rejected request exits with 1:

```sh
$ echo '{"schema":"population-epidemiology.request/9","model":{"kind":"moran","population":3,"mutant_a":1,"mutant_b":1,"updates":2},"candidates":[{"id":"a","parameters":[]}]}' | target/release/population-epidemiology run -
{"detail":"unsupported request schema","status":"rejected"}
$ echo $?
1
$ target/release/population-epidemiology
{"detail":"usage: population-epidemiology describe | run FILE (use - for stdin)","status":"rejected"}
```

`models/request.json` is a traced SIR request: `run models/request.json`
exits with 0 and prints one record that includes the full FSP result. The
CLI prints its JSON with sorted keys.

## 11. Reproducibility and provenance

- `model_sha256` is the SHA-256 of the exact model bytes. The Python
  wrapper encodes with `json.dumps(model, sort_keys=True,
  separators=(",", ":"), allow_nan=False)`, so the same dictionary gives
  the same hash. The CLI hashes the model text as written in the
  request, so whitespace or key order there changes the hash.
- SIR full records add `crn_model_sha256` and `crn_run_sha256`, the
  hashes of the CRN-craft documents for each rate binding. They are the
  same whether the prepared projection or the CRN-craft solver ran.
- Results do not depend on the worker count or on which SIR kernel ran.
  Fixed-work replays (no deadline) are deterministic for a given build.
- `version` in every result and `build_sha256` in campaign rows identify
  the code that produced a record.
- The CRN-craft version is pinned (`=0.5.0`); see the [source notes](sources.md).

```python
encoded = json.dumps(case["model"], sort_keys=True, separators=(",", ":"), allow_nan=False)
digest = hashlib.sha256(encoded.encode()).hexdigest()
assert digest == prepared.batch([{"id": "h", "parameters": [0.7, 0.4]}])["model_sha256"]
rows = [{"id": str(i), "parameters": [0.5 + 0.05 * i, 0.4]} for i in range(6)]
assert prepared.batch(rows, workers=1) == prepared.batch(rows, workers=3)
print(digest[:16], prepared.evaluate([0.7, 0.4])["metrics"]["crn_model_sha256"][:16])
```

## 12. Performance

Measured on the 0.3.0 release host: 2 vCPUs (Intel Xeon Processor @
2.10GHz, KVM), Rust 1.97 and CPython 3.12. Numbers come from
`examples/breakdown.rs`, saved as
`evidence/performance-0.3.0/breakdown.json`, as medians of 41 alternating
rounds. The [architecture guide](architecture.md#performance) has the
full 0.2.0 to 0.3.0 comparison.

| One call | Matrix | Lean | Full |
| --- | ---: | ---: | ---: |
| 1 SIR candidate (12 people, 5 times) | 14.7 µs | 15.5 µs | 28.1 µs |
| 96 SIR candidates, 1 worker, per candidate | 11.2 µs | 11.4 µs | 21.7 µs |
| 48 Moran cases (N = 18, 65 updates), per case | 47.5 µs | 47.9 µs | 51.5 µs |
| 48 calendar cases (350 segments), per case | 0.02 µs | 0.27 µs | 0.81 µs |

- **Use `score_matrix` in optimization loops.** It skips candidate and
  record JSON in both directions. From Python it also skips `json.dumps`
  and `json.loads`, which cost more than the SIR evaluation itself for
  small models.
- **Full records cost about twice as much for SIR.** Most of the extra
  time goes to building and hashing the CRN-craft provenance documents and
  serializing the diagnostics. Request them for replay and inspection,
  not inside a search.
- **Traces are larger still.** A traced SIR record took 46.5 µs with the
  prepared projection.
- **Avoid monitored budgets in tight SIR loops.** A `wall_ms` or
  `cancel_file` budget routes SIR through the CRN-craft solver: 154.6 µs
  instead of 16.1 µs for one lean candidate. Moran and calendar are not
  affected this way.
- **Workers.** On this 2-vCPU host the 96-candidate SIR batch took
  2.75 ms with one worker, 1.71 ms with two and 1.90 ms with four (full
  records). Use at most one worker per available CPU, batch at least a
  few dozen candidates per call, and keep `workers=1` for single
  candidates, which always run inline.
- **Large invalid batches are JSON bound.** Decoding 1,024 candidates with
  63 numbers each took 1.14 ms of a 1.46 ms call. The matrix path scored
  the same rows in 0.02 ms.

## 13. Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| `ValueError: invalid SIR population, time grid, observations or tolerance` | check N ≤ 30, 2–64 increasing times in [0, 120], one observation per time within [0, N], `train_count` between 1 and the number of times minus 1 with a training observation, and tolerance in [1e-13, 1e-5] |
| every SIR row is `infeasible` | the first sample time must be 0; the `detail` string holds CRN-craft's diagnostics, including the `CRN-E020` error for the save grid |
| `work_limit` on a large Moran model | steps are (N+1)(N+2)/2 × (updates + 1); raise `max_steps` (up to 100,000,000) |
| `invalid` with `expected beta and gamma per time unit` | SIR candidates need exactly two parameters |
| `invalid` with `neutral Moran has no candidate parameters` or `calendar rates come from the model` | pass `[]`; change the model to change rates |
| `ValueError: Out of range float values are not JSON compliant` | a NaN or infinity in `batch`, `scores` or `evaluate`; filter it or use `score_matrix`, which classifies it as `invalid` |
| `ValueError: ragged parameter rows` | every `score_matrix` row needs the same length |
| `ValueError: at most 16 traced candidates per batch` | split traced replays into several batches |
| a campaign re-evaluates after reinstalling | the build hash changed; that is intended |
| the wheel does not install | the wheel targets Linux x86-64 and CPython 3.10+; build from the source distribution |
| timings vary between runs | use medians of repeated measurements and keep workers ≤ CPUs; shared hosts are noisy |

## 14. Example programs

[`examples/python/`](../examples/python/README.md) holds application and
optimization scripts. Their recorded outputs are in
`examples/python/outputs/`.

| Script | What it shows |
| --- | --- |
| `sir_rate_landscape.py` | a 41 x 41 rate grid in one `score_matrix` call, the R0 valley and a full replay |
| `sir_design_comparison.py` | identifiability of beta, gamma and R0 under three sampling designs, sensitivities and tolerance diagnostics |
| `moran_fixation_study.py` | neutral fixation for populations 4 to 32, the joint law and step budgets |
| `calendar_schedule_comparison.py` | on-rate schedules with equal exposure, a closed-form check and per-segment sensitivity |
| `optimize_scipy_de.py` | SciPy vectorized differential evolution through `score_matrix`, Nelder-Mead polish and traced replay |
| `optimize_fcmaes_cmaes.py` | fcmaes CMA-ES ask/tell with one `score_matrix` call per generation on a noisy 20-person outbreak |
| `optimize_fcmaes_minimize.py` | fcmaes `bitecpp.minimize` and `cmaescpp.minimize` with a graded failure penalty |
