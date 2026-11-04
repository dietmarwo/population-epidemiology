# population-epidemiology performance implementation log

3 November 2026. Version 0.2.0, based on the prior 0.1.0 delivery
at `64515114da7d`. The ZIP contains the updated source, runnable
examples, an ordered patch series, [tutorial](docs/tutorial.md),
[architecture](docs/architecture.md), raw evidence and build artifacts.
The [earlier implementation record](docs/implementation-0.1.0.md)
documents the original research scope and remaining external work.

## Ten incremental changes

Each row compares the indicated stage with its immediate predecessor.
Times are median milliseconds per prepared call. A ratio below 1.00×
means the new median was slower for that measurement. Small differences
vary with thread scheduling, code layout and host noise.

| Round | Commit | Change | Case | Before ms | After ms | Before / after |
| ---: | --- | --- | --- | ---: | ---: | ---: |
| 1 | `982e547` | Evaluate singleton and one-worker batches on the calling thread, removing worker startup and join. | `one` | 0.370 | 0.250 | 1.48× |
| 2 | `bfb3cd6` | Replace worker-local row collections and final sorting with one result slot per candidate index. Keep dynamic work assignment. | `parallel` | 5.212 | 5.299 | 0.98× |
| 3 | `3cfd11b` | Validate unique candidate IDs through a pre-sized hash set. The large invalid batch isolates ID checks and response overhead from most domain computation. | `invalid` | 2.592 | 2.386 | 1.09× |
| 4 | `1fcc0de` | Serialize a typed prepared batch directly and keep the public CLI result adapter. | `trace` | 2.980 | 2.374 | 1.26× |
| 5 | `6d4dfbd` | Borrow candidate IDs until the batch is serialized, eliminating output-side string clones. | `invalid` | 1.738 | 1.558 | 1.12× |
| 6 | `6d4a7ef` | Inline the unmonitored work counter. Keep cancellation, deadline and exhausted-budget checks on the guarded path. | `domain` | 6.412 | 5.936 | 1.08× |
| 7 | `0bc81cf` | Add Rust and Python score-only entry points, preserving status, step counts, score and record order while skipping diagnostic construction. This row compares the earlier full path with the new lean path. | `lean-serial` | 14.073 | 14.700 | 0.96× |
| 8 | `4260532` | Build triangular Moran states and transition probabilities once during preparation, preserving logical budget charges. | `domain` | 5.942 | 4.375 | 1.36× |
| 9 | `0309099` | Replace per-state transition allocations with contiguous transition rows and offsets. Preserve update and accumulation order. | `domain` | 4.115 | 4.486 | 0.92× |
| 10 | `6ca965e` | Precompute the exact rate-calendar boundary law and bulk-charge unmonitored segment visits. Preserve per-segment budget checks under monitoring. | `lean-domain2` | 0.158 | 0.013 | 12.31× |

Each stage passed locked Rust release tests and paired output checks
for the supplied workloads. The result comparisons
ignore the package-version string and require matching record order,
status, detail, steps, provenance and numeric values within tolerance.
Lean comparisons also normalize metrics to null. The benchmark commit
precedes the ten rounds. Documentation and version changes follow them.

## Cumulative native results

| Workload | 0.1.0 ms | Final ms | Before / after |
| --- | ---: | ---: | ---: |
| Single prepared evaluation | 0.284 | 0.154 | 1.84× |
| 96 full records, one worker | 16.414 | 14.940 | 1.10× |
| 96 full records, four workers | 5.404 | 4.874 | 1.11× |
| 12 traced records | 3.016 | 2.267 | 1.33× |
| 1,024 invalid candidates | 2.570 | 1.549 | 1.66× |
| 48 domain cases | 7.394 | 4.811 | 1.54× |
| 48 secondary domain cases | 4.369 | 3.097 | 1.41× |
| 96 lean scores | 17.092 | 14.781 | 1.16× |
| 48 lean domain scores | 7.550 | 4.538 | 1.66× |
| 48 lean secondary scores | 4.453 | 0.013 | 344.38× |

The lean API was introduced in 0.2.0. The baseline uses full
records for the three lean rows. The paired comparisons alternate process order,
exclude startup and model preparation, warm up each prepared evaluator,
and measure JSON decode, evaluation, checks and JSON encoding. They
use five observations per side on an AMD EPYC 9V74 with Rust 1.88.0,
using logical CPUs 0–3. These are fixed test workloads, with four
workers only in the parallel row. Ratios from individual stages do
not multiply into the cumulative ratio because the timings were rerun.

The supplemental workloads use 48 neutral Moran cases at population 18 and 65 updates, and 48 fixed rate-calendar cases at 350 segments. SIR still runs the vendored CRN-craft FSP path per candidate, with a fresh projection per rate binding. The standalone registry crate remains unavailable because its CRN dependency is vendored by path. Source and Python builds use the pinned CRN core in the ZIP.

## Installed Python optimization

Seven runs per configuration used the same seed, bounds, starting point,
optimizer budget and final traced replay. CPython 3.12 and SciPy 1.17.0
ran on one logical CPU. Times include preparation, objective evaluation
and traced winner replay, excluding Python imports and wheel compilation.
The total median is calculated per run, separately from the column medians.

| Configuration | Prepare ms | Search ms | Replay ms | Total ms | Baseline / total |
| --- | ---: | ---: | ---: | ---: | ---: |
| v0.1.0-full | 0.070 | 23.000 | 0.460 | 23.772 | 1.00× |
| v0.2.0-full | 0.072 | 19.032 | 0.365 | 19.529 | 1.22× |
| v0.2.0-scores | 0.073 | 17.908 | 0.372 | 18.469 | 1.29× |

The run made **48** objective calls and classified
**0** failures. The selected vector was
`[0.7, 0.4]` with replay score **1.99962631623e-24**.
All configurations agreed on evaluation and failure counts, final
parameters and full replay score within the recorded tolerance.

## Release checks

- Rust formatting, Clippy with Python enabled and warnings denied,
  release tests, CLI request and native sweep passed.
- The independent domain validation passed against the release CLI.
- The wheel built from the extracted source distribution and installed
  into a fresh CPython 3.12 environment. Its ABI is CPython 3.10+ and
  its bundled target is Linux x86-64.
- The installed lean and full paths agreed on results, including mixed
  valid/invalid candidates, worker count, fixed work limits and
  cancellation. The tutorial's 4 Python blocks ran.
- SQLite cancellation/resume retained 21
  attempts and recovered 12 completed rows.
  Two full manifest runs against the same database produced identical output.
- The source ZIP and sdist are the Rust source delivery for the vendored population companion.

Run the packaged `validation/performance_pair.py` with saved native
benchmark binaries to reproduce comparisons. To reconstruct the ten
rounds in Git, apply `patches/*.patch` in filename order to the source
from the 0.1.0 ZIP. A baseline may first be committed locally. The patch
series includes the benchmark and release documentation commits around
the ten implementation commits.
