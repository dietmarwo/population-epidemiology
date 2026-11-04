# population-epidemiology performance implementation log

4 November 2026. Version 0.3.0, based on the 0.2.0 delivery imported
as commit `9e427df`. The ZIP contains the updated source, runnable
Python and Rust examples with their outputs, an ordered patch series,
the new [user manual](docs/user-manual.md), the updated
[architecture guide](docs/architecture.md) and [tutorial](docs/tutorial.md),
raw evidence and build artifacts. The
[0.2.0 log](population-epidemiology-performance-log-0.2.0.md) is kept as
delivered. Its numbers come from another host and are not compared with
these.

All measurements below ran on **2 vCPUs, Intel(R) Xeon(R) Processor @
2.10GHz (KVM), Rust 1.97.0 and CPython 3.12**, each inside a shared CPU
lock that kept builds of other projects from overlapping the timings.

## Ten incremental changes

Each row compares the indicated round with the previous kept round's
benchmark binary. Times are median milliseconds per prepared call over
11 alternating samples of `validation/performance_pair.py`. The workload
and decision rule of every round were written to
`evidence/performance-0.3.0/plan.md` before it was timed. A round was
kept at 1.03x or more with a passing output check. A ratio below 1.00x
would mean the new median was slower.

| Round | Commit | Change | Case | Before ms | After ms | Before / after |
| ---: | --- | --- | --- | ---: | ---: | ---: |
| 1 | `61db6b1` | Add the typed score-matrix API: `Prepared::score_matrix` and `MatrixScore`, a GIL-free Python `score_matrix`, a shared dispatch path, a matrix benchmark mode and matrix flags for the pair and optimization scripts. Lean `scores_json` before, matrix after. | `serial` | 14.611 | 13.288 | 1.10x |
| 2 | `e2b042a` | Prepare the SIR state space, generator edges and column order once and bind the two rates per candidate, reproducing the vendored FSP arithmetic bit for bit. Traced, monitored and edge-case candidates fall back to the vendored solver. | matrix `serial` | 16.029 | 2.212 | 7.25x |
| 3 | `5a2ad07` | Gather each Moran target's incoming transitions in ascending source order instead of zeroing and scattering. | `domain` | 4.432 | 3.394 | 1.31x |
| 4 | `155cb29` | Charge unmonitored Moran row visits with one bulk step and keep per-row checks under monitoring. | lean `domain` | 3.292 | 3.051 | 1.08x |
| 5 | `5c471dc` | Carry metrics as serialized JSON and render the calendar's model-only metrics once at preparation. | `domain2` | 2.653 | 0.035 | 76.78x |
| 6 | `f820ce1` | Trace SIR candidates through the prepared projection and assemble the vendored solver's full result. | `trace` | 2.308 | 1.237 | 1.87x |
| 7 | `8a1a637` | Serialize SIR metrics from typed structs in sorted key order, without a `Value` tree. | `trace` | 1.090 | 0.817 | 1.33x |
| 8 | `900fce1` | Store SIR generator columns as three fixed slots, padded with zero weights, with the sink handled once. | matrix `serial` | 3.025 | 1.532 | 1.97x |
| 9 | `76cd01e` | Update Moran laws with a nine-slot stencil over a zero-padded triangular grid. | lean `domain` | 3.673 | 2.346 | 1.57x |
| 10 | `09246d9` | Run the SIR series as a three-stream stencil on the (S, I) grid, after checking the discovery order. | matrix `serial` | 1.841 | 1.597 | 1.15x |

Each round passed `cargo fmt --check`, Clippy with the `python` feature
and warnings denied, locked release tests and the paired output check on
every workload and API mode, comparing record order, status, detail,
steps, provenance, scores and metrics with the existing tolerances.
Lean comparisons normalize metrics to null, and matrix comparisons
normalize metrics and detail. Every SIR and Moran kernel change is
bit-identical by construction and covered by unit tests against the
vendored solver or the previous kernel. `validation/validate.py` also
passed against the release CLI for rounds 2, 3, 5, 6, 9 and 10, and the
traced CLI request printed the same bytes as the 0.2.0 CLI after rounds
5, 6, 7 and 10. The logical step counts did not change.

The host was shared, and absolute times drifted by up to about 20%
between periods. Two consecutive 7-sample diagnostic pairs of the same
round 4 and round 5 binaries gave 0.78x and 1.07x. That is why each "before"
column is re-measured against the previous kept binary instead of reusing
the last "after" value, and why ratios from individual rounds do not
multiply into the cumulative ratio.

The registry release of CRN-craft 0.2.0 was checked first: `cargo info
crn-craft` and the crates.io API report that the crate does not exist.
`vendor/crn-craft` stays as delivered, with no vendored patch. The
prepared SIR path is code in this crate that calls CRN-craft's public API.

## Rejected attempts

| After round | Change | Case | Before ms | After ms | Before / after | Evidence |
| ---: | --- | --- | ---: | ---: | ---: | --- |
| 4 | Template the CRN-craft model and run documents at preparation and assemble their bytes per full SIR record instead of building two `json!` trees | `serial` | 4.319 | 4.215 | 1.02x | `evidence/performance-0.3.0/rejected/attempt-document-template-serial.json` |

The attempt was reverted and does not count as a round. The later time
breakdown attributes about 12.5 µs per full SIR record to provenance and
metrics, so the attempt may have been hit by host noise. It is the
largest known remaining saving for full SIR records.

## Cumulative native results

`validation/performance_pair.py --samples 11`, the 0.2.0 benchmark binary
built from the import commit against round 10. The `parallel` workload
keeps its four workers on this 2-vCPU host. The two matrix rows use the
0.2.0 lean path as the baseline, because 0.2.0 has no matrix API.
Evidence: `evidence/performance-0.3.0/cumulative/`.

| Workload | Case | 0.2.0 ms | 0.3.0 ms | Before / after |
| --- | --- | ---: | ---: | ---: |
| Single prepared SIR evaluation | `one` | 0.191 | 0.031 | 6.09x |
| 96 full records, one worker | `serial` | 14.643 | 2.598 | 5.64x |
| 96 full records, four workers on 2 vCPUs | `parallel` | 9.466 | 1.943 | 4.87x |
| 12 traced records | `trace` | 2.289 | 0.593 | 3.86x |
| 1,024 invalid candidates | `invalid` | 1.640 | 1.936 | 0.85x |
| 48 Moran cases | `domain` | 4.088 | 2.468 | 1.66x |
| 48 calendar cases | `domain2` | 2.710 | 0.051 | 52.87x |
| 96 lean scores | lean `serial` | 14.088 | 1.456 | 9.68x |
| 48 lean Moran scores | lean `domain` | 4.710 | 3.160 | 1.49x |
| 48 lean calendar scores | lean `domain2` | 0.013 | 0.013 | 1.02x |
| 96 matrix scores (0.2.0 lean) | matrix `serial` | 14.266 | 1.291 | 11.05x |
| 48 matrix Moran scores (0.2.0 lean) | matrix `domain` | 4.454 | 2.468 | 1.80x |

No round targeted the `invalid` workload, which decodes 64,512 JSON
numbers per call. Two further 11-sample pairs of the same two binaries
gave 1.02x and 0.95x (`cumulative/invalid-repeat-1.json` and `-2.json`),
so it is unchanged within this host's noise. The lean calendar path was
already 0.013 ms in 0.2.0.

`examples/breakdown.rs` splits the 0.3.0 calls into decoding, evaluation
and encoding (`evidence/performance-0.3.0/breakdown.json`). Per candidate
in the 96-candidate SIR batch, the matrix path took 11.2 µs, lean 11.4 µs
and full 21.7 µs. One SIR candidate through the vendored solver, which
monitored budgets select, took 154.6 µs against 16.1 µs through the
prepared projection.

## Installed Python optimization

Seven interleaved runs per configuration used the same seed, bounds,
starting point, optimizer budget and final traced replay, on one logical
CPU (`taskset -c 0`) with CPython 3.12, NumPy 2.3.5 and SciPy 1.17.0. The
0.2.0 wheel came from the 0.2.0 `packages/`. The 0.3.0 wheel was built
with maturin from the extracted 0.3.0 source distribution. Its native
module is byte-identical to the released wheel's. Each configuration ran
in its own fresh venv. Times include preparation, objective evaluation
and traced winner replay, excluding Python imports and wheel
compilation. The total median is calculated per run, separately from the
column medians. Evidence: `evidence/performance-0.3.0/python-optimization.json`.

| Configuration | Prepare ms | Search ms | Replay ms | Total ms | Baseline / total |
| --- | ---: | ---: | ---: | ---: | ---: |
| v0.2.0-full | 0.092 | 17.068 | 0.479 | 17.694 | 1.00x |
| v0.2.0-scores | 0.091 | 15.481 | 0.441 | 16.011 | 1.11x |
| v0.3.0-full | 0.294 | 7.366 | 0.255 | 7.909 | 2.24x |
| v0.3.0-scores | 0.295 | 5.564 | 0.288 | 6.125 | 2.89x |
| v0.3.0-matrix | 0.350 | 5.639 | 0.329 | 6.367 | 2.78x |

Every run made **48** objective calls and classified **0** failures. The
selected vector was `[0.7, 0.4]` with replay score
**1.9996263162291143e-24**. All configurations agreed on evaluation and
failure counts, final parameters and replay score. SciPy's one-at-a-time
differential evolution spends about 117 µs per evaluation in 0.3.0, most
of it outside the native call, so the lean and matrix objectives are
equal within noise in this search. Preparation grew to about 0.3 ms
because it now validates the CRN-craft documents and builds the SIR
projection.

## Release checks

All steps ran inside the CPU lock; the summary is
`evidence/release-0.3.0/release-checks.json`.

- `cargo fmt --check`, Clippy with `--features python` and warnings
  denied, and the 13 locked release tests passed
  (`evidence/release-0.3.0/rust-checks.log`).
- The CLI request (`evidence/cli.json`) equals the 0.2.0 output except
  for the version string. The native sweep (`evidence/native-sweep.json`)
  equals the 0.2.0 output, and `describe` reports 0.3.0
  (`evidence/capabilities.json`).
- The independent domain validation passed against the release CLI. It
  regenerated `models/example.json` byte for byte, and
  `evidence/validation.json` equals the 0.2.0 record.
- The wheel was built from the extracted source distribution and
  installed into a fresh CPython 3.12 environment. Its ABI is CPython
  3.10+ and its target is Linux x86-64.
- `validation/installed.py` passed (`evidence/installed.json`): ordered
  replay with one and four workers, fixed work limits, and SQLite
  cancellation and resume, which retained 21 attempts and recovered 12
  completed rows.
- `validation/release_checks.py` passed (`evidence/tutorial-check.json`).
  The installed lean, full and matrix paths agreed on mixed valid and
  invalid candidates, worker count, work limits and cancellation; NumPy
  matrix input matched lists; the tutorial's 5 and the user manual's 18
  Python blocks ran.
- Two full manifest runs against the same database produced identical
  output (`evidence/campaign-replay.json`).
- Every example under `examples/python/` ran against the installed
  wheel; outputs are in `examples/python/outputs/`.
- `validation/check_crn_ssa.py` ran with an existing CRN-craft 0.2.0
  Python wheel and reproduced `evidence/crn-ssa.json` exactly.
- `cargo package` still fails, because `crn-craft` is a path dependency
  without a crates.io release (`evidence/release-0.3.0/cargo-package.txt`).
  As in 0.2.0, the source ZIP and sdist are the Rust source delivery, with
  the pinned CRN core vendored.

## Reproduction

The round evidence pairs benchmark binaries saved on the release host
as `roundNN/performance` for rounds 0 to 10; they are not in the ZIP. To
rebuild them, apply `patches/0.3.0/*.patch` in filename order to the
0.2.0 source, committing after each patch, and run
`cargo build --release --locked --example performance` at each round
commit. Then pair two binaries, for example:

```sh
python validation/performance_pair.py round09/performance round10/performance \
    models/example.json serial --samples 11 --before-matrix --after-matrix
python validation/performance_pair.py round00/performance round10/performance \
    models/example.json serial --samples 11 --before-scores --after-matrix
cargo run --release --locked --example breakdown -- models/example.json
taskset -c 0 python validation/optimize_performance.py population_epidemiology --matrix
```

`patches/0.3.0/` holds the ten rounds and the documentation and release
commits from the import commit, without the binary `packages/` files. The
0.2.0 series moved to `patches/0.2.0/` and applies to the 0.1.0 source.
`SHA256SUMS` lists every file in the ZIP.
