# population-epidemiology 0.3.0 performance plan

Host: 2 vCPUs, Intel(R) Xeon(R) Processor @ 2.10GHz (KVM), Rust 1.97.0,
CPython 3.12. Every timing runs inside `/home/claude/craft5/cpu.lock`.
Each entry below is written before the round is timed. The decision rule
is the brief's: keep a round when the declared workload's median speedup
under `validation/performance_pair.py --samples 11` against the previous
kept binary is at least 1.03x and the paired output check passes.
Otherwise revert it and record it under `rejected/`.

## Profile before round 1 (0.2.0 code on this host)

Single SIR candidate (example fixture, 12 individuals, 5 sample times),
phase times from a scratch harness that calls the vendored CRN-craft API:

| Phase | us |
| --- | ---: |
| Build model JSON (`json!` + `to_string`) | 5.6 |
| Build run JSON incl. model SHA-256 | 4.3 |
| `check_bytes` (hash, parse, validate model and run) | 32.7 |
| `parse_model` + `parse_run` (again inside `solve`) | 25.8 |
| `CompiledNetwork::compile` | 8.6 |
| Output-length serialization of `FspResult` | 13.6 |
| `fsp::solve` total | 159.3 |
| `scores_json`, one candidate | 172.7 |
| `evaluate_json`, one candidate | 202.3 |
| `scores_json`, one invalid candidate | 0.8 |

The projection has 75 states and 119 edges and charges 25,936 FSP work
units, so the uniformization arithmetic is roughly 25 us of the 159 us.
The rest is per-candidate JSON construction, validation, parsing,
compilation, BTreeMap state discovery and result-size serialization.
The JSON API boundary is about 1 us per candidate.

## Round 1: typed score-matrix API (required by the brief)

- Change: `Prepared::score_matrix`, `MatrixScore`, a Python native
  `score_matrix` that releases the GIL, the Python wrapper, a `matrix`
  API mode in `examples/performance.rs` (`CRAFT_MATRIX=1`), matrix flags in
  `validation/performance_pair.py` and `--matrix` in
  `validation/optimize_performance.py`.
- Workload (fixed by the brief): `serial`, lean `scores_json` before,
  matrix after.
- Expected effect: about 1.00x to 1.01x. Serial is 96 SIR candidates at
  about 170 us each, and the matrix API only removes the roughly 1 us per
  candidate JSON boundary.
- Decision: the brief's rule. Because the brief requires this API as
  round 1 and later steps depend on it, a ratio below 1.03x does not
  revert it. In that case the round is recorded as a mandated round below
  the threshold, it does not count toward the ten threshold rounds, and
  one extra round follows so that ten kept rounds meet the threshold.
- Result: lean 14.611 ms before, matrix 13.288 ms after, 1.10x
  (`round01-matrix-serial.json`). Kept under the brief's rule, so no extra
  round is needed. The gain is larger than the ~1 us per candidate JSON
  boundary seen in the single-candidate profile. An alternating
  single-candidate check showed equal per-call cost, so the batch gain
  comes from not decoding, holding and encoding 96 candidate records
  around the allocation-heavy SIR evaluations; host noise in diagnostic
  reruns was about +/-5%.
