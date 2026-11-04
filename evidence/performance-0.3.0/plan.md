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

## Round 2: prepared SIR projection

- Change: for SIR, build the reachable states, the generator edges and
  their column order once in `Prepared::new` (a new `src/sir_fsp.rs`). A
  candidate binds its two rates to the prepared structure and runs the
  same uniformization arithmetic as `crn_craft::fsp::solve`, in the same
  order, so means, diagnostics and the FSP work count are bit-identical.
  Traced candidates, monitored budgets (deadline or cancel file),
  non-normal rates and every runtime limit fall back to the vendored
  solver, which then produces the identical result or diagnostic. Model
  inputs that CRN-craft rejects (for example a save grid not starting at
  0) disable the prepared path at preparation. Full records still build
  the CRN model and run documents to report their SHA-256 values.
- Workload: `serial` through the matrix API (round 1 binary, matrix,
  against this round, matrix).
- Expected effect: per-candidate cost falls from about 140 us to about
  25 us, roughly 4x to 6x.
- Decision: the brief's rule (at least 1.03x and the paired check on all
  workloads). Results should be bit-identical; the independent domain
  validation runs for this round anyway because the FSP arithmetic moves
  into this crate.
- Result: 16.029 ms before, 2.212 ms after, 7.25x
  (`round02-matrix-serial.json`). Kept. `validation/validate.py` passed
  against the round 2 CLI and reproduced `models/example.json` and
  `evidence/validation.json` exactly.

## Profile after round 2 (round 2 binary, one run each)

| Workload | ms |
| --- | ---: |
| one (full) | 0.041 |
| serial (full) | 3.869 |
| lean serial / matrix serial | 2.199 / 2.149 |
| parallel (4 workers on 2 vCPUs) | 2.762 |
| trace (12 traced, vendored solver) | 2.179 |
| invalid | 1.277 |
| domain (Moran, full) | 4.427 |
| lean domain / matrix domain | 4.570 / 4.344 |
| domain2 (calendar, full) | 2.597 |

Moran (`domain`) is now the largest workload: 190 states, about 1,500
transitions and 65 updates, about 92 us per case.

## Round 3: gather the Moran update by target state

- Change: store Moran transitions by target state (incoming entries in
  ascending source order, then the source's parent/death order) and
  compute each next probability as one sum, instead of zeroing the next
  law and scattering every source's contributions. The sequence of
  floating-point additions per target is unchanged; zero-mass sources add
  +0.0, which leaves a non-negative sum unchanged. One logical step per
  state row per update is still charged, interleaved with the rows.
- Workload: `domain` (48 Moran cases, full records), round 2 against
  this round.
- Expected effect: 1.2x to 1.5x from removing the fill and the
  read-modify-write scatter.
- Decision: the brief's rule. Results should be bit-identical; a unit
  test compares the gather kernel with the scatter kernel bit for bit.
- Result: 4.432 ms before, 3.394 ms after, 1.31x (`round03-domain.json`).
  Kept. The unit test confirms bit-identical laws; `validate.py` also
  passed.
