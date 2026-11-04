"""Fit beta and gamma with SciPy differential evolution over whole populations.

SciPy's ``vectorized=True`` mode hands the objective every trial vector of
a generation at once, as a (dimension, population) array. One
``score_matrix`` call scores the generation without per-candidate JSON.
Non-complete rows (here, rates outside the supported [0, 5] range because
the search box deliberately reaches 5.5) get a penalty. A Nelder-Mead
polish and a traced full replay of the winner follow.
"""

from time import perf_counter

import numpy as np
from scipy.optimize import differential_evolution, minimize

from population_epidemiology import Prepared, load_example

PENALTY = 1e6
case = load_example()
prepared = Prepared(case["model"])
calls = {"generations": 0, "rows": 0, "failures": 0}


def objective(x):
    rows = np.atleast_2d(x.T)  # one row per trial vector
    scores, statuses, _ = prepared.score_matrix(rows)
    calls["generations"] += 1
    calls["rows"] += len(rows)
    failed = [s != "complete" for s in statuses]
    calls["failures"] += sum(failed)
    return np.where(failed, PENALTY, scores)


start = perf_counter()
result = differential_evolution(objective, [(0.0, 5.5), (0.0, 5.5)], seed=11, popsize=15,
                                maxiter=60, tol=1e-12, polish=False, vectorized=True,
                                updating="deferred")
search = perf_counter() - start
print(f"DE: {calls['rows']} rows scored in {calls['generations']} score_matrix calls, "
      f"{calls['failures']} penalized rows, {search * 1e3:.0f} ms")
print(f"DE best {np.round(result.x, 6).tolist()} train MSE {result.fun:.3e}")

polish = minimize(lambda x: float(objective(x.reshape(2, 1))[0]), result.x, method="Nelder-Mead",
                  options={"xatol": 1e-9, "fatol": 1e-28, "maxiter": 400})
print(f"Nelder-Mead polish {np.round(polish.x, 9).tolist()} train MSE {polish.fun:.3e} "
      f"({polish.nfev} evaluations)")

row = prepared.evaluate(polish.x.tolist(), trace=True)
assert row["status"] == "complete", row
metrics = row["metrics"]
print("traced replay: holdout MSE", f"{metrics['holdout_mse']:.3e}",
      "| FSP states", len(metrics["fsp"]["states"]), "| FSP work", metrics["fsp_work"])
print("truth used to generate the packaged observations:", case["parameters"])
assert np.allclose(polish.x, case["parameters"], atol=1e-4)
