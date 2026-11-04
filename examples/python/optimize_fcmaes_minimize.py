"""Fit the noisy 20-person outbreak with fcmaes standard minimizers and penalties.

``bitecpp.minimize`` (BiteOpt) and ``cmaescpp.minimize`` call a scalar
objective once per candidate. The objective scores one row with
``score_matrix`` and maps every non-complete status to a fixed penalty, so
the search box can extend past the supported rate range. Both runs use
the same budget and seed and are replayed with full records.

A flat penalty gives the optimizer no direction back into the supported
range, so the penalty grows with the distance outside [0, 5] per rate.
"""

from time import perf_counter

import numpy as np
from loguru import logger
from scipy.optimize import Bounds

from fcmaes import bitecpp, cmaescpp
from population_epidemiology import Prepared

logger.remove()
truth = [1.1, 0.35]
times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0]
base = {"kind": "sir", "susceptible": 16, "infected": 3, "recovered": 1,
        "sample_times": times, "train_count": 7, "tolerance": 1e-11}
means = Prepared(base | {"observations": [0.0] * len(times)}).evaluate(truth)["metrics"]["infected_means"]
noise = np.random.default_rng(5).normal(0.0, 0.15, len(times))
observations = np.clip(np.array(means) + noise, 0.0, 20.0)
observations[0] = means[0]
prepared = Prepared(base | {"observations": observations.tolist()})

PENALTY = 1e3
LIMIT = 5.0  # rates must be in [0, 5] per time unit
counts = {}


def objective(x):
    scores, statuses, _ = prepared.score_matrix([x])
    counts["calls"] = counts.get("calls", 0) + 1
    if statuses[0] != "complete":
        counts["failures"] = counts.get("failures", 0) + 1
        excess = np.sum(np.maximum(np.asarray(x) - LIMIT, 0.0) + np.maximum(-np.asarray(x), 0.0))
        return PENALTY * (1.0 + excess)
    return scores[0]


bounds = Bounds([0.0, 0.0], [6.0, 6.0])
for name, run in (
    ("bitecpp.minimize", lambda: bitecpp.minimize(objective, bounds, x0=[2.0, 2.0],
                                                  max_evaluations=1500,
                                                  rg=np.random.default_rng(7))),
    ("cmaescpp.minimize", lambda: cmaescpp.minimize(objective, bounds, x0=[2.0, 2.0],
                                                    input_sigma=0.3, popsize=16,
                                                    max_evaluations=1500,
                                                    rg=np.random.default_rng(7))),
):
    counts.clear()
    start = perf_counter()
    result = run()
    elapsed = perf_counter() - start
    row = prepared.evaluate(np.asarray(result.x).tolist())
    print(f"{name}: {result.nfev} evaluations, {counts.get('failures', 0)} penalized, "
          f"{elapsed * 1e3:.0f} ms")
    print(f"  beta {result.x[0]:.4f} gamma {result.x[1]:.4f} train MSE {row['score']:.4f} "
          f"holdout MSE {row['metrics']['holdout_mse']:.4f}")
    assert row["status"] == "complete" and abs(row["score"] - result.fun) < 1e-12
print("truth:", truth)
