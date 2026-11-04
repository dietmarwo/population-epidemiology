"""Fit a larger, noisy outbreak with fcmaes CMA-ES through its ask/tell interface.

The data are 11 noisy readings of the expected infected count of a
20-person SIR outbreak (true beta 1.1, gamma 0.35). The first 7 points
train the fit and the last 4 are held out. Each CMA-ES generation is
scored with one ``score_matrix`` call on two native workers. The search
box reaches past the supported rate range, so any out-of-range row comes
back ``invalid`` and receives a penalty instead of a score.
"""

from time import perf_counter

import numpy as np
from loguru import logger
from scipy.optimize import Bounds

from fcmaes import cmaes
from population_epidemiology import Prepared

logger.remove()
truth = [1.1, 0.35]
times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0]
base = {"kind": "sir", "susceptible": 16, "infected": 3, "recovered": 1,
        "sample_times": times, "train_count": 7, "tolerance": 1e-11}
means = Prepared(base | {"observations": [0.0] * len(times)}).evaluate(truth)["metrics"]["infected_means"]
noise = np.random.default_rng(5).normal(0.0, 0.15, len(times))
observations = np.clip(np.array(means) + noise, 0.0, 20.0)
observations[0] = means[0]  # the initial count is known exactly
prepared = Prepared(base | {"observations": observations.tolist()})

PENALTY = 1e6
rng = np.random.default_rng(23)
# With normalize=True (the default), input_sigma is relative to the box
# mapped onto [-1, 1]; a much larger value clips most samples to its faces.
optimizer = cmaes.Cmaes(bounds=Bounds([0.0, 0.0], [6.0, 6.0]), x0=np.array([2.0, 2.0]),
                        input_sigma=0.3, popsize=16, max_evaluations=1600, rg=rng,
                        randn=lambda *shape: rng.standard_normal(shape))
evaluations = failures = generations = 0
start = perf_counter()
while True:
    population = optimizer.ask()
    scores, statuses, _ = prepared.score_matrix(population, workers=2)
    values = np.array([s if st == "complete" else PENALTY for s, st in zip(scores, statuses)])
    failures += int(np.sum(values == PENALTY))
    evaluations += len(population)
    generations += 1
    if optimizer.tell(values) != 0 or evaluations >= 1600:
        break
elapsed = perf_counter() - start
result = optimizer.result()
x = np.asarray(result.x)
print(f"CMA-ES: {evaluations} evaluations in {generations} generations, "
      f"{failures} penalized rows, {elapsed * 1e3:.0f} ms")
print(f"estimate beta {x[0]:.4f} gamma {x[1]:.4f} (truth {truth}), R0 {x[0] / x[1]:.3f}")
row = prepared.evaluate(x.tolist())
print(f"train MSE {row['score']:.4f} holdout MSE {row['metrics']['holdout_mse']:.4f}")
truth_row = prepared.evaluate(truth)
print(f"at the truth: train MSE {truth_row['score']:.4f} holdout MSE "
      f"{truth_row['metrics']['holdout_mse']:.4f}")
print("observed:", np.round(observations, 3).tolist())
print("fitted:  ", np.round(row["metrics"]["infected_means"], 3).tolist())
assert row["score"] <= truth_row["score"] + 1e-9
