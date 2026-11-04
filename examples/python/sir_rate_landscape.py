"""Which infection and recovery rates explain the packaged outbreak, and how sharply?

Scans a 41 x 41 grid of (beta, gamma) through one ``score_matrix`` call,
reports the best grid points with their basic reproduction number
R0 = beta / gamma, prints a coarse log10-loss map and replays the winner
with a full record to read its holdout error and fitted infected means.
"""

from time import perf_counter

import numpy as np

from population_epidemiology import Prepared, load_example

case = load_example()
model = case["model"]
prepared = Prepared(model)
print("model:", {k: model[k] for k in ("susceptible", "infected", "recovered", "train_count")})
print("sample times:", model["sample_times"])

betas = np.linspace(0.3, 1.2, 41)
gammas = np.linspace(0.1, 0.8, 41)
grid = np.array([(b, g) for b in betas for g in gammas])
start = perf_counter()
scores, statuses, steps = prepared.score_matrix(grid)
elapsed = perf_counter() - start
scores = np.array(scores)
assert set(statuses) == {"complete"}, set(statuses)
print(f"evaluated {len(grid)} rate pairs in {elapsed * 1e3:.1f} ms "
      f"({elapsed / len(grid) * 1e6:.1f} us per candidate, {steps[0]} logical steps each)")

order = np.argsort(scores)
print("\nbest grid points (training MSE of the expected infected count):")
print("  beta    gamma   R0      train MSE")
for i in order[:8]:
    b, g = grid[i]
    print(f"  {b:.4f}  {g:.4f}  {b / g:.3f}  {scores[i]:.3e}")

# Coarse map: log10 training loss on every fourth grid line.
print("\nlog10(train MSE), rows beta (top 0.30 to bottom 1.20), columns gamma 0.10 to 0.80:")
loss = np.log10(np.maximum(scores.reshape(len(betas), len(gammas)), 1e-30))
for r in range(0, len(betas), 4):
    print(f"  beta {betas[r]:.2f} " + " ".join(f"{loss[r, c]:6.1f}" for c in range(0, len(gammas), 4)))

best = grid[order[0]]
row = prepared.evaluate(best.tolist())
metrics = row["metrics"]
print("\nfull replay of the best grid point", best.tolist())
print("  status", row["status"], "score", f"{row['score']:.3e}", "holdout MSE", f"{metrics['holdout_mse']:.3e}")
print("  fitted infected means:", [round(m, 4) for m in metrics["infected_means"]])
print("  observations:         ", [None if o is None else round(o, 4) for o in model["observations"]])
print("  largest series tail bound:", max(e["series_tail_upper_bound"] for e in metrics["error"]))
# Near-optimal pairs share R0: the valley runs along beta / gamma.
near = order[: max(5, int(0.01 * len(grid)))]
r0 = grid[near, 0] / grid[near, 1]
print(f"\nR0 across the best 1% of the grid: {r0.min():.3f} to {r0.max():.3f}")
