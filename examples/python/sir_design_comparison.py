"""How well do different sampling designs pin down beta and gamma?

Builds synthetic observations from known rates for three sampling designs,
then measures how much the training loss rises when each rate, or both
along the R0 = beta / gamma direction, moves by 10%. A design that barely
changes the loss leaves that direction poorly identified. It also reports
local sensitivities of the expected infected count and the FSP numerical
diagnostics at two tolerances.
"""

import numpy as np

from population_epidemiology import Prepared

truth = np.array([0.9, 0.3])
population = {"susceptible": 14, "infected": 2, "recovered": 0}
designs = {
    "early (0-2)": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
    "spread (0-10)": [0.0, 2.0, 4.0, 6.0, 8.0, 10.0],
    "late sparse": [0.0, 6.0, 9.0, 12.0, 15.0, 18.0],
}


def sir(times, observations, tolerance=1e-11):
    return {"kind": "sir", **population, "sample_times": times,
            "observations": observations, "train_count": 4, "tolerance": tolerance}


for name, times in designs.items():
    # Synthetic data: the model's own expected infected counts at the truth.
    reference = Prepared(sir(times, [0.0] * len(times)))
    means = reference.evaluate(truth.tolist())["metrics"]["infected_means"]
    prepared = Prepared(sir(times, means))
    moves = {
        "beta +10%": truth * [1.1, 1.0],
        "gamma +10%": truth * [1.0, 1.1],
        "both +10% (same R0)": truth * 1.1,
        "beta +10%, gamma -10%": truth * [1.1, 0.9],
    }
    scores, statuses, _ = prepared.score_matrix([truth, *moves.values()])
    assert set(statuses) == {"complete"}
    print(f"\n{name}: times {times}")
    print(f"  loss at truth {scores[0]:.2e}; peak expected infected {max(means):.2f}")
    for label, value in zip(moves, scores[1:], strict=True):
        print(f"  {label:24s} train MSE {value:.3e}")
    # Local sensitivity of the expected infected count by central differences.
    h = 1e-4
    rows = [truth + [h, 0], truth - [h, 0], truth + [0, h], truth - [0, h]]
    batch = prepared.batch([{"id": str(i), "parameters": r.tolist()} for i, r in enumerate(rows)])
    m = [np.array(r["metrics"]["infected_means"]) for r in batch["records"]]
    print("  d mean / d beta :", np.round((m[0] - m[1]) / (2 * h), 3).tolist())
    print("  d mean / d gamma:", np.round((m[2] - m[3]) / (2 * h), 3).tolist())

# Numerical diagnostics at a loose and a tight supported tolerance.
times = designs["spread (0-10)"]
for tolerance in (1e-5, 1e-13):
    row = Prepared(sir(times, [0.0] * len(times), tolerance)).evaluate(truth.tolist())
    tail = max(e["series_tail_upper_bound"] for e in row["metrics"]["error"])
    print(f"\ntolerance {tolerance:g}: largest series tail bound {tail:.2e}, "
          f"FSP work {row['metrics']['fsp_work']}")
