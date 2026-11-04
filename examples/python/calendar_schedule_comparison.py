"""Which switching schedule leaves a two-state system most likely in the on state?

Compares constant, front-loaded, back-loaded and pulsed on-rate calendars
over ten time units with the same total on-rate exposure (rate x time),
reports the final and time-averaged on-probability, checks the constant
schedule against the closed form and measures how the final probability
responds to each quarter's on-rate.
"""

import math

from population_epidemiology import Prepared

OFF = 0.5
QUARTERS = 4
SPAN = 10.0


def calendar(on_rates, initial=0.0, pieces=25):
    """Split each quarter into equal segments so the boundary trace is smooth."""
    duration = SPAN / QUARTERS / pieces
    segments = [{"duration": duration, "on_rate": rate, "off_rate": OFF}
                for rate in on_rates for _ in range(pieces)]
    return {"kind": "calendar", "initial_on_probability": initial, "segments": segments}


def summary(model):
    row = Prepared(model).evaluate([])
    boundaries = row["metrics"]["boundary_time_and_on_probability"]
    # Trapezoid rule over the boundary trace, starting from the initial probability.
    points = [(0.0, model["initial_on_probability"])] + [tuple(b) for b in boundaries]
    area = sum((t1 - t0) * (p0 + p1) / 2 for (t0, p0), (t1, p1) in zip(points, points[1:]))
    return row["score"], area / SPAN, row["steps"]


schedules = {
    "constant": [1.0, 1.0, 1.0, 1.0],
    "front-loaded": [2.5, 1.0, 0.5, 0.0],
    "back-loaded": [0.0, 0.5, 1.0, 2.5],
    "pulsed": [2.0, 0.0, 2.0, 0.0],
}
print("schedule       final P(on)  mean P(on)  steps   (on-rate x time = 10 for each)")
for name, rates in schedules.items():
    assert abs(sum(rates) * SPAN / QUARTERS - 10.0) < 1e-12
    final, mean, steps = summary(calendar(rates))
    print(f"{name:14s} {final:.6f}     {mean:.6f}    {steps}")

# Closed form for a constant schedule: p(t) = p_eq (1 - exp(-(on + off) t)).
final, _, _ = summary(calendar([1.0] * 4))
expected = 1.0 / (1.0 + OFF) * (1 - math.exp(-(1.0 + OFF) * SPAN))
print(f"\nconstant schedule: {final:.15f} vs closed form {expected:.15f}")
assert abs(final - expected) < 1e-12

# Sensitivity of the final probability to each quarter's on-rate.
base = schedules["pulsed"]
h = 1e-6
print("\nd final P(on) / d on-rate, pulsed schedule:")
for q in range(QUARTERS):
    up = list(base)
    down = list(base)
    up[q] += h
    down[q] = max(down[q] - h, 0.0)
    slope = (summary(calendar(up))[0] - summary(calendar(down))[0]) / (up[q] - down[q])
    print(f"  quarter {q + 1}: {slope:+.6f}")
