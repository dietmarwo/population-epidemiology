"""When does neutral drift fix a mutant, and how does that depend on population size?

For populations of 4 to 32 with a quarter of the individuals carrying each
of two neutral mutants, reports the probability that mutant A, mutant B or
the wild type has reached fixation after a given number of update attempts.
Neutral theory says each type eventually fixes with probability equal to
its initial frequency (0.25, 0.25 and 0.5 here); the long runs approach
that limit. The score is P(A fixed) + P(B fixed). The joint law from a
traced evaluation gives the wild-type fixation probability. Large runs
need an explicit ``max_steps`` because each update charges one logical
step per state.
"""

from population_epidemiology import Prepared

for n in (4, 8, 16, 32):
    a0 = b0 = n // 4
    states = (n + 1) * (n + 2) // 2
    print(f"\npopulation {n} ({states} states), initial A = B = {a0}")
    print("  updates   P(A fixed)  P(B fixed)  P(wild fixed)  P(segregating)  steps")
    for updates in (n, 4 * n, 16 * n, 64 * n, 256 * n):
        if updates > 10_000:
            continue
        model = {"kind": "moran", "population": n, "mutant_a": a0, "mutant_b": b0,
                 "updates": updates}
        row = Prepared(model).evaluate([], trace=True, max_steps=10_000_000)
        assert row["status"] == "complete", row
        m = row["metrics"]
        fixed_a, fixed_b = m["marginal_a"][n], m["marginal_b"][n]
        law = dict(zip(map(tuple, m["joint"]["states"]), m["joint"]["probabilities"]))
        fixed_wild = law[(0, 0)]
        segregating = 1 - fixed_a - fixed_b - fixed_wild
        print(f"  {updates:7d}   {fixed_a:.6f}    {fixed_b:.6f}    {fixed_wild:.6f}       "
              f"{segregating:.6f}        {row['steps']}")
        assert abs(row["score"] - (fixed_a + fixed_b)) < 1e-15
        assert abs(m["mass"] - 1) < 1e-12
        assert row["steps"] == states * (updates + 1)

# The default budget of 1,000,000 logical steps stops the largest case early.
model = {"kind": "moran", "population": 32, "mutant_a": 8, "mutant_b": 8, "updates": 8192}
row = Prepared(model).evaluate([])
print(f"\npopulation 32, 8192 updates, default budget: {row['status']} after {row['steps']} "
      f"steps ({row['detail']})")

# The joint law of a small population, as returned by a traced evaluation.
row = Prepared({"kind": "moran", "population": 3, "mutant_a": 1, "mutant_b": 1,
                "updates": 3}).evaluate([], trace=True)
joint = row["metrics"]["joint"]
print("\njoint law, population 3 after 3 updates (state (a, b): probability):")
for state, p in zip(joint["states"], joint["probabilities"], strict=True):
    if p > 0:
        print(f"  {tuple(state)}: {p:.6f}")
