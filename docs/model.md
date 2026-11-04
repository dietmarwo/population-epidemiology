# Population models and numerical contracts

The SIR pilot uses fixed population N and stochastic counts. Infection
has intensity beta × S × I / N and recovery has intensity gamma × I.
Both rates are per input time unit. CRN-craft FSP enumerates reachable
count states and propagates the finite distribution. Bounds cover the
entire fixed population. Its result exposes physical leakage, Poisson
series bounds and the roundoff mass residual. It preserves unnormalized
probabilities. Trace mode includes the original FSP result.

The fitting target is the model's expected infected count at declared
sample times. Training uses mean squared error over unmasked observations.
Holdout remains separate. Synthetic observations come from a dense
matrix exponential built independently of CRN-craft. A checked SSA
comparison uses a predeclared Hoeffding bound. This study uses no
empirical observation model.

The public Moran fixture uses three labels including two mutant types.
One discrete attempt selects a uniform parent and a uniform replacement
site, allowing self-replacement and null attempts. A finite recurrence
returns both mutant marginals and optionally the joint law. This is a
transparent neutral fixture separate from the unavailable research workload.

The calendar model is a two-state continuous-time chain with an exact
piecewise-constant propagation formula. At a boundary the preceding
segment has finished and the next rates take effect. Boundary probability
is continuous. It supplies a small schedule oracle and the SIR model
uses constant rates.
