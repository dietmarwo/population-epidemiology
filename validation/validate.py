from collections import Counter
import itertools
import sys

import numpy as np
from scipy.linalg import expm

from tools import good, run, save

binary=sys.argv[1]
sir={"kind":"sir","susceptible":9,"infected":2,"recovered":1,
     "sample_times":[0.,1.,2.,3.,5.],"observations":[2,3,3,None,2],"train_count":3,"tolerance":1e-11}
def oracle(model,beta,gamma):
    initial=tuple(model[k] for k in ("susceptible","infected","recovered"))
    n=sum(initial)
    space=[(s,i,n-s-i) for s in range(n+1) for i in range(n-s+1)]
    lookup={v:i for i,v in enumerate(space)}
    q=np.zeros((len(space),len(space)))
    for j,(s,i,r) in enumerate(space):
        if s and i:q[j,lookup[s-1,i+1,r]]=beta*s*i/n
        if i:q[j,lookup[s,i-1,r+1]]=gamma*i
        q[j,j]=-q[j].sum()
    assert np.max(abs(q.sum(axis=1)))<1e-14
    law=[expm(q*t)[lookup[initial]] for t in model["sample_times"]]
    return space,law
truth=[.7,.4]
space,law=oracle(sir,*truth)
means=[float(p@np.array([s[1] for s in space])) for p in law]
sir["observations"]=means
metrics=good(run(binary,sir,[truth],trace=True)[0])
fsp=metrics["fsp"]
lookup={tuple(s):i for i,s in enumerate(space)}
max_delta=0
for obs,ref in zip(fsp["snapshots"],law,strict=True):
    delta=max(abs(prob-ref[lookup[tuple(state)]]) for state,prob in zip(fsp["states"],obs["probabilities"],strict=True))
    max_delta=max(max_delta,delta)
    assert delta<1e-10 and obs["physical_leakage_estimate"]==0
    assert abs(sum(obs["probabilities"])-1)<1e-10
assert max_delta < 1e-10
params=[[b,g] for b in [.3,.5,.7,.9] for g in [.2,.4,.6]]
fit=run(binary,sir,params,budget={"workers":4})
best=min(fit,key=lambda r:r["score"])
assert params[int(best["id"])]==truth and good(best)["holdout_mse"]<1e-18
tight=good(run(binary,sir|{"tolerance":1e-13},[truth])[0])
assert tight["holdout_mse"]<1e-22
moran={"kind":"moran","population":3,"mutant_a":1,"mutant_b":1,"updates":3}
dist=Counter({(0,1,2):1})
for _ in range(moran["updates"]):
    nxt=Counter()
    for state,n in dist.items():
        for parent,dead in itertools.product(range(3),repeat=2):
            y=list(state)
            y[dead]=state[parent]
            nxt[tuple(y)]+=n
    dist=nxt
denom=9**moran["updates"]
m=good(run(binary,moran,[[]],trace=True)[0])
for cls,field in [(1,"marginal_a"),(2,"marginal_b")]:
    ref=[sum(n for state,n in dist.items() if state.count(cls)==i)/denom for i in range(4)]
    assert max(abs(a-b) for a,b in zip(ref,m[field],strict=True))<1e-14
long=good(run(binary,moran|{"updates":1200},[[]])[0])
assert abs(long["marginal_a"][-1]-1/3)<1e-12
calendar={"kind":"calendar","initial_on_probability":0,"segments":[{"duration":1,"on_rate":.3,"off_rate":0},
         {"duration":2,"on_rate":0,"off_rate":.4}]}
cal=good(run(binary,calendar,[[]])[0])
assert abs(cal["boundary_time_and_on_probability"][-1][1]-(1-np.exp(-.3))*np.exp(-.8))<1e-14
save({"model":sir,"parameters":truth},{"sir_state_count":len(fsp["states"]),"max_fsp_state_error":max_delta,
     "fsp_budget":metrics["error"],"fit_candidates":len(params),"winner":truth,
     "holdout_mse":good(best)["holdout_mse"],"tight_holdout_mse":tight["holdout_mse"],
     "moran_mass":m["mass"],"neutral_fixation":long["marginal_a"][-1],"calendar":cal})
