"""Optional comparison with the CRN-craft Python wheel from the source baseline."""

import json
from pathlib import Path

import numpy as np
from crn_craft import simulate_ssa, ssa_spec

root=Path(__file__).resolve().parents[1]
case=json.loads((root/"models/example.json").read_text())
m=case["model"]
beta,gamma=case["parameters"]
n=sum(m[k] for k in ("susceptible","infected","recovered"))
model={"schema":"crn-craft.model/0.1","name":"small_sir","convention":"stochastic_counts","omega":1.0,
    "parameters":{"infection":beta/n,"gamma":gamma},
    "species":[{"name":name,"initial_count":m[key]} for name,key in
               zip("SIR",["susceptible","infected","recovered"],strict=True)],
    "reactions":[{"name":"infection","reactants":{"S":1,"I":1},"products":{"I":2},
                 "rate":{"kind":"mass_action","coefficient":"infection"}},
                {"name":"recovery","reactants":{"I":1},"products":{"R":1},
                 "rate":{"kind":"mass_action","coefficient":"gamma"}}]}
r=simulate_ssa(model,ssa_spec(model,aggregator="direct",fast=False,replicates=4096,threads=4,
                              seed=4242,end=m["sample_times"][-1],save_times=m["sample_times"]))
assert r["termination"]=="completed",r["termination"]
mean=np.asarray(r["statistics"]["mean"])[:,1]
ref=np.array(m["observations"])
# Bounded-variable Hoeffding union bound for all saved infected-count means.
bound=n*np.sqrt(np.log(2*len(ref)/.001)/(2*4096))
delta=float(np.max(abs(mean-ref)))
assert delta<bound
(root/"evidence/crn-ssa.json").write_text(json.dumps({"replicates":4096,"seed":4242,
     "alpha":.001,"bound":float(bound),"max_mean_difference":delta},indent=2)+"\n")
