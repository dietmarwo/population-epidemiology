"""Prepared global search with SciPy DE. Optional demonstration, no runtime dependency."""

import importlib
import json
import sys

from scipy.optimize import differential_evolution

package = importlib.import_module(sys.argv[1])
ex = package.load_example()
prepared = package.Prepared(ex["model"])
x0 = ex["parameters"]
if not x0:
    raise SystemExit("this fixture has no search parameters")
bounds = [(max(0.0, v*.6), max(.05, v*1.4)) for v in x0]
failed = 0
def objective(params):
    global failed
    result = prepared.evaluate(params.tolist())
    if result["status"] != "complete":
        failed += 1
        return 1e20
    return result["score"]
out = differential_evolution(objective,bounds,seed=17,maxiter=5,popsize=4,polish=False,x0=x0,tol=0)
final = prepared.evaluate(out.x.tolist(),trace=True)
if final["status"] != "complete":
    raise RuntimeError(final)
print(json.dumps({"evaluations":int(out.nfev),"infeasible":failed,"best":out.x.tolist(),"score":final["score"]}))
