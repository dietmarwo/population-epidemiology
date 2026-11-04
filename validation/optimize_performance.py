"""Time one complete, seeded optimization plus a full solution replay."""

import argparse
import importlib
import json
from time import perf_counter

from scipy.optimize import differential_evolution


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("module")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--scores", action="store_true")
    modes.add_argument("--matrix", action="store_true")
    args = parser.parse_args()
    package = importlib.import_module(args.module)
    example = package.load_example()
    x0 = example["parameters"]
    bounds = [(max(0.0, x * .6), max(.05, x * 1.4)) for x in x0]
    failures = 0

    def objective(parameters):
        nonlocal failures
        if args.matrix:
            scores, statuses, _ = prepared.score_matrix([parameters.tolist()])
            if statuses[0] != "complete":
                failures += 1
                return 1e20
            return scores[0]
        if args.scores:
            row = prepared.scores([{"id": "one", "parameters": parameters.tolist()}])["records"][0]
        else:
            row = prepared.evaluate(parameters.tolist())
        if row["status"] != "complete":
            failures += 1
            return 1e20
        return row["score"]

    start = perf_counter()
    prepared = package.Prepared(example["model"])
    initialized = perf_counter()
    best = differential_evolution(
        objective, bounds, seed=17, maxiter=5, popsize=4, polish=False, tol=0, x0=x0,
    )
    searched = perf_counter()
    row = prepared.evaluate(best.x.tolist(), trace=True)
    finished = perf_counter()
    if row["status"] != "complete":
        raise RuntimeError(row)
    print(json.dumps({
        "mode": "matrix" if args.matrix else "scores" if args.scores else "full",
        "total_seconds": finished - start,
        "prepare_seconds": initialized - start,
        "search_seconds": searched - initialized,
        "replay_seconds": finished - searched,
        "evaluations": int(best.nfev),
        "failures": failures,
        "best": best.x.tolist(),
        "score": row["score"],
    }))


if __name__ == "__main__":
    main()
