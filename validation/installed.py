"""Test an installed wheel from a directory outside its source tree."""

import importlib
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
from time import perf_counter

name = sys.argv[1]
package = importlib.import_module(name)
campaign = importlib.import_module(f"{name}.campaign")
example = package.load_example()
model, parameters = example["model"], example["parameters"]
prepared = package.Prepared(model)
one = prepared.evaluate(parameters, trace=True)
assert one["status"] == "complete", one

try:
    package.Prepared(model | {"unexpected": True})
except ValueError:
    pass
else:
    raise RuntimeError("unknown model field accepted")

cases = [{"id": f"candidate-{i}", "parameters": parameters} for i in range(12)]
serial = prepared.batch(cases, workers=1)
parallel = prepared.batch(list(reversed(cases)), workers=4)
assert serial["records"] == parallel["records"][::-1]
assert all(r["status"] == "complete" for r in serial["records"])
if name != "thermo_craft" or one["steps"] > 1:
    assert prepared.evaluate(parameters, max_steps=1)["status"] == "work_limit"
with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    flag = folder / "cancel"
    scenario = [{"id": "nominal", "model": model}]
    db = folder / "campaign.db"
    first = campaign.run(db, scenario, cases[:3], workers=2)
    flag.touch()
    interrupted = campaign.run(db, scenario, cases, workers=2, cancel_file=str(flag))
    assert interrupted[:3] == first
    assert {r["record"]["status"] for r in interrupted[3:]} == {"cancelled"}
    flag.unlink()
    resumed = campaign.run(db, scenario, cases, workers=4)
    again = campaign.run(db, scenario, cases, workers=1)
    assert resumed == again
    with sqlite3.connect(db) as conn:
        assert conn.execute("SELECT count(*) FROM results").fetchone()[0] == 12
        assert conn.execute("SELECT count(*) FROM attempts").fetchone()[0] == 21

def measure(fn, repeats=3):
    times = []
    for _ in range(repeats):
        start = perf_counter()
        fn()
        times.append(perf_counter()-start)
    return times

sample = [{"id":f"timed-{i}","parameters":parameters} for i in range(96)]
evidence = {
    "installed": True, "ordered_replay": True, "resumed_results": 12,
    "retained_attempts": 21, "introductory_score": one["score"],
    "batch_candidates":len(sample),
    "cold_prepare_and_run_seconds":measure(lambda:package.Prepared(model).batch(sample)),
    "prepared_worker_1_seconds":measure(lambda:prepared.batch(sample,workers=1)),
    "prepared_worker_4_seconds":measure(lambda:prepared.batch(sample,workers=4)),
}
print(json.dumps(evidence, indent=2))
