"""Check the installed lean API and run the Python snippets in the tutorial."""

import importlib
import json
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

module, tutorial = sys.argv[1:]
package = importlib.import_module(module)
example = package.load_example()
prepared = package.Prepared(example["model"])
candidates = [{"id": "ok", "parameters": example["parameters"]},
              {"id": "bad", "parameters": [1] * 63}]
with TemporaryDirectory() as tmp:
    flag = Path(tmp) / "cancel"
    flag.touch()
    for budget in ({}, {"workers": 4}, {"max_steps": 1}, {"cancel_file": str(flag)}):
        full = prepared.batch(candidates, **budget)
        lean = prepared.scores(candidates, **budget)
        assert lean["model_sha256"] == full["model_sha256"]
        assert lean["records"] == [r | {"metrics": None} for r in full["records"]]
        if "cancel_file" in budget:
            assert lean["records"][0]["status"] == "cancelled"
    if module == "population_epidemiology":
        calendar = package.Prepared({"kind": "calendar", "initial_on_probability": 0.5,
                                     "segments": [{"duration": 1, "on_rate": 0, "off_rate": 0}] * 3})
        query = [{"id": "calendar", "parameters": []}]
        limited = calendar.scores(query, max_steps=2)["records"][0]
        assert (limited["status"], limited["steps"]) == ("work_limit", 2)
        canceled = calendar.scores(query, cancel_file=str(flag))["records"][0]
        assert (canceled["status"], canceled["steps"]) == ("cancelled", 0)
        complete = calendar.scores(query, max_steps=3)["records"][0]
        assert (complete["status"], complete["steps"], complete["score"]) == ("complete", 3, 0.5)
    if module == "rail_craft":
        closed = package.Prepared({
            "resources": [{"name": "one", "clearance": 0}], "horizon": 20,
            "trains": [{"name": "train", "release": 0, "planned": 0, "priority": 0,
                        "route": [{"resource": 0, "duration": 1}]}],
            "closures": [{"resource": 0, "start": a, "end": b} for a, b in ((0, 4), (2, 7), (7, 8))],
        })
        row = closed.evaluate([0, 0], trace=True)
        assert row["status"] == "complete"
        assert row["metrics"]["intervals"][0]["entry"] == 8
    if module == "robotics_control_craft":
        long = package.Prepared(example["model"] | {"dt_s": 1e-4})
        assert long.evaluate(example["parameters"])["status"] == "complete"
        traced = long.evaluate(example["parameters"], trace=True)
        assert traced["status"] == "invalid" and traced["steps"] == 10_001
try:
    prepared.scores([candidates[0] | {"trace": True}])
except ValueError:
    pass
else:
    raise AssertionError("scores accepted a traced candidate")

snippets = re.findall(r"^```python\n(.*?)^```", Path(tutorial).read_text(), re.M | re.S)
assert len(snippets) >= 4
namespace = {}
for i, snippet in enumerate(snippets, 1):
    exec(compile(snippet, f"tutorial-{i}", "exec"), namespace)
print(json.dumps({"status": "passed", "tutorial_python_blocks": len(snippets),
                  "lean_api_and_budget_parity": True}))
