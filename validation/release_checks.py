"""Check the installed lean and matrix APIs and run the Python blocks of the tutorial
and the user manual.

Usage: release_checks.py MODULE TUTORIAL [MANUAL]. MANUAL defaults to
user-manual.md next to TUTORIAL. Each document's blocks run in order in one
namespace.
"""

import importlib
import json
import math
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

module, tutorial = sys.argv[1:3]
manual = sys.argv[3] if len(sys.argv) > 3 else str(Path(tutorial).with_name("user-manual.md"))
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
        # The typed matrix API agrees with the lean records row by row.
        for candidate, row in zip(candidates, lean["records"], strict=True):
            scores, statuses, steps = prepared.score_matrix([candidate["parameters"]], **budget)
            assert (statuses, steps) == ([row["status"]], [row["steps"]]), (budget, row)
            if row["status"] == "complete":
                assert scores == [row["score"]], (budget, row)
            else:
                assert math.isnan(scores[0]), (budget, row)
        rows = [example["parameters"], [x * 1.1 for x in example["parameters"]],
                [9.0] * len(example["parameters"])]
        listed = prepared.scores([{"id": str(i), "parameters": r} for i, r in enumerate(rows)],
                                 **budget)["records"]
        scores, statuses, steps = prepared.score_matrix(rows, **budget)
        assert statuses == [r["status"] for r in listed] and steps == [r["steps"] for r in listed]
        assert all(s == r["score"] if r["status"] == "complete" else math.isnan(s)
                   for s, r in zip(scores, listed, strict=True))
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
try:
    prepared.score_matrix([example["parameters"], example["parameters"][:1]])
except ValueError:
    pass
else:
    raise AssertionError("score_matrix accepted ragged rows")
try:
    import numpy
except ImportError:
    numpy = None
if numpy is not None:
    array = numpy.array([example["parameters"]] * 3)
    assert prepared.score_matrix(array) == prepared.score_matrix(array.tolist())



def run_blocks(path, label, minimum):
    snippets = re.findall(r"^```python\n(.*?)^```", Path(path).read_text(), re.M | re.S)
    assert len(snippets) >= minimum, (path, len(snippets))
    namespace = {}
    for i, snippet in enumerate(snippets, 1):
        exec(compile(snippet, f"{label}-{i}", "exec"), namespace)
    return len(snippets)


snippets = run_blocks(tutorial, "tutorial", 4)
manual_blocks = run_blocks(manual, "manual", 15)
print(json.dumps({"status": "passed", "tutorial_python_blocks": snippets,
                  "manual_python_blocks": manual_blocks,
                  "lean_api_and_budget_parity": True, "matrix_api_parity": True,
                  "numpy_matrix_input": numpy is not None}))
