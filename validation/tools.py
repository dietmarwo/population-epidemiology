"""Shared test helpers. Numerical reference logic lives in each domain."""

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAME = next(line.split('"')[1] for line in (ROOT / "Cargo.toml").read_text().splitlines() if line.startswith("name = "))


def run(binary, model, params, *, trace=False, budget=None):
    candidates = [{"id": str(i), "parameters": x, "trace": trace} for i, x in enumerate(params)]
    req = {"schema": f"{NAME}.request/0.1", "model": model, "candidates": candidates, "budget": budget or {}}
    result = subprocess.run([binary, "run", "-"], input=json.dumps(req), capture_output=True, text=True, timeout=90, check=False)
    if result.returncode not in (0, 2):
        raise RuntimeError(result.stderr)
    rows = json.loads(result.stdout)["records"]
    assert [r["id"] for r in rows] == [c["id"] for c in candidates]
    return rows


def good(row):
    assert row["status"] == "complete", row
    return row["metrics"]


def save(example, report):
    (ROOT / "models" / "example.json").write_text(json.dumps(example, indent=2) + "\n")
    (ROOT / "evidence" / "validation.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
