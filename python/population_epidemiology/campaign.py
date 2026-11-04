"""Local, resumable candidate/scenario evaluation with atomic result commits."""

import hashlib
import json
import sqlite3
from pathlib import Path

from . import Prepared, _encode
from . import _native


def digest(value):
    return hashlib.sha256(_encode(value).encode()).hexdigest()


def run(database, scenarios, candidates, *, workers=1, chunk_size=16, **budget):
    """Evaluate a fixed scenario/candidate set and return rows in manifest order.

    Failed and cancelled attempts remain in the journal. Fixed-work complete,
    invalid and infeasible results may be replayed. Work settings are part of
    the identity. A later attempt can retry cancelled or time-limited records.
    """
    if not 1 <= chunk_size <= 256:
        raise ValueError("invalid chunk size")
    for records in (scenarios, candidates):
        ids = [r["id"] for r in records]
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate scenario or candidate identifier")
    build = hashlib.sha256(Path(_native.__file__).read_bytes()).hexdigest()
    # Runtime interrupt controls do not change a completed numerical result.
    settings = {k: v for k, v in budget.items() if k != "cancel_file"}
    result = []
    with sqlite3.connect(database, timeout=30) as conn:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("CREATE TABLE IF NOT EXISTS results (key TEXT PRIMARY KEY, payload TEXT NOT NULL)")
        conn.execute("CREATE TABLE IF NOT EXISTS attempts (seq INTEGER PRIMARY KEY, key TEXT, payload TEXT)")
        for scenario in scenarios:
            model = scenario["model"]
            prepared = Prepared(model)
            keys = [
                digest([build, scenario["id"], model, c, settings]) for c in candidates
            ]
            rows = [None] * len(keys)
            for start in range(0, len(keys), chunk_size):
                missing = []
                for i in range(start, min(start + chunk_size, len(keys))):
                    hit = conn.execute("SELECT payload FROM results WHERE key=?", (keys[i],)).fetchone()
                    if hit is None:
                        missing.append(i)
                    else:
                        rows[i] = json.loads(hit[0])
                if not missing:
                    continue
                evaluated = prepared.batch([candidates[i] for i in missing], workers=workers, **budget)
                with conn:
                    for i, record in zip(missing, evaluated["records"], strict=True):
                        row = {"scenario": scenario["id"], "key": keys[i], "build_sha256": build,
                               "model_sha256": evaluated["model_sha256"], "record": record}
                        payload = _encode(row)
                        conn.execute("INSERT INTO attempts (key, payload) VALUES (?,?)", (keys[i], payload))
                        if record["status"] in {"complete", "invalid", "infeasible"}:
                            conn.execute("INSERT OR IGNORE INTO results VALUES (?,?)", (keys[i], payload))
                        rows[i] = row
            result.extend(rows)
    return result
