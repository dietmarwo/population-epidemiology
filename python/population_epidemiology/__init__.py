"""In-process prepared evaluation with one native crossing per batch."""

import json
from importlib.resources import files

from ._native import Prepared as _Prepared


def _encode(value):
    return json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":"))


class Prepared:
    def __init__(self, model):
        self._native = _Prepared(_encode(model))

    def batch(self, candidates, **budget):
        return json.loads(self._native.batch(_encode(candidates), _encode(budget)))

    def scores(self, candidates, **budget):
        """Return ordered objectives and statuses for untraced candidates."""
        return json.loads(self._native.scores(_encode(candidates), _encode(budget)))

    def score_matrix(self, rows, **budget):
        """Score untraced parameter rows without JSON candidates or records.

        ``rows`` is a list of equal-length lists or a 2-D NumPy array, one row
        per candidate. Returns ``(scores, statuses, steps)`` as three lists in
        row order. A score is NaN unless its status is ``"complete"``. The
        budget keywords are the same as for ``scores``.
        """
        if hasattr(rows, "ndim"):
            if rows.ndim != 2:
                raise ValueError("score_matrix expects a 2-D array of parameter rows")
            rows = rows.tolist()
        scores, statuses, steps = self._native.score_matrix(rows, _encode(budget))
        return scores, statuses, steps

    def evaluate(self, parameters, *, trace=False, **budget):
        return self.batch([{"id": "one", "parameters": parameters, "trace": trace}], **budget)[
            "records"
        ][0]


def load_example():
    """Return the packaged model and a replayable parameter vector."""
    return json.loads(files(__package__).joinpath("example.json").read_text())
