"""Run or resume the local campaign manifest against the installed module."""

import importlib
import json
from pathlib import Path
import sys

module, config, database = sys.argv[1:]
manifest = json.loads(Path(config).read_text())
campaign = importlib.import_module(f"{module}.campaign")
print(json.dumps(campaign.run(database, manifest["scenarios"], manifest["candidates"], workers=4), sort_keys=True))
