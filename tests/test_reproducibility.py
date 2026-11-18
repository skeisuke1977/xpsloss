"""A saved configuration must be directly reusable with strict input validation."""

import json
from pathlib import Path
import unittest

from xpsloss.config import load_config, parse_config


class ReproducibilityTests(unittest.TestCase):
    def test_saved_configuration_can_be_reloaded(self):
        example = Path(__file__).resolve().parents[1] / "examples" / "synthetic_bulk.json"
        original = load_config(example)
        serialized = json.dumps(original.public_dict(), allow_nan=False)
        reloaded = parse_config(json.loads(serialized))
        self.assertEqual(original.public_dict(), reloaded.public_dict())
        self.assertEqual(original.grid().tolist(), reloaded.grid().tolist())


if __name__ == "__main__":
    unittest.main()
