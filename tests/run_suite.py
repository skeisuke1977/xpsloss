"""Run the same unittest discovery in the source tree and unpacked sdist."""
from pathlib import Path
import sys
import unittest


if __name__ == "__main__":
    # Discovery by script path must also expose the source-tree benchmark modules.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parent))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.testsRun > 0 and result.wasSuccessful() else 1)
