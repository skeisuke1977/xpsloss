"""Run the same unittest discovery in the source tree and unpacked sdist."""
import os
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def use_source_tree() -> Path:
    """Select this tree for discovery and for children that inherit the environment."""
    init = SRC / "xpsloss" / "__init__.py"
    if not init.is_file():
        print(f"Missing source package: {ascii(str(init))}", file=sys.stderr)
        raise SystemExit(1)

    # Benchmark modules are imported from ROOT, and xpsloss from SRC.
    paths = [str(SRC), str(ROOT)]
    if sys.path[:2] != paths:
        sys.path[:0] = paths
    pythonpath = os.environ.get("PYTHONPATH", "")
    if pythonpath.split(os.pathsep, 1)[0] != str(SRC):
        os.environ["PYTHONPATH"] = str(SRC) + os.pathsep + pythonpath if pythonpath else str(SRC)

    import xpsloss

    loaded = Path(xpsloss.__file__).resolve()
    if not loaded.is_relative_to(SRC.resolve()):
        print(f"Unexpected xpsloss at {ascii(str(loaded))}: it was imported before "
              f"the runner selected {ascii(str(SRC))}. Run tests/run_suite.py in a new "
              "process that does not import xpsloss first (check .pth and sitecustomize files).",
              file=sys.stderr)
        raise SystemExit(1)
    print(f"xpsloss under test: {ascii(str(loaded))}", file=sys.stderr)
    return loaded


if __name__ == "__main__":
    use_source_tree()
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parent))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.testsRun > 0 and result.wasSuccessful() else 1)
