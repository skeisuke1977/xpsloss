"""Check which xpsloss the suite and its subprocesses actually import."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


RUNNER = Path(__file__).resolve().with_name("run_suite.py")
SRC = (RUNNER.parents[1] / "src").resolve()
REMEDY = "Run tests/run_suite.py, set PYTHONPATH=src, or use pip install -e ."
SHOW_PATH = "import json, xpsloss; print(json.dumps(xpsloss.__file__))"
LOAD_RUNNER = """
import importlib.util
import sys

spec = importlib.util.spec_from_file_location("run_suite", sys.argv[1])
run_suite = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run_suite)
"""


class SuiteSourceTests(unittest.TestCase):
    def test_same_process(self):
        try:
            import xpsloss
        except ModuleNotFoundError as exc:
            if exc.name != "xpsloss":
                raise
            self.fail(REMEDY)
        loaded = Path(xpsloss.__file__).resolve()
        self.assertTrue(loaded.is_relative_to(SRC), f"Loaded {ascii(str(loaded))}. {REMEDY}")

    def test_inheriting_child(self):
        child = subprocess.run([sys.executable, "-c", SHOW_PATH], capture_output=True)
        self.assertEqual(child.returncode, 0, f"{child.stderr!r}. {REMEDY}")
        loaded = Path(json.loads(child.stdout)).resolve()
        self.assertTrue(loaded.is_relative_to(SRC), f"Child loaded {ascii(str(loaded))}. {REMEDY}")


class RunnerWithDecoyTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        self.cwd = root / "empty"
        self.cwd.mkdir()
        decoy = root / "decoy-\u2603"
        self.init = decoy / "xpsloss" / "__init__.py"
        self.init.parent.mkdir(parents=True)
        self.init.write_text("", encoding="ascii")
        self.env = os.environ | {"PYTHONPATH": str(decoy)}

    def run_with_decoy(self, code, *args):
        return subprocess.run([sys.executable, "-c", code, *args], cwd=self.cwd,
                              env=self.env, capture_output=True)

    def test_runner_selects_source_for_child_and_grandchild(self):
        control = self.run_with_decoy(SHOW_PATH)
        self.assertEqual(control.returncode, 0, control.stderr)
        self.assertEqual(Path(json.loads(control.stdout)).resolve(), self.init.resolve())

        child = self.run_with_decoy(LOAD_RUNNER + """
import json
import subprocess

run_suite.use_source_tree()
import xpsloss
grandchild = subprocess.run([sys.executable, "-c", sys.argv[2]], check=True, stdout=subprocess.PIPE)
print(json.dumps([xpsloss.__file__, json.loads(grandchild.stdout)]))
""", str(RUNNER), SHOW_PATH)
        self.assertEqual(child.returncode, 0, child.stderr)
        child_path, grandchild_path = map(Path, json.loads(child.stdout))
        self.assertTrue(child_path.resolve().is_relative_to(SRC), f"Child loaded {ascii(str(child_path))}")
        self.assertTrue(grandchild_path.resolve().is_relative_to(SRC),
                        f"Grandchild loaded {ascii(str(grandchild_path))}")

    def test_runner_rejects_preloaded_decoy(self):
        child = self.run_with_decoy(LOAD_RUNNER + "\nimport xpsloss\nrun_suite.use_source_tree()", str(RUNNER))
        self.assertEqual(child.returncode, 1, child.stderr)
        self.assertIn(ascii(str(self.init.resolve())), child.stderr.decode("ascii"))
