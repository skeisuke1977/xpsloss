"""Round-trip user text through the CLI on redirected cp932 streams."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
NAME = "synthetic_\u00e9\u2603\U0001f600\u96ea"


class CLIEncodingTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.directory = Path(tmp.name)
        example = ROOT / "examples" / "synthetic_bulk.json"
        config = json.loads(example.read_text(encoding="utf-8"))
        config["material"]["name"] = NAME
        (self.directory / "input.json").write_text(json.dumps(config, ensure_ascii=False), encoding="utf-8")

    def call(self, *args):
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONUTF8="0", PYTHONIOENCODING="cp932")
        return subprocess.run([sys.executable, "-m", "xpsloss", *args], cwd=self.directory,
                              env=env, capture_output=True)

    def test_validate_returns_original_material_name(self):
        proc = self.call("validate", "input.json")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stderr, b"")
        result = json.loads(proc.stdout.decode("ascii"))
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["config"]["material"]["name"], NAME)

    def test_calculation_returns_original_path_and_writes_utf8_manifest(self):
        proc = self.call("bulk", "input.json", "--out", NAME)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stderr, b"")
        result = json.loads(proc.stdout.decode("ascii"))
        output = (self.directory / NAME).resolve()
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["output_dir"], str(output))
        self.assertTrue((output / "bulk.csv").is_file())
        manifest = (output / "run.json").read_text(encoding="utf-8")
        self.assertIn(NAME, manifest)
        self.assertEqual(json.loads(manifest)["config"]["material"]["name"], NAME)

    def test_error_contains_original_missing_filename(self):
        missing = NAME + ".json"
        proc = self.call("validate", missing)
        self.assertEqual(proc.returncode, 2)
        self.assertEqual(proc.stdout, b"")
        error = json.loads(proc.stderr.decode("ascii"))
        self.assertEqual(error["status"], "error")
        self.assertEqual(error["error_type"], "FileNotFoundError")
        self.assertIn(missing, error["message"])
