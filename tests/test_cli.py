import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from xpsloss.config import load_config, parse_config


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "synthetic_bulk.json"


class CLITests(unittest.TestCase):
    def call(self, *args):
        env = dict(os.environ)
        env["PYTHONPATH"] = str(ROOT / "src")
        return subprocess.run([sys.executable, "-m", "xpsloss", *map(str, args)],
                              env=env, capture_output=True, text=True, check=False)

    def test_capabilities_and_validation(self):
        r = self.call("capabilities", "--json")
        self.assertEqual(r.returncode, 0, r.stderr)
        cap = json.loads(r.stdout)
        self.assertFalse(cap["queels_xps_equivalence"])
        self.assertIn("core_hole", cap["unsupported_physics"])
        r = self.call("validate", EXAMPLE)
        self.assertEqual(r.returncode, 0, r.stderr)
        result = json.loads(r.stdout)
        self.assertEqual(result["input_sha256"], hashlib.sha256(EXAMPLE.read_bytes()).hexdigest())

    def test_bulk_end_to_end_and_overwrite_protection(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "bulk"
            r = self.call("bulk", EXAMPLE, "--out", output)
            self.assertEqual(r.returncode, 0, r.stderr)
            result = json.loads(r.stdout)
            self.assertEqual(result["point_count"], 160)
            manifest = json.loads((output / "run.json").read_text())
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["units"]["diimfp_per_A_eV"], "angstrom^-1 eV^-1")
            data = np.loadtxt(output / "bulk.csv", delimiter=",", skiprows=1)
            np.testing.assert_allclose(data[:, 0], np.arange(0.5, 80.01, 0.5), rtol=0, atol=1e-12)
            self.assertTrue(np.all(np.isfinite(data)))
            self.assertTrue(np.all(data[:, 1:] >= 0))
            before = (output / "bulk.csv").read_bytes()
            r = self.call("bulk", EXAMPLE, "--out", output)
            self.assertEqual(r.returncode, 2)
            self.assertEqual(json.loads(r.stderr)["error_type"], "FileExistsError")
            self.assertEqual(before, (output / "bulk.csv").read_bytes())

    def test_elf_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "elf"
            r = self.call("elf", EXAMPLE, "--out", output)
            self.assertEqual(r.returncode, 0, r.stderr)
            manifest = json.loads((output / "run.json").read_text())
            self.assertTrue(manifest["diagnostics"]["optical_limit"])
            data = np.loadtxt(output / "elf.csv", delimiter=",", skiprows=1)
            self.assertEqual(data.shape, (160, 2))

    def test_unknown_and_invalid_config_rejected(self):
        base = json.loads(EXAMPLE.read_text())
        mutations = [
            lambda d: d.update(model="surface"),
            lambda d: d.update(surfaec=True),
            lambda d: d["material"].update(name=""),
            lambda d: d["material"].update(oscillators=[]),
            lambda d: d["material"]["oscillators"][0].update(strength_eV2=-1),
            lambda d: d["material"]["oscillators"][0].update(energy_eV=float("inf")),
            lambda d: d["material"]["oscillators"][0].update(width_eV=True),
            lambda d: d["loss_grid_eV"].update(min=0),
            lambda d: d["loss_grid_eV"].update(max=1000),
            lambda d: d["loss_grid_eV"].update(step=0.4),
            lambda d: d["integration"].update(epsrel=0),
            lambda d: d["integration"].update(limit=True),
        ]
        for i, mutate in enumerate(mutations):
            data = copy.deepcopy(base)
            mutate(data)
            with self.subTest(mutation=i):
                with self.assertRaises(ValueError):
                    parse_config(data)

    def test_nonfinite_and_duplicate_json_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            for text in ['{"model": "a", "model": "b"}', '{"x": NaN}', '{"x": Infinity}']:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    load_config(path)
            r = self.call("validate", path)
            self.assertEqual(r.returncode, 2)
            self.assertEqual(json.loads(r.stderr)["status"], "error")


if __name__ == "__main__":
    unittest.main()
