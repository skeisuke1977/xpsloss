"""V14(b): full fixed-depth API and CLI comparison with main fd31dcd."""

from dataclasses import asdict, fields
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile

import numpy as np

import xpsloss
from xpsloss import Material, Oscillator, effective_cross_section_fixed_depth
from xpsloss.config import parse_config
from test_m22_migration import cases


ROOT = Path(__file__).resolve().parents[1]
BASE = "fd31dcd"
GA = 0.2*np.arange(1, 500)
DEF = dict(epsabs=1e-12, epsrel=1e-9, limit=200)


def canonical(x):
    if isinstance(x, np.ndarray):
        return [x.dtype.str, list(x.shape), hashlib.sha256(x.tobytes()).hexdigest()]
    if isinstance(x, (float, np.floating)):
        return float(x).hex()
    if isinstance(x, dict):
        return {k: canonical(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [canonical(v) for v in x]
    return x


def capture(material, energy, grid, mode, charge, settings):
    try:
        r = effective_cross_section_fixed_depth(material, energy, grid, 10., epsilon_momentum=mode,
                                                theta_deg=0., hole_charge=charge, **settings)
    except (ValueError, xpsloss.NumericalConvergenceError) as exc:
        return {"type": type(exc).__name__, "message": str(exc), "attrs": canonical(vars(exc))}
    return {f.name: canonical(getattr(r, f.name)) for f in fields(r)}


def snapshot(folder):
    records = {}
    for name, obj in cases():
        cfg = parse_config(obj)
        for mode in ("k_parallel", "k_total"):
            settings = dict(epsabs=cfg.epsabs, epsrel=cfg.epsrel, limit=cfg.limit)
            records[name + "/" + mode] = capture(cfg.material, cfg.electron_energy_eV, cfg.grid(), mode, 1., settings)
    assert len(records) == 40
    table = json.loads((ROOT/"benchmarks/numerics/m3a/expected_states_m3a.json").read_bytes())
    materials = {k: Material(k, v["gap_eV"], tuple(Oscillator(*o) for o in v["oscillators"]))
                 for k, v in table["materials"].items()}
    both = ("k_parallel", "k_total")
    edge_grid = np.r_[GA[:10], 1.12]
    for tag, mat, grid, modes, charges in (
        ("si", "Si", GA, both, (0, 1)), ("s0", "S0", GA, ("k_total",), (1,)),
        ("z0", "Z0", edge_grid, both, (1,)), ("z1", "Z1", GA, ("k_total",), (1,)),
        ("n2", "N2", 4., both, (0, 1)), ("n1", "N1", 1., ("k_total",), (1,)),
        ("p1", "P1", 4., ("k_total",), (1,)), ("p0", "P0", 4., ("k_total",), (1,)),
        ("ue", "UE", 2.5, both, (1,)), ("fail_n1", "N1", 4., both, (1,)),
        ("fail_p1", "P1", 1.999, both, (1,)), ("fail_m1", "M1", 2.5, both, (1,)),
        ("fail_ue", "UE", 1.5, ("k_parallel",), (1,)), ("fail_edge", "SP", 1., both, (1,)),
        ("fail_sign", "SP", 2.5, both, (1,)), ("fail_si", "Si", edge_grid, both, (1,)),
    ):
        for mode in modes:
            for q in charges:
                record = capture(materials[mat], 3. if mat == "UE" else 1000., grid, mode, q, DEF)
                if tag.startswith("fail") != ("type" in record):
                    raise AssertionError((tag, record.keys()))
                records[f"API/{tag}/{mode}/{q}"] = record

    for mat, lo, hi, modes in (("Si", .2, 99.8, ("k_total",)), ("Z0", 1.12, 1.12, both),
                                ("P1", 1.999, 1.999, ("k_total",))):
        for mode in modes:
            obj = dict(material=asdict(materials[mat]), electron_energy_eV=1000.,
                       loss_grid_eV=dict(min=lo, max=hi, step=.2), integration=DEF)
            cfg = folder/"input.json"
            cfg.write_bytes(json.dumps(obj).encode("ascii"))
            with tempfile.TemporaryDirectory(dir=folder) as tmp:
                out = Path(tmp)/"output"
                proc = subprocess.run([sys.executable, "-m", "xpsloss", "xps-fixed-depth", str(cfg),
                                       "--depth-A", "10", "--epsilon-momentum", mode, "--out", str(out)],
                                      capture_output=True, check=False)
                stdout, stderr = (json.loads(s) if s else None for s in (proc.stdout, proc.stderr))
                record = dict(exit=proc.returncode, stdout=stdout, stderr=stderr)
                if proc.returncode != (2 if mat == "P1" else 0):
                    raise AssertionError(record)
                if stdout:
                    stdout.pop("output_dir")
                    manifest = json.loads((out/"run.json").read_bytes())
                    if manifest.pop("package_version") != xpsloss.__version__:
                        raise AssertionError("wrong package version")
                    for key in ("created_utc", "elapsed_seconds"):
                        manifest.pop(key)
                    record.update(manifest=canonical(manifest),
                                  csv=hashlib.sha256((out/"xps_fixed_depth.csv").read_bytes()).hexdigest())
                elif out.exists():
                    raise AssertionError("partial output")
                records[f"CLI/{mat}/{mode}"] = record
    return records


class FixedDepthMigrationTests(unittest.TestCase):
    def test_all_fixed_depth_fields_exceptions_and_cli(self):
        if not (ROOT/".git").exists():
            self.skipTest("V14(b) requires the Git checkout, unavailable in an exported sdist")
        if not shutil.which("git") or subprocess.run(
                ["git", "cat-file", "-e", BASE], cwd=ROOT, capture_output=True).returncode:
            if "CI" in os.environ:
                self.fail(f"V14(b) requires Git history with {BASE}")
            self.skipTest(f"Git history with {BASE} unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            raw = subprocess.check_output(["git", "archive", "--format=zip", BASE, "src/xpsloss"], cwd=ROOT)
            with zipfile.ZipFile(io.BytesIO(raw)) as archive:
                archive.extractall(folder)
            snapshots = []
            for source in (folder/"src", ROOT/"src"):
                dest = folder/"snapshot.json"
                subprocess.run([sys.executable, str(Path(__file__).resolve()), "--snapshot", str(dest)],
                               cwd=ROOT, env=os.environ | {"PYTHONPATH": str(source)}, check=True)
                snapshots.append(json.loads(dest.read_bytes()))
            self.assertEqual(snapshots[0].keys(), snapshots[1].keys())
            for key in snapshots[0]:
                with self.subTest(key=key):
                    self.assertEqual(snapshots[0][key], snapshots[1][key])


if __name__ == "__main__" and sys.argv[1:2] == ["--snapshot"]:
    dest = Path(sys.argv[2])
    dest.write_text(json.dumps(snapshot(dest.parent)), encoding="utf8")
