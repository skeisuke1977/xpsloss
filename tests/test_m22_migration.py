"""V16: compare actual CLI bytes and binary64 API arrays in one environment."""

from dataclasses import asdict
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from xpsloss import cli, effective_cross_section_fixed_depth
from xpsloss.elfsource import parse_elf
from xpsloss.xps_fixed_depth import TERM_NAMES


ROOT = Path(__file__).resolve().parents[1]
RELEASE = "0167e380855f158358ecc2c2284e9cf6297f9914"
COLUMNS = (*TERM_NAMES, "no_hole", "hole", "total", "quadrature_error")


def cases():
    inputs = [ROOT / "examples/synthetic_bulk.json", *sorted((ROOT / "benchmarks/inputs").glob("*.bulk.json"))]
    for path in inputs:
        obj = json.loads(path.read_text(encoding="utf8"))
        if obj["material"]["gap_eV"]:
            obj["material"].update(name="synthetic gapless " + obj["material"]["name"], gap_eV=0.0)
            obj.pop("provenance", None)
        yield "shipped/" + obj["material"]["name"], obj
    raw = (ROOT / "benchmarks/sources/zenodo6024064_v1.0/ELF.txt").read_bytes()
    labels = re.findall(
        r"(?m)^([^@\[\]\r\n]+?)\s*\[[0-9, ]+\][ \t\r]*$", raw.decode("latin1").split("\nReferences", 1)[0]
    )
    for label in labels:
        name = label.strip()
        if name == "PPY - polypyrrole":  # Original print has no '@'; gapped.
            continue
        (mat,) = parse_elf(raw, material=name).materials
        if mat.gap_eV:
            continue
        lines = [{k: float(v) for k, v in asdict(o).items()} for o in mat.oscillators]
        yield (
            "ELF/" + name,
            {
                "material": {"name": name, "gap_eV": 0.0, "oscillators": lines},
                "electron_energy_eV": 1000.0,
                "loss_grid_eV": {"min": 1.0, "max": 100.0, "step": 1.0},
                "integration": {"epsabs": 1e-12, "epsrel": 1e-9, "limit": 200},
            },
        )
    yield (
        "synthetic unresolvable precision",
        {
            "material": {
                "name": "synthetic",
                "gap_eV": 0.0,
                "oscillators": [dict(energy_eV=15.0, strength_eV2=180.0, width_eV=2.0, dispersion=0.0)],
            },
            "electron_energy_eV": 1000.0,
            "loss_grid_eV": {"min": 15.0, "max": 15.0, "step": 1.0},
            "integration": {"epsabs": 1e-30, "epsrel": 1e-15, "limit": 30},
        },
    )


def snapshot():
    result = {}
    for name, obj in cases():
        for mode in ("k_parallel", "k_total"):
            with tempfile.TemporaryDirectory() as tmp:
                p, out = Path(tmp) / "input.json", Path(tmp) / "output"
                p.write_text(json.dumps(obj), encoding="utf8")
                record = {}

                def capture(*args, **kwargs):
                    r = effective_cross_section_fixed_depth(*args, **kwargs)
                    record["binary64_SHA256"] = {
                        k: hashlib.sha256(getattr(r, k).tobytes()).hexdigest() for k in COLUMNS
                    }
                    return r

                def status(message, *args):
                    record["status"] = message["status"]
                    if "error_type" in message:
                        record["error_type"] = message["error_type"]

                with (
                    patch.object(cli, "effective_cross_section_fixed_depth", side_effect=capture),
                    patch.object(cli, "emit", side_effect=status),
                ):
                    record["exit"] = cli.main(
                        [
                            "xps-fixed-depth",
                            str(p),
                            "--out",
                            str(out),
                            "--depth-A",
                            "10",
                            "--epsilon-momentum",
                            mode,
                        ]
                    )
                if record["exit"] == 0:
                    record["CSV_SHA256"] = hashlib.sha256(
                        (out / "xps_fixed_depth.csv").read_bytes()
                    ).hexdigest()
                result[name + " " + mode] = record
    return result


def compare_release(folder):
    raw = subprocess.check_output(["git", "archive", "--format=zip", RELEASE, "src/xpsloss"], cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        archive.extractall(folder / "release")
    snapshots = []
    for key, source in [("0167e38", folder / "release/src"), ("current", ROOT / "src")]:
        env = os.environ | {"PYTHONPATH": str(source)}
        dest = folder / (key + ".json")
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--snapshot", str(dest)],
            env=env,
            cwd=ROOT,
            check=True,
        )
        snapshots.append(json.loads(dest.read_text(encoding="utf8")))
    return snapshots


class MigrationTests(unittest.TestCase):
    def test_gapless_binary64_csv_and_failure_states(self):
        if not (ROOT / ".git").exists():
            self.skipTest("release-object comparison requires a Git checkout")
        if not shutil.which("git") or subprocess.run(
            ["git", "cat-file", "-e", RELEASE], cwd=ROOT, capture_output=True
        ).returncode:
            message = f"V16 requires release {RELEASE} in this checkout"
            if "CI" in os.environ:
                self.fail(message)
            self.skipTest(message)
        with tempfile.TemporaryDirectory() as folder:
            release, current = compare_release(Path(folder))
        self.assertEqual(len([k for k in release if k.startswith("ELF/")]), 30)
        self.assertEqual(len([k for k in release if k.startswith("shipped/")]), 8)
        self.assertEqual({r["status"] for r in release.values()}, {"complete", "error", "not_converged"})
        self.assertEqual(current, release)


if __name__ == "__main__" and sys.argv[1:2] == ["--snapshot"]:
    Path(sys.argv[2]).write_text(json.dumps(snapshot(), indent=2, sort_keys=True), encoding="utf8")
