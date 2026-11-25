"""V15: actual GC, effective settings, provenance, CSV and byte encoding."""

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from xpsloss import cli
from xpsloss.config import load_config
from xpsloss import xps_depth_averaged as avg
from test_m3a_physics import GRIDS, ROOT, average, evaluate, losses, material, table


def invocation(row, folder):
    cfg = folder/"input.json"
    out = folder/(table()["cli_strings"]["ja_output_dir"] if row.get("output_dir_name") == "ja" else "output")
    obj = dict(material=asdict(material(row["material"])), electron_energy_eV=row["electron_energy_eV"],
               loss_grid_eV=row["cli_config"]["loss_grid_eV"])
    prec = row["cli_config"]["integration"]
    if prec is not None:
        obj["integration"] = table()["settings"][prec]
    cfg.write_bytes(json.dumps(obj).encode("ascii"))
    parsed = load_config(cfg)
    np.testing.assert_array_equal(parsed.grid(), losses(row))
    for k, v in table()["settings"][row["settings"]].items():
        if getattr(parsed, k) != v:
            raise AssertionError(k)
    args = ["xps-depth-averaged", str(cfg), "--weight-length-A", str(row["weight_length_A"]["value"]),
            "--epsilon-momentum", row["epsilon_momentum"], "--theta-deg", str(row["theta_deg"]),
            "--hole-charge", str(row["hole_charge"]["value"]), "--out", str(out)]
    if row["weight_length_source"] is not None:
        args.extend(["--weight-length-source", row["weight_length_source"]])
    return args, cfg, out


def run_cli_capture(row):
    with tempfile.TemporaryDirectory() as tmp:
        args, cfg, out = invocation(row, Path(tmp))
        captured = []
        def calculate(*args, **kwargs):
            result = average(*args, **kwargs)
            captured.append((result, kwargs))
            return result
        with patch.object(cli, "effective_cross_section_depth_averaged", side_effect=calculate), patch.object(cli, "emit"):
            if cli.main(args):
                raise AssertionError(row["id"])
        (result, kwargs), = captured
        manifest = json.loads((out/"run.json").read_bytes())
        for k, value in table()["settings"][row["settings"]].items():
            if kwargs[k] != value or manifest["diagnostics"]["effective_settings"][k] != value:
                raise AssertionError((k, kwargs, manifest["diagnostics"]["effective_settings"]))
        if manifest["input_sha256"] != hashlib.sha256(cfg.read_bytes()).hexdigest():
            raise AssertionError("input bytes differ")
        csv = np.genfromtxt(out/"xps_depth_averaged.csv", delimiter=",", names=True)
        names = ("loss_eV", "total", "no_hole", "hole", *avg.fixed.TERM_NAMES, "quadrature_error")
        if csv.dtype.names != names:
            raise AssertionError(csv.dtype.names)
        for name in names:
            if csv[name].tobytes() != getattr(result, name).tobytes():
                raise AssertionError(name)
        if row["weight_length_source"] and row["weight_length_source"].encode("utf8") not in (out/"run.json").read_bytes():
            raise AssertionError("missing UTF-8 source")
        return result, manifest


class AverageCLITests(unittest.TestCase):
    def test_v15_ms4_cli_api_and_manifest(self):
        scope = cli.capabilities()["xps_depth_averaged"]
        self.assertEqual(scope["model_id"], "simonsen1997_depth_average_external_v0")
        # Review finding 5: the fixed-depth limits on sub-gap delta terms and broadening carry over.
        for name in ("sub_gap_surface_mode_delta_terms", "gap_broadening", "oblique_emission"):
            self.assertIn(name, scope["unsupported_physics"])
        for row in table()["runs"]:
            if row["family"] != "MS4" or row["entry"] != "CLI":
                continue
            with self.subTest(run=row["id"]):
                r, m = run_cli_capture(row)
                independent_call = evaluate(row)
                settings = table()["settings"][row["settings"]]
                self.assertEqual(m["diagnostics"]["effective_settings"], settings)
                self.assertEqual(m["model_id"], scope["model_id"])
                self.assertEqual(m["fixed_depth_model_id"], cli.FIXED_DEPTH_MODEL_ID)
                self.assertEqual(m["dielectric_model_id"], cli.DIELECTRIC_MODEL_ID)
                self.assertEqual(m["weight"], dict(form="x*exp(-x/L)/L^2", normalization=1, length_A=20.,
                    length_kind="external_constant", length_source=row["weight_length_source"],
                    depth_integration="closed_form_laplace", elastic_scattering=False,
                    depth_distribution="uniform_semi_infinite", loss_window_for_length=None))
                self.assertFalse(m["applicability"]["assumptions"]["elastic_scattering"])
                self.assertEqual(m["warnings"], list(r.warnings))
                self.assertEqual(any("low-loss" in w for w in r.warnings), row["material"] != "Z")
                diag = m["diagnostics"]
                errors, ratios = [], []
                for c in avg.COLUMNS:
                    self.assertEqual(getattr(r, c).tobytes(), getattr(independent_call, c).tobytes())
                    self.assertEqual(m["units"][c], "angstrom^-1 eV^-1")
                    errors.append(r.component_errors[c])
                    ratios.append(errors[-1]/np.maximum(settings["epsabs"], settings["epsrel"]*abs(getattr(r, c))))
                    np.testing.assert_array_equal(diag["component_errors"][c], errors[-1])
                    np.testing.assert_array_equal(diag["error_budget_ratios"][c], ratios[-1])
                self.assertEqual(m["units"]["loss_eV"], "eV")
                self.assertEqual(m["units"]["quadrature_error"], "angstrom^-1 eV^-1")
                np.testing.assert_array_equal(r.quadrature_error, np.max(errors, axis=0))
                c, i = np.unravel_index(np.argmax(ratios), np.shape(ratios))
                self.assertEqual(diag["max_budget_ratio"], dict(value=ratios[c][i], column=avg.COLUMNS[c], index=i))
                self.assertEqual(diag["direct_sum"], r.diagnostics["direct_sum"])
                self.assertEqual(diag["direct_sum_max"], r.diagnostics["direct_sum_max"])
                self.assertEqual(any(d is not None for d in diag["direct_sum"]), row["material"] == "Si")

    def test_v15_api_default(self):
        kwargs = dict(weight_length_A=20., epsilon_momentum="k_total")
        args = (material("Al"), 1000., GRIDS["GC"])
        default = average(*args, **kwargs)
        explicit = average(*args, **kwargs, **table()["settings"]["DEF"])
        self.assertEqual(default.diagnostics["effective_settings"], table()["settings"]["DEF"])
        for c in (*avg.COLUMNS, "quadrature_error"):
            self.assertEqual(getattr(default, c).tobytes(), getattr(explicit, c).tobytes())
        for c in avg.COLUMNS:
            self.assertEqual(default.component_errors[c].tobytes(), explicit.component_errors[c].tobytes())

    def test_cli_bytes_failures_and_existing_output(self):
        for row in table()["runs"]:
            if row["entry"] != "CLI" or not (row.get("output_dir_name") == "ja" or row["family"] in ("NC2", "RJ3")):
                continue
            with self.subTest(run=row["id"]), tempfile.TemporaryDirectory() as tmp:
                args, _, out = invocation(row, Path(tmp))
                if row["family"] == "RJ5":
                    out.mkdir()
                    (out/"marker").write_bytes(bytes(range(256)))
                    before = hashlib.sha256((out/"marker").read_bytes()).digest()
                proc = subprocess.run([sys.executable, "-m", "xpsloss", *args], cwd=ROOT, capture_output=True)
                self.assertEqual(proc.returncode, row["expect"]["exit_code"])
                self.assertTrue(proc.stdout.isascii() and proc.stderr.isascii())
                record = json.loads(proc.stderr if proc.returncode else proc.stdout)
                self.assertEqual(record["status"], row["expect"]["status"])
                if proc.returncode == 0:
                    self.assertEqual(record["output_dir"], str(out.resolve()))
                elif row["family"] == "RJ5":
                    self.assertEqual(record["error_type"], "FileExistsError")
                    self.assertIn(str(out), record["message"])
                    self.assertEqual(list(out.iterdir()), [out/"marker"])
                    self.assertEqual(hashlib.sha256((out/"marker").read_bytes()).digest(), before)
                else:
                    self.assertFalse(out.exists())
                    self.assertEqual(record["reason"], row["expect"]["reason"])
                    if row["family"] == "NC2":
                        self.assertEqual(record["error_type"], "DepthAverageConvergenceError")
                        self.assertEqual((record["loss_eV"], record["index"], record["subdivisions"]), (5., 0, 0))
                    else:
                        self.assertEqual(record["error_type"], "FixedDepthGapDomainError")
                        self.assertEqual(record["details"]["points"][0]["index"], 0)
