"""Replay CLI artifacts, schema validation and conditional scope display."""

import contextlib
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from xpsloss import DeltaPrimary, cli, loss_window_area, single_loss_spectrum
from xpsloss import xps_single_loss as api
from xpsloss.config import load_config


class SingleLossCLITests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.config = self.root/"input.json"
        self.spec = self.root/"spectrum.json"
        self.out = self.root/"out"
        self.config.write_text(json.dumps(dict(
            material=dict(name="synthetic", gap_eV=0, oscillators=[
                dict(energy_eV=20, strength_eV2=200, width_eV=8, dispersion=0)]),
            electron_energy_eV=1000, loss_grid_eV=dict(min=1, max=100, step=1),
            integration=dict(epsabs=1e-12, epsrel=1e-9, limit=200))), encoding="utf-8")
        self.primary = dict(kind="delta", energies_eV=[1000], intensities=[1])
        self.spec_input = dict(schema="xpsloss-single-loss-spectrum-v0", axis="kinetic_energy_eV",
                               primary=self.primary, output_edges_eV={"min": 960, "max": 1000, "step": 8})
        self.spec.write_text(json.dumps(self.spec_input), encoding="utf-8")

    def invoke(self, spectrum, *extra, window=("1", "50")):
        args = ["xps-single-loss" if spectrum else "xps-loss-window-area", str(self.config),
                "--loss-window-eV", *window, "--weight-length-A", "20", "--epsilon-momentum", "k_total",
                "--weight-length-source", "example", "--out", str(self.out)]
        if spectrum:
            args += ["--spectrum", str(self.spec)]
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = cli.main([*args, *extra])
        data = stderr.getvalue() if status else stdout.getvalue()
        self.assertTrue(data.isascii())
        return status, json.loads(data)

    def replay(self, spectrum):
        status, data = self.invoke(spectrum, "--kernel-limit", "201", window=("15", "25"))
        self.assertEqual(status, 0, data)
        run = json.loads((self.out/"run.json").read_text(encoding="utf-8"))
        self.assertEqual(run["status"], "complete")
        self.assertFalse(run["spectrum_definition"]["config_loss_grid_used"])
        self.assertEqual(run["effective_settings_sources"]["kernel_limit"], "cli")
        self.assertEqual(run["effective_settings_sources"]["kernel_epsrel"], "config")
        self.assertEqual(run["effective_settings_sources"]["epsrel"], "default")
        self.assertEqual(run["effective_settings"], run["diagnostics"]["effective_settings"])
        self.assertEqual(len(run["effective_settings"]), 8 if spectrum else 6)
        m = load_config(self.config).material
        definition = run["spectrum_definition"]
        kw = {k: definition[k] for k in ["loss_window_eV", "hole_charge", "epsilon_momentum", "theta_deg"]}
        kw.update(weight_length_A=definition["weight"]["length_A"], **run["effective_settings"])
        if spectrum:
            prim = run["primary"]
            result = single_loss_spectrum(m, 1000, DeltaPrimary(prim["energies_eV"], prim["intensities"]),
                                          run["output_edges_eV"], **kw)
            expected_edges = np.arange(960, 1001, 8, dtype="<f8")
            self.assertEqual(run["spectrum_input"]["output_edges_sha256_float64_le"],
                             hashlib.sha256(expected_edges.tobytes()).hexdigest())
            self.assertEqual(run["spectrum_input"]["sha256"], hashlib.sha256(self.spec.read_bytes()).hexdigest())
            with (self.out/"single_loss_spectrum.csv").open(encoding="ascii", newline="") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([r["region"] for r in rows], ["low", "0", "1", "2", "3", "4", "high"])
            self.assertEqual((rows[0]["edge_low_eV"], rows[-1]["edge_high_eV"]), ("", ""))
            for i, row in enumerate(rows):
                outside = result.outside.get(row["region"])
                for q in api.QUANTITIES:
                    for suffix in ("", "_error"):
                        k = q+suffix
                        expected = outside[k] if outside else getattr(result, k)[i-1]
                        self.assertEqual(float(row[k]).hex(), float(expected).hex())
                for sign in ("plus", "minus"):
                    for c in api.COMPONENTS:
                        for suffix in ("", "_error"):
                            k = f"cutoff_{sign}_{c}"+suffix
                            expected = (outside["cutoff_sensitivity"+suffix][sign][c] if outside else
                                        getattr(result, "cutoff_sensitivity_bins"+suffix)[sign][c][i-1])
                            self.assertEqual(float(row[k]).hex(), float(expected).hex())
        else:
            result = loss_window_area(m, 1000, **kw)
            self.assertIsNone(run["primary"])
            self.assertIsNone(run["spectrum_input"])
            self.assertIsNone(run["output_edges_eV"])
        with (self.out/"window_area.csv").open(encoding="ascii", newline="") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual([row["column"] for row in rows], list(api.COLUMNS))
        for row in rows:
            for k, v in list(row.items())[1:]:
                expected = getattr(result, k.removesuffix("_per_A"))[row["column"]]
                self.assertEqual(float(v).hex(), expected.hex())
        for k, v in run["results"].items():
            self.assertEqual(v, cli._single_loss_cli.json_value(getattr(result, k)))
        self.assertEqual(run["comparison"], run["diagnostics"]["comparison"])

    def test_area_manifest_and_csv_replay(self):
        self.replay(False)

    def test_spectrum_manifest_and_both_csvs_replay(self):
        self.replay(True)

    def test_comparison_requires_every_key(self):
        cond = {"electron_energy_eV": 1000, "output_edges_eV": [1, 2], "primary": self.primary}
        record = dict(entry="single_loss_spectrum", conditions=cond, key_status={})
        table = dict(version="test", registration_sha256="", reference_sha256=None, records=[record])
        with patch.object(api, "_load_comparison_table", return_value=table):
            def compare(conditions):
                return api._comparison(conditions, conditions["output_edges_eV"], "fixed_poles_enumerated")
            keys = compare(cond)["key_status"]
            record["key_status"] = {g: dict.fromkeys(keys[g], dict(status="pass", reason=None)) for g in keys}
            self.assertEqual(compare(cond)["spectrum"], "matched")
            del record["key_status"]["spectrum"]["outside.low.single_loss_total"]
            c = compare(cond)
            self.assertEqual((c["areas"], c["spectrum"]), ("matched", "not_matched"))
            self.assertEqual(c["key_status"]["spectrum"]["outside.low.single_loss_total"]["status"], "no_reference")
            self.assertEqual(compare({**cond, "output_edges_eV": [1, 3]})["spectrum"], "kernel_conditions_checked")
            self.assertEqual(compare({**cond, "electron_energy_eV": 900})["areas"], "not_matched")

    def test_complete_table_resource_format(self):
        self.assertEqual(set(api._load_comparison_table()),
                         {"version", "registration_sha256", "reference_sha256", "records"})

    def test_failure_json_and_absent_output(self):
        rc, err = self.invoke(False, "--max-cells", "10")
        self.assertEqual((rc, err["error_type"]), (2, "ValueError"))
        self.assertNotIn("reason", err)
        rc, err = self.invoke(True, window=("1", "nan"))
        self.assertEqual((rc, err["reason"], err["loss_window_eV"]), (2, "invalid_window", [1., "nan"]))
        for invalid in [{**self.spec_input, "axis": "binding_energy_eV"},
                        {**self.spec_input, "primary": {**self.primary, "edges_eV": [1, 2]}},
                        {**self.spec_input, "output_edges_eV": {"min": 0, "max": 1, "step": .3}}]:
            self.spec.write_text(json.dumps(invalid), encoding="utf-8")
            rc, err = self.invoke(True)
            self.assertEqual((rc, err["reason"]), (2, "invalid_schema"))
        self.assertFalse(self.out.exists())

