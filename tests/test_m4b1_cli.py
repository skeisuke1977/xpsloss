"""CLI round trip for integrated detector bins, order areas and error records."""

import contextlib
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from xpsloss import cli, DeltaPrimary, MULTIPLE_LOSS_MODEL_ID, MultipleLossConvergenceError, multiple_loss_spectrum
from xpsloss.config import load_config
from xpsloss._single_loss_types import COLUMNS, COMPONENTS


class MultipleCLI(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.out = self.root / "run"
        self.cfg = self.root / "material.json"
        self.spec = self.root / "spectrum.json"
        self.cfg.write_text(
            json.dumps(
                dict(
                    material=dict(
                        name="synthetic",
                        gap_eV=0,
                        oscillators=[dict(energy_eV=20, strength_eV2=200, width_eV=8, dispersion=0)],
                    ),
                    electron_energy_eV=1000,
                    loss_grid_eV=dict(min=1, max=100, step=1),
                    integration=dict(epsrel=1e-9, epsabs=1e-12, limit=200),
                )
            )
        )
        self.input = dict(
            schema="xpsloss-single-loss-spectrum-v0",
            axis="kinetic_energy_eV",
            primary=dict(kind="delta", energies_eV=[1000], intensities=[1]),
            output_edges_eV=[997.5, 998.5, 999.5, 1000.5],
        )
        self.spec.write_text(json.dumps(self.input))

    def invoke(self, *extra):
        out, err = io.StringIO(), io.StringIO()
        args = [
            "xps-multiple-loss",
            str(self.cfg),
            "--spectrum",
            str(self.spec),
            "--loss-window-eV",
            "1",
            "50",
            "--weight-length-A",
            "20",
            "--epsilon-momentum",
            "k_total",
            "--model",
            MULTIPLE_LOSS_MODEL_ID,
            "--out",
            str(self.out),
            *extra,
        ]
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = cli.main(args)
        return status, json.loads((err if status else out).getvalue())

    def test_csv_and_manifest(self):
        status, info = self.invoke("--max-order", "display_support")
        self.assertEqual(status, 0, info)
        manifest = json.loads((self.out / "run.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["truncation"]["computed"], 2)
        self.assertEqual(manifest["truncation"]["state"], "complete_on_display")
        self.assertEqual(len(manifest["effective_settings"]), 13)
        self.assertEqual(manifest["spectrum_input"]["sha256"], hashlib.sha256(self.spec.read_bytes()).hexdigest())
        self.assertEqual(manifest["effective_settings_sources"]["kernel_epsrel"], "config")
        cfg = load_config(self.cfg)
        definition = manifest["spectrum_definition"]
        replay = multiple_loss_spectrum(
            cfg.material,
            cfg.electron_energy_eV,
            DeltaPrimary(manifest["primary"]["energies_eV"], manifest["primary"]["intensities"]),
            manifest["output_edges_eV"],
            loss_window_eV=definition["loss_window_eV"],
            weight_length_A=definition["weight"]["length_A"],
            epsilon_momentum=definition["epsilon_momentum"],
            theta_deg=definition["theta_deg"],
            hole_charge=definition["hole_charge"],
            multiplier=definition["multiplier"]["convention"],
            model=manifest["model_id"],
            max_order=definition["max_order"],
            **manifest["effective_settings"],
        )
        self.assertEqual(list(replay.warnings), manifest["warnings"])
        self.assertEqual(replay.spectrum_error.tolist(), manifest["results"]["spectrum_error"])
        with (self.out / "multiple_loss_spectrum.csv").open(encoding="ascii", newline="") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual([r["region"] for r in rows], ["low", "0", "1", "2", "high"])
        for o in manifest["results"]["orders"]:
            for c in COMPONENTS:
                for suffix in ("", "_error"):
                    key = c + suffix
                    expected = [o["outside"]["low"][key], *o[key], o["outside"]["high"][key]]
                    self.assertEqual(
                        [float(r[f"order{o['order']}_{key}"]).hex() for r in rows], [v.hex() for v in expected]
                    )
        with (self.out / "order_areas.csv").open(encoding="ascii", newline="") as f:
            areas = list(csv.DictReader(f))
        self.assertEqual([row["column"] for row in areas], list(COLUMNS) * 2)
        for row in areas:
            order = replay.orders[int(row["order"]) - 1]
            for key in ("area", "area_error", "cutoff_sensitivity", "cutoff_sensitivity_error"):
                self.assertEqual(float(row[key]).hex(), getattr(order, key)[row["column"]].hex())
        self.assertEqual(self.invoke("--max-order", "1")[1]["error_type"], "FileExistsError")

    def test_rejection_records(self):
        status, err = self.invoke("--max-order", "65")
        self.assertEqual((status, err["reason"], err["requested_order"]), (2, "order_limit_exceeded", 65))
        self.input["axis"] = "binding_energy_eV"
        self.spec.write_text(json.dumps(self.input))
        status, err = self.invoke("--max-order", "1")
        self.assertEqual(
            (status, err["error_type"], err["reason"]), (2, "MultipleLossInputError", "invalid_schema")
        )
        self.assertFalse(self.out.exists())

    def test_capabilities(self):
        caps = cli.capabilities()
        self.assertEqual(caps["supported"][-1], "xps_multiple_loss")
        self.assertEqual(caps["multiple_scattering"]["model_id"], MULTIPLE_LOSS_MODEL_ID)
        self.assertIn("binding_energy_axis", caps["multiple_scattering"]["unsupported_physics"])

    def test_inherited_budget_message_keeps_machine_readable_suggestion(self):
        why = MultipleLossConvergenceError("inherited_error_exceeds_order_budget", order=2,
            stage="order_fold", key="bin[2].hole", side="absolute", suggestion="spectrum_epsabs_fraction")
        with patch("xpsloss.xps_multiple_loss.multiple_loss_spectrum", side_effect=why):
            status, err = self.invoke("--max-order", "2")
        self.assertEqual((status, err["status"], err["reason"]), (2, "not_converged", why.reason))
        self.assertEqual(err["message"], str(why))
        self.assertIn("lowering epsrel", err["message"])
        self.assertIn("1e-7", err["message"])
        self.assertEqual((err["side"], err["suggestion"]), ("absolute", "spectrum_epsabs_fraction"))
        self.assertFalse(self.out.exists())
