"""Measured-spectrum JSON, reproducibility records and integrated-bin CSV."""

from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from xpsloss import cli, MULTIPLE_LOSS_MODEL_ID
from test_m4b2_axis import axis


class MeasuredSpectrumCLI(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.cfg, self.path = self.root/"config.json", self.root/"input.json"
        self.cfg.write_text(json.dumps(dict(material=dict(name="synthetic zero", gap_eV=0., oscillators=[
            dict(energy_eV=20., strength_eV2=0., width_eV=8., dispersion=0.)]), electron_energy_eV=32.,
            loss_grid_eV=dict(min=1., max=8., step=1.))), encoding="utf8")
        self.spec = dict(schema="xpsloss-measured-spectrum-v0", axis="binding_energy_eV",
                         measurement=asdict(axis()), primary=dict(kind="delta", binding_energies_eV=[-1., 0.],
                         intensities=[1., 3.]), output_binding_edges_eV=[-2., -1., 0., 2.])

    def invoke(self, spec, *, cmd="xps-multiple-loss", out="out", extra=()):
        raw = spec if isinstance(spec, str) else json.dumps(spec)
        self.path.write_text(raw, encoding="utf8")
        args = [cmd, str(self.cfg), "--spectrum", str(self.path), "--out", str(self.root/out),
                "--loss-window-eV", "1", "8", "--weight-length-A", "20", "--epsilon-momentum", "k_total", *extra]
        if cmd == "xps-multiple-loss":
            args += ["--model", MULTIPLE_LOSS_MODEL_ID, "--max-order", "1"]
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = cli.main(args)
        if code:
            self.assertFalse((self.root/out).exists())
            self.assertEqual(stdout.getvalue(), "")
            self.assertTrue(stderr.getvalue().isascii())
        return code, json.loads((stderr if code else stdout).getvalue())

    def measured(self):
        # In the file, even a null key belonging to the unselected kind is exclusive.
        self.spec["measurement"]["reference"] = dict(kind="photon_and_offset", kinetic_energy_offset_eV=0.)
        return self.spec

    def test_schema_before_meaning_before_arrays(self):
        spec = self.measured()
        spec["measurement"]["photon_energy_eV"] = "not_recorded"
        spec["primary"]["intensities"] = [-1.]
        for write in (lambda: dict(spec, observed_counts=[]),
                      lambda: dict(spec, output_binding_edges_eV=dict(min=0., max=1.)),
                      lambda: json.dumps(spec).replace('"axis":', '"axis":0,"axis":'),
                      lambda: json.dumps(spec).replace('"schema"', r'"\u0073chema"')[:-1],
                      lambda: json.dumps(spec)[:-1]):
            code, err = self.invoke(write())
            self.assertEqual((code, err["reason"], err["error_type"]), (2, "invalid_schema", "MultipleLossInputError"))
            self.assertTrue(err["detail"])
        # Only the schema identifier selects the measured reader: an unknown identifier (or none)
        # goes to the kinetic-energy reader and fails as in the M4b-1 baseline (no detail field).
        for other in (dict(spec, schema="unknown"), {k: v for k, v in spec.items() if k != "schema"}):
            code, err = self.invoke(other)
            self.assertEqual((code, err["reason"], err["error_type"]), (2, "invalid_schema", "MultipleLossInputError"))
            self.assertNotIn("detail", err)
        code, err = self.invoke(spec)
        self.assertEqual((code, err["reason"], err["field"], err["detail"]),
                         (2, "invalid_energy_axis", "photon_energy_eV", "not_recorded_not_allowed"))
        spec["measurement"]["photon_energy_eV"] = 32.
        self.assertEqual(self.invoke(spec)[1]["reason"], "invalid_primary")

    def test_reference_key_presence_and_command_scope(self):
        spec = self.measured()
        ref = spec["measurement"]["reference"]
        ref["binding_anchor_eV"] = None
        code, err = self.invoke(spec)
        self.assertEqual((code, err["field"], err["detail"]), (2, "reference", "exclusive"))
        del ref["binding_anchor_eV"]
        ref["kinetic_energy_offset_eV"] = None
        code, err = self.invoke(spec)
        self.assertEqual((code, err["field"], err["detail"]), (2, "reference.kinetic_energy_offset_eV", "missing"))
        code, err = self.invoke(spec, cmd="xps-single-loss")
        self.assertEqual((code, err["error_type"], err["reason"]), (2, "SingleLossInputError", "invalid_schema"))

    def test_edges_and_collision_error(self):
        spec = self.measured()
        spec["output_binding_edges_eV"] = dict(min=-2., max=2., step=.1)
        self.assertEqual(self.invoke(spec, extra=["--max-cells", "40"])[1]["reason"], "cell_limit_exceeded")
        spec["output_binding_edges_eV"] = [0., 2., 1.]
        self.assertEqual(self.invoke(spec)[1]["reason"], "invalid_output_edges")
        spec["output_binding_edges_eV"] = [0., 1., 2.]
        spec["measurement"]["photon_energy_eV"] = 8192.
        spec["primary"] = dict(kind="delta", binding_energies_eV=[1.+2.**-45], intensities=[0.])
        code, err = self.invoke(spec)
        self.assertEqual((code, err["reason"], err["roles"]), (2, "axis_conversion_not_injective",
                         ["output_edges_eV[1]", "primary.binding_energies_eV[0]"]))
        self.assertEqual(err["binding_eV"], [1., 1.+2.**-45])

    def test_two_directions_and_edge_digests(self):
        spec = self.measured()
        for name, edges in (("up", [-2., -1., 0., 1., 2.]), ("down", [2., 1., 0., -1., -2.])):
            spec["output_binding_edges_eV"] = dict(min=-2., max=2., step=1.) if name=="up" else edges
            status, reply = self.invoke(spec, out=name)
            self.assertEqual(status, 0, reply)
            folder = self.root/name
            m = json.loads((folder/"run.json").read_text(encoding="utf8"))
            a = m["energy_axis"]
            self.assertEqual(a, m["results"]["energy_axis"])
            self.assertEqual(m["measurement"], spec["measurement"])
            self.assertEqual(a["input_binding_edges_eV"], edges)
            self.assertEqual(m["spectrum_definition"]["axis"], "binding_energy_eV")
            self.assertEqual(set(a["input_primary"]), {"kind", "binding_energies_eV", "intensities", "intensity_unit",
                             "description", "included_broadening", "truncation", "intrinsic_satellites"})
            self.assertEqual(len(m["effective_settings"]), 13)
            digest = lambda arr: hashlib.sha256(np.array(arr, dtype="<f8").tobytes()).hexdigest()
            ke = sorted(32.-b for b in edges)
            self.assertEqual(m["spectrum_input"], dict(file=str(self.path), sha256=hashlib.sha256(self.path.read_bytes()).hexdigest(),
                             schema=spec["schema"], binding_output_edges_sha256_float64_le=digest(edges),
                             kinetic_output_edges_sha256_float64_le=digest(ke)))
            with (folder/"multiple_loss_spectrum.csv").open(encoding="ascii", newline="") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            self.assertEqual(reader.fieldnames[:7], ["region", "binding_edge_low_eV", "binding_edge_high_eV",
                             "kinetic_edge_low_eV", "kinetic_edge_high_eV", "zero_loss", "zero_loss_error"])
            ends = ["low_binding", "high_binding"] if name=="up" else ["high_binding", "low_binding"]
            self.assertEqual([r["region"] for r in rows], [ends[0], "0", "1", "2", "3", ends[1]])
            outside = {r["region"]: list(r.values())[1:5] for r in (rows[0],rows[-1])}
            self.assertEqual(outside, dict(low_binding=["","-2","34",""], high_binding=["2","","","30"]))
            for i, row in enumerate(rows[1:-1]):
                lo, hi = sorted(edges[i:i+2])
                self.assertEqual(list(map(float,list(row.values())[1:5])), [lo,hi,32.-hi,32.-lo])
            mass = [1., 3., 0., 0.]
            self.assertEqual([float(r["zero_loss"]) for r in rows[1:-1]], mass if name=="up" else mass[::-1])

    def test_legacy_manifest_and_capabilities(self):
        legacy = dict(schema="xpsloss-single-loss-spectrum-v0", axis="kinetic_energy_eV",
                      primary=dict(kind="delta", energies_eV=[32.], intensities=[1.]), output_edges_eV=[30.,32.,34.])
        status, reply = self.invoke(legacy)
        self.assertEqual(status, 0, reply)
        m = json.loads((self.root/"out"/"run.json").read_text(encoding="utf8"))
        self.assertNotIn('"energy_axis":', json.dumps(m))
        self.assertNotIn("measurement", m)
        self.assertEqual(set(m["spectrum_input"]), {"file", "sha256", "output_edges_sha256_float64_le"})
        self.assertEqual(set(m["results"]["outside"]["low"]), {"spectrum", "spectrum_error"})
        self.assertEqual(set(m["results"]["outside"]), {"low", "high"})
        cap = cli.capabilities()["multiple_scattering"]
        self.assertEqual(cap["energy_axes"], ["kinetic_energy_eV", "binding_energy_eV"])
        self.assertNotIn("binding_energy_axis", cap["unsupported_physics"])
        self.assertEqual(cap["measured_spectrum_schema"], "xpsloss-measured-spectrum-v0")
        # A nested discriminator in an invalid KE document must retain the
        # original KE error shape, even on the JSON decoder's failure path.
        malformed = '{"nested":{"schema":"xpsloss-measured-spectrum-v0"},' + json.dumps(legacy)[1:-1]
        code, err = self.invoke(malformed, out="invalid")
        self.assertEqual((code, err["reason"]), (2, "invalid_schema"))
        self.assertNotIn("detail", err)
