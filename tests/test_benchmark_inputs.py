"""Check the fixed input records against their source bytes and the bulk schema."""

from codecs import BOM_UTF8
from dataclasses import asdict
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from xpsloss import cli
from xpsloss.config import load_config, parse_config
from xpsloss.elfsource import parse_elf, read_elf


ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "benchmarks" / "inputs"
SOURCES = ROOT / "benchmarks" / "sources" / "zenodo6024064_v1.0"
ABSOLUTE = ("si_elfdb_v1.0", "au_elfdb_v1.0", "al_elfdb_v1.0")


def read_record(stem):
    return json.loads((INPUTS / f"{stem}.json").read_text(encoding="utf-8"), parse_float=Decimal)


def raw_gapless_rows(raw):
    """Read the published decimal tokens directly as exact rational numbers."""
    body = raw.decode("iso-8859-1").split("\nReferences", 1)[0]
    blocks = re.split(r"(?m)^([^\n]+?)[ \t]+\[[0-9,]+\][ \t\r]*$", body)
    for name, block in zip(blocks[1::2], blocks[2::2], strict=True):
        rows = [line.split() for line in block.splitlines() if line.strip() not in ("", "@")]
        index = Fraction(rows[2][0])
        moment = sum(Fraction(strength) / Fraction(energy)**2 for energy, strength, _, _ in rows[3:])
        yield name, (index, moment)


class ELFReaderTests(unittest.TestCase):
    def test_preserves_exact_decimals_and_original_bytes(self):
        raw = (b"\xe9 in the header\r\n@\r\n0.001 -2\r\n\r\n0\r\n999\r\n"
               b"2.89 8.58685901234567890123E-02 0.81 1\r\n")
        source = parse_elf(raw)
        entry, = source.materials
        self.assertEqual(source.raw_file_sha256, hashlib.sha256(raw).hexdigest())
        self.assertEqual(entry.marker_line, 2)
        self.assertEqual(entry.last_line, 7)
        self.assertEqual((entry.ac, entry.bc, entry.refractive_index),
                         (Decimal("0.001"), Decimal("-2"), Decimal("999")))
        self.assertEqual(entry.oscillators[0].strength_eV2, Decimal("8.58685901234567890123E-02"))

    def test_requires_marker_and_oscillators(self):
        for raw, error in [
            (b"comment\n1 -2\n0\n1\n15 225 1.5 1\n", "missing '@'"),
            (b"@\n1 -2\n0\n1\n", "at least one oscillator"),
        ]:
            with self.subTest(raw=raw), self.assertRaisesRegex(ValueError, error):
                parse_elf(raw)

    def test_rejects_wrong_column_counts(self):
        rows = [b"@", b"1 -2", b"0", b"1", b"15 225 1.5 1"]
        for index, invalid in [(1, b"1"), (2, b"0 1"), (3, b"1 2"), (4, b"15 225 1.5 1 2")]:
            with self.subTest(index=index), self.assertRaisesRegex(ValueError, "numeric columns"):
                malformed = rows[:]
                malformed[index] = invalid
                parse_elf(b"\n".join(malformed))

    def test_rejects_nonfinite_and_nonnumeric_values(self):
        for invalid in (b"NaN", b"-Infinity", b"+inf", b"bad", b"1_000"):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "finite decimal"):
                parse_elf(b"@\n1 -2\n0\n1\n15 " + invalid + b" 1.5 1\n")
        with self.assertRaisesRegex(ValueError, "line 5 .*nonfinite after float conversion"):
            parse_elf(b"@\n1 -2\n0\n1\n15 1e999999 1.5 1\n")

    def test_rejects_out_of_range_material_and_oscillator_values(self):
        rows = [b"@", b"1 -2", b"0", b"1", b"15 225 1.5 1"]
        for index, invalid, field in [
            (2, b"-1", "gap must be nonnegative"),
            (3, b"0", "index must be positive"),
            (3, b"-1", "index must be positive"),
            (4, b"0 225 1.5 1", "energy and width must be positive"),
            (4, b"-15 225 1.5 1", "energy and width must be positive"),
            (4, b"15 225 0 1", "energy and width must be positive"),
            (4, b"15 225 -1.5 1", "energy and width must be positive"),
            (4, b"15 -225 1.5 1", "strength must be nonnegative"),
        ]:
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, f"line {index + 1}.*{field}"):
                malformed = rows[:]
                malformed[index] = invalid
                parse_elf(b"\n".join(malformed))

    def test_cr_only_and_utf8_bom_keep_line_numbers_and_hash(self):
        body = "Data for oscillators that describe an ELF function of Café\r@\r1 -2\r0\r1\r15 225 1.5 1\r"
        for raw in (body.encode("iso-8859-1"), BOM_UTF8 + body.encode("utf-8")):
            with self.subTest(raw=raw):
                source = parse_elf(raw, material="Café")
                entry, = source.materials
                self.assertEqual((entry.name, entry.marker_line, entry.last_line), ("Café", 2, 6))
                self.assertEqual(source.raw_file_sha256, hashlib.sha256(raw).hexdigest())

    def test_bom_errors_identify_encoding_or_data(self):
        for raw, error in [
            ("@\r".encode("utf-16"), "unsupported UTF-16 or UTF-32 BOM"),
            (BOM_UTF8 + b"\xff", "invalid UTF-8"),
            (BOM_UTF8 + b"@\r1 -2\r0\r1\r15 225 1\r", "line 5 .*numeric columns"),
        ]:
            with self.subTest(raw=raw), self.assertRaisesRegex(ValueError, error):
                parse_elf(raw)

    def test_rejects_unknown_material_in_both_formats(self):
        for filename in ("ELF_Si.txt", "ELF.txt"):
            with self.subTest(filename=filename), self.assertRaisesRegex(ValueError, "unknown material 'CoO'"):
                read_elf(SOURCES / filename, material="CoO")

    def test_reads_all_blocks_in_a_well_formed_compilation(self):
        raw = (b"First[1, 2]\n@\n1 -2\n0\n1\n15 225 1.5 1\n"
               b"Second [3]\n@\n1 -2\n0.5\n2\n20 300 2 0.1\n"
               b"References\n[1] First citation.\n[2] Second citation.\n[3] Third citation.\n")
        first, second = parse_elf(raw).materials
        self.assertEqual((first.name, first.reference_ids, first.last_line), ("First", (1, 2), 6))
        self.assertEqual((second.name, second.reference_ids, second.last_line), ("Second", (3,), 12))

    def test_requires_defined_reference_ids_when_footer_is_present(self):
        raw = b"A[77]\n@\n1 -2\n0\n1\n15 225 1.5 1\n"
        self.assertEqual(parse_elf(raw, material="A").materials[0].reference_ids, (77,))
        with self.assertRaisesRegex(ValueError, r"line 1 \(A\): undefined reference ids \[77\]"):
            parse_elf(raw + b"References\n[1] Citation.\n", material="A")

    def test_compilation_reports_missing_ppy_marker(self):
        for name in (None, "PPY - polypyrrole"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "line 415 .*missing '@'"):
                read_elf(SOURCES / "ELF.txt", material=name)


class BenchmarkRecordTests(unittest.TestCase):
    def test_gapless_sum_rule_table_matches_all_raw_entries(self):
        raw = (SOURCES / "ELF.txt").read_bytes()
        published = dict(raw_gapless_rows(raw))
        self.assertEqual(len(published), 48)
        # Complete only the missing PPY separator for full-reader coverage.
        lines = raw.splitlines(keepends=True)
        self.assertEqual(lines[414].strip(), b"PPY - polypyrrole [20]")
        self.assertEqual(lines[415].split(), [b"0.001", b"-2"])
        lines.insert(415, b"@\r\n")
        parsed = {
            entry.name: (
                Fraction(entry.refractive_index),
                sum(Fraction(o.strength_eV2) / Fraction(o.energy_eV)**2 for o in entry.oscillators),
            )
            for entry in parse_elf(b"".join(lines)).materials
        }
        self.assertEqual(parsed, published)

        doc = (ROOT / "docs" / "BENCHMARK_INPUTS.md").read_text(encoding="utf-8")
        cells = [[col.strip() for col in line.strip("|").split("|")]
                 for line in doc.splitlines() if line.startswith("|")]
        table = [row for row in cells if len(row) == 9 and row[1].isdigit()]
        self.assertEqual([row[0] for row in table], list(published))

        finite_index, index_999 = {}, {}
        for name, _, _, index_cell, _, target_cell, _, moment_cell, delta_cell in table:
            n, s0 = published[name]
            target = 1 - n**-2
            delta = s0 - target
            with self.subTest(material=name):
                self.assertEqual(Fraction(index_cell), n)
                self.assertEqual(target_cell, f"{float(target):.6f}")
                self.assertEqual(moment_cell, f"{float(s0):.9f}")
                self.assertEqual(delta_cell, f"{float(delta):+.8e}")
            (index_999 if n == 999 else finite_index)[name] = abs(delta)

        self.assertEqual((len(finite_index), len(index_999)), (34, 14))
        close_finite = {name for name, diff in finite_index.items() if diff <= Fraction("1e-5")}
        self.assertEqual(close_finite, {
            "Al2O3", "CuO", "HfO2", "NiO", "PA - polyacetylene", "PE - polyethylene",
            "PMMA - polymethyl methacrylate", "PS - polystyrene", "TiO2", "ZrO2",
        })
        close_finite_relaxed = {name for name, diff in finite_index.items() if diff <= Fraction("4.1e-5")}
        self.assertEqual(close_finite_relaxed, close_finite | {"SiO2"})
        self.assertEqual({name for name, diff in index_999.items() if diff > Fraction("2.2e-7")}, {"Mo", "Pd"})
        close_total = sum(diff <= Fraction("4.1e-5") for diff in (finite_index | index_999).values())
        self.assertEqual(close_total, 23)
        self.assertGreater(finite_index["Si"], Fraction("0.04"))

    def test_records_match_source_bytes_and_parameters(self):
        for stem in ABSOLUTE:
            with self.subTest(stem=stem):
                record = read_record(stem)
                metadata, material = record["source"], record["material"]
                location = metadata["location"]
                path = ROOT / location["file"]
                source = read_elf(path, material=material["name"])
                entry, = source.materials
                self.assertEqual(metadata["raw_file_sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
                self.assertEqual(source.raw_file_sha256, metadata["raw_file_sha256"])
                self.assertEqual(entry.label_line, location["label_line"])
                self.assertEqual(entry.marker_line, location["marker_line"])
                self.assertEqual([entry.marker_line + 1, entry.last_line], location["data_lines_inclusive"])
                self.assertEqual(list(entry.reference_ids), location.get("original_reference_ids", []))
                for field in ("name", "gap_eV", "refractive_index", "ac", "bc"):
                    self.assertEqual(getattr(entry, field), material[field])
                self.assertEqual(list(map(asdict, entry.oscillators)), material["oscillators"])
                self.assertEqual(record["strength_basis"], "absolute_eV2")
                self.assertEqual(record["source_run_conditions"], dict.fromkeys((
                    "electron_energy_eV", "emission_angle", "depth_A", "loss_grid_eV",
                    "lower_loss_limit_eV", "normalization")))

    def test_relative_strengths_cannot_be_used_in_bulk(self):
        candidate = json.loads((ROOT / "examples/synthetic_bulk.json").read_text(encoding="utf-8"))
        oscillator = candidate["material"]["oscillators"][0]
        oscillator["strength_relative"] = oscillator.pop("strength_eV2")
        with self.assertRaisesRegex(ValueError, "oscillator 0: missing keys .*strength_eV2"):
            parse_config(candidate)

        oscillator["strength_eV2"] = oscillator.pop("strength_relative")
        parse_config(candidate)  # The renamed values satisfy the numerical schema.
        candidate["provenance"] = {
            "source_id": "synthetic",
            "reference_input": "examples/synthetic_bulk.json",
            "run_conditions_source": "synthetic test",
            "strength_basis": "relative",
        }
        with self.assertRaisesRegex(ValueError, "provenance.strength_basis must be 'absolute_eV2'"):
            parse_config(candidate)

    def test_absolute_bulk_configs_validate_and_match_records(self):
        for stem in ABSOLUTE:
            with self.subTest(stem=stem):
                config_path = INPUTS / f"{stem}.bulk.json"
                config = load_config(config_path)
                public = config.public_dict()
                raw_config = json.loads(config_path.read_text(encoding="utf-8"), parse_float=Decimal)
                record = read_record(stem)
                self.assertEqual(public["model"], "bulk_drude_lindhard_v0")
                self.assertEqual(record["bulk_usage"]["run_config"], config_path.name)
                material = record["material"]
                self.assertEqual(raw_config["material"],
                                 {key: material[key] for key in ("name", "gap_eV", "oscillators")})
                self.assertEqual(public["provenance"]["reference_input"], f"benchmarks/inputs/{stem}.json")
                self.assertEqual(public["provenance"]["source_id"], record["source"]["id"])
                self.assertEqual(public["provenance"]["raw_file_sha256"], record["source"]["raw_file_sha256"])
                self.assertEqual(public["provenance"]["strength_basis"], record["strength_basis"])
                self.assertIn("ignores ac, bc, and refractive index n", public["provenance"]["scope_note"])
                self.assertIn("not a QUEELS-XPS comparison", public["provenance"]["scope_note"])

    def test_provenance_round_trips_through_public_dict(self):
        config = load_config(INPUTS / f"{ABSOLUTE[0]}.bulk.json")
        public = config.public_dict()
        saved = json.loads(json.dumps(public, allow_nan=False))
        self.assertEqual(parse_config(saved).public_dict(), public)
        del saved["provenance"]["raw_file_sha256"]
        del public["provenance"]["raw_file_sha256"]
        saved["provenance"]["reference_input"] = "unavailable/source.json"
        public["provenance"]["reference_input"] = saved["provenance"]["reference_input"]
        self.assertEqual(parse_config(saved).public_dict(), public)

    def test_bulk_manifest_keeps_provenance_scope(self):
        config_path = INPUTS / "al_elfdb_v1.0.bulk.json"
        expected = load_config(config_path).public_dict()["provenance"]
        with TemporaryDirectory() as tmp, patch("xpsloss.cli.emit"):
            out = Path(tmp) / "al"
            args = cli.build_parser().parse_args(["bulk", str(config_path), "--out", str(out)])
            self.assertEqual(cli.run(args), 0)
            manifest = json.loads((out / "run.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["config"]["provenance"], expected)

    def test_provenance_rejects_wrong_shape_and_values(self):
        valid = json.loads((INPUTS / f"{ABSOLUTE[0]}.bulk.json").read_text(encoding="utf-8"))
        provenance = valid["provenance"]
        for bad in (None, [], {}, {**provenance, "unexpected": "value"},
                    {**provenance, "source_id": " "}, {**provenance, "reference_input": 5},
                    {**provenance, "run_conditions_source": False},
                    {key: value for key, value in provenance.items() if key != "strength_basis"},
                    {**provenance, "strength_basis": "relative"},
                    {**provenance, "scope_note": ""},
                    {**provenance, "scope_note": "lone\udfff"},
                    {**provenance, "raw_file_sha256": "0" * 63},
                    {**provenance, "raw_file_sha256": "g" * 64},
                    {**provenance, "raw_file_sha256": None}):
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, "provenance"):
                parse_config({**valid, "provenance": bad})


if __name__ == "__main__":
    unittest.main()
