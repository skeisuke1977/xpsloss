"""Check the fixed input records against their source bytes and the bulk schema."""

from dataclasses import asdict
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import unittest

from xpsloss.config import load_config, parse_config
from xpsloss.elfsource import parse_elf, read_elf


ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "benchmarks" / "inputs"
SOURCES = ROOT / "benchmarks" / "sources" / "zenodo6024064_v1.0"
ABSOLUTE = ("si_elfdb_v1.0", "au_elfdb_v1.0", "al_elfdb_v1.0")
MANUAL = "si_queels_manual_v2.1"


def read_record(stem):
    return json.loads((INPUTS / f"{stem}.json").read_text(encoding="utf-8"), parse_float=Decimal)


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

    def test_rejects_unknown_material_in_both_formats(self):
        for filename in ("ELF_Si.txt", "ELF.txt"):
            with self.subTest(filename=filename), self.assertRaisesRegex(ValueError, "unknown material 'CoO'"):
                read_elf(SOURCES / filename, material="CoO")

    def test_reads_all_blocks_in_a_well_formed_compilation(self):
        raw = (b"First [1, 2]\n@\n1 -2\n0\n1\n15 225 1.5 1\n"
               b"Second [3]\n@\n1 -2\n0.5\n2\n20 300 2 0.1\n"
               b"References\n[1] Citation.\n")
        first, second = parse_elf(raw).materials
        self.assertEqual((first.name, first.reference_ids, first.last_line), ("First", (1, 2), 6))
        self.assertEqual((second.name, second.reference_ids, second.last_line), ("Second", (3,), 12))

    def test_compilation_reports_missing_ppy_marker(self):
        for name in (None, "PPY - polypyrrole"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "line 415 .*missing '@'"):
                read_elf(SOURCES / "ELF.txt", material=name)


class BenchmarkRecordTests(unittest.TestCase):
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
                self.assertTrue(all(value is None for value in record["source_run_conditions"].values()))

    def test_manual_relative_strengths_cannot_be_used_in_bulk(self):
        manual = json.loads((INPUTS / f"{MANUAL}.json").read_text(encoding="utf-8"))
        self.assertEqual(manual["strength_basis"], "relative")
        self.assertIsNone(manual["source_run_conditions"]["normalization"])
        self.assertIsNone(manual["bulk_usage"]["run_config"])
        self.assertFalse((INPUTS / f"{MANUAL}.bulk.json").exists())
        with self.assertRaises(ValueError):
            load_config(INPUTS / f"{MANUAL}.json")
        # Even projecting onto the material fields expected by bulk leaves
        # explicitly relative strengths, which strict config validation rejects.
        material = {key: manual["material"][key] for key in ("name", "gap_eV", "oscillators")}
        settings = manual["source_run_conditions"]
        candidate = {
            "material": material,
            "electron_energy_eV": settings["electron_energy_eV"],
            "loss_grid_eV": settings["loss_grid_eV"],
        }
        with self.assertRaisesRegex(ValueError, "oscillator 0: missing keys .*strength_eV2"):
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

    def test_provenance_round_trips_through_public_dict(self):
        config = load_config(INPUTS / f"{ABSOLUTE[0]}.bulk.json")
        public = config.public_dict()
        saved = json.loads(json.dumps(public, allow_nan=False))
        self.assertEqual(parse_config(saved).public_dict(), public)
        del saved["provenance"]["raw_file_sha256"]
        del public["provenance"]["raw_file_sha256"]
        self.assertEqual(parse_config(saved).public_dict(), public)

    def test_provenance_rejects_wrong_shape_and_values(self):
        valid = json.loads((INPUTS / f"{ABSOLUTE[0]}.bulk.json").read_text(encoding="utf-8"))
        provenance = valid["provenance"]
        for bad in (None, [], {}, {**provenance, "unexpected": "value"},
                    {**provenance, "source_id": " "}, {**provenance, "reference_input": 5},
                    {**provenance, "run_conditions_source": False},
                    {**provenance, "raw_file_sha256": "0" * 63},
                    {**provenance, "raw_file_sha256": "g" * 64},
                    {**provenance, "raw_file_sha256": None}):
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, "provenance"):
                parse_config({**valid, "provenance": bad})


if __name__ == "__main__":
    unittest.main()
