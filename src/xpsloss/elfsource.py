"""Read the QUEELS ELF data-set format without rescaling its values."""

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
from pathlib import Path
import re


_SINGLE_HEADER = re.compile(r"Data for oscillators that describe an ELF function of (.+)")
_MATERIAL_LABEL = re.compile(r"([^@\[\]]+?)\s+\[([0-9]+(?:\s*,\s*[0-9]+)*)\]")
_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


@dataclass(frozen=True)
class ELFOscillator:
    energy_eV: Decimal
    strength_eV2: Decimal
    width_eV: Decimal
    dispersion: Decimal


@dataclass(frozen=True)
class ELFMaterial:
    name: str | None
    ac: Decimal
    bc: Decimal
    gap_eV: Decimal
    refractive_index: Decimal
    oscillators: tuple[ELFOscillator, ...]
    label_line: int | None
    marker_line: int
    last_line: int
    reference_ids: tuple[int, ...] = ()


@dataclass(frozen=True)
class ELFSource:
    raw_file_sha256: str
    materials: tuple[ELFMaterial, ...]


def _numbers(row, columns, field):
    line, content = row
    tokens = content.split()
    if len(tokens) != columns:
        raise ValueError(f"line {line} ({field}): expected {columns} numeric columns")
    values = []
    for token in tokens:
        if not _NUMBER.fullmatch(token):
            raise ValueError(f"line {line} ({field}): expected a finite decimal, got {token!r}")
        try:
            values.append(Decimal(token))
        except InvalidOperation as exc:
            raise ValueError(f"line {line} ({field}): invalid decimal {token!r}") from exc
    return values


def _read_block(rows, name, label_line, reference_ids=()):
    if not rows or rows[0][1] != "@":
        raise ValueError(f"line {label_line} ({name}): missing '@' before material data")
    marker_line = rows[0][0]
    if len(rows) < 5:
        raise ValueError(f"line {marker_line}: expected ac/bc, gap, index and at least one oscillator")
    ac, bc = _numbers(rows[1], 2, "ac/bc")
    gap, = _numbers(rows[2], 1, "energy gap")
    index, = _numbers(rows[3], 1, "refractive index")
    oscillators = tuple(ELFOscillator(*_numbers(row, 4, "oscillator")) for row in rows[4:])
    return ELFMaterial(name, ac, bc, gap, index, oscillators, label_line,
                       marker_line, rows[-1][0], reference_ids)


def _read_compilation(rows, labels, material):
    if any(content == "@" for _, content in rows[:labels[0][0]]):
        raise ValueError("missing material label before '@'")
    names = [label[1] for label in labels]
    if len(set(names)) != len(names):
        raise ValueError("duplicate material name in ELF compilation")
    if material is not None and material not in names:
        raise ValueError(f"unknown material {material!r}")
    ends = [start for start, _, _ in labels[1:]] + [len(rows)]
    return tuple(
        _read_block(rows[start + 1:end], name, rows[start][0], reference_ids)
        for (start, name, reference_ids), end in zip(labels, ends)
        if material is None or material == name
    )


def parse_elf(raw: bytes, *, material: str | None = None) -> ELFSource:
    """Parse ELF bytes in ISO-8859-1 (also covers US-ASCII).

    Without a name, validate all material blocks. An exact material name
    limits validation to its block, but the hash always covers the whole file.
    Rows before the first marker in a single-material file are comments.
    Decimal values retain source precision, including any sentinel values.
    """
    rows = [(line, content.strip())
            for line, content in enumerate(raw.decode("iso-8859-1").split("\n"), 1)
            if content.strip()]
    marker = next((i for i, (_, row) in enumerate(rows) if row == "@"), len(rows))
    preamble = rows[:marker]
    single_headers = [(line, match[1]) for line, row in preamble
                      if (match := _SINGLE_HEADER.fullmatch(row))]
    if len(single_headers) > 1:
        raise ValueError("multiple single-material headers")
    if not single_headers:
        footer = next((i for i, (_, row) in enumerate(rows) if row == "References"), len(rows))
        labels = [(i, match[1], tuple(int(x) for x in match[2].split(",")))
                  for i, (_, row) in enumerate(rows[:footer])
                  if (match := _MATERIAL_LABEL.fullmatch(row))]
        if labels:
            entries = _read_compilation(rows[:footer], labels, material)
            return ELFSource(hashlib.sha256(raw).hexdigest(), entries)
    if marker == len(rows):
        raise ValueError("missing '@' before ELF data")
    label_line, name = single_headers[0] if single_headers else (None, None)
    if material is not None and material != name:
        raise ValueError(f"unknown material {material!r}")
    entry = _read_block(rows[marker:], name, label_line)
    return ELFSource(hashlib.sha256(raw).hexdigest(), (entry,))


def read_elf(path, *, material: str | None = None) -> ELFSource:
    """Read and hash a source file without altering its bytes."""
    return parse_elf(Path(path).read_bytes(), material=material)
