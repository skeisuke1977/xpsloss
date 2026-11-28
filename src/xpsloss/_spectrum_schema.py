"""Compact spectrum grids and the measured-spectrum input schema (M4b-2 §3)."""

from dataclasses import fields
from fractions import Fraction
import json

import numpy as np

from ._binding_axis import (
    PRIMARY_METADATA, RECORD_FIELDS, REFERENCE_FIELDS, BindingBinnedPrimary, BindingDeltaPrimary,
    BindingEnergyAxis, EnergyAxisReference, invalid, validate_axis,
)
from ._order_types import MultipleLossInputError
from ._single_loss_types import SingleLossInputError, vector
from .config import keys
from .model import finite_number

SCHEMA = "xpsloss-measured-spectrum-v0"
AXIS = "binding_energy_eV"


def declares_measured_schema(raw):
    """Read the top-level discriminator of a rejected file, including escaped keys.

    Parsing one member at a time also works when a later member is malformed.
    A nested schema in a legacy file must never change its error record.
    """
    text = raw.decode("utf-8-sig", errors="replace")
    decoder = json.JSONDecoder()

    def after_space(i):
        while i < len(text) and text[i] in " \t\r\n":
            i += 1
        return i

    i = after_space(0)
    delimiter = "{"
    try:
        while text.startswith(delimiter, i):
            key, i = decoder.raw_decode(text, after_space(i + 1))
            i = after_space(i)
            if not isinstance(key, str) or not text.startswith(":", i):
                break
            value, i = decoder.raw_decode(text, after_space(i + 1))
            if key == "schema":
                return value == SCHEMA
            delimiter, i = ",", after_space(i)
    except ValueError:
        pass
    return False


def schema_error(message):
    detail = str(message).encode("ascii", "backslashreplace").decode("ascii")
    exc = MultipleLossInputError("invalid_schema", detail=detail)
    exc.detail = detail
    return exc


def _grid(spec, *, binding):
    keys(spec, ("min", "max", "step"), label="output_edges_eV")
    lo, hi, step = (finite_number(spec[k], k) for k in ("min", "max", "step"))
    if not lo < hi or (not binding and lo < 0) or step <= 0:
        raise ValueError("invalid output edge range")
    # Count the decimal grid before allocating the binary64 array.
    intervals = (Fraction(str(hi))-Fraction(str(lo)))/Fraction(str(step))
    if intervals.denominator != 1:
        raise ValueError("nonintegral output edge count")
    return lo, hi, int(intervals)


def edges(spec, max_cells, *, binding=False, compact=None):
    if isinstance(spec, dict):
        lo, hi, intervals = compact if compact is not None else _grid(spec, binding=binding)
        if intervals+2 > max_cells:
            raise SingleLossInputError("cell_limit_exceeded",
                                       message_details=dict(cells_at_least=intervals+2, max_cells=max_cells))
        spec = np.linspace(lo, hi, intervals+1)
    return vector(spec, "invalid_output_edges", increasing=True, binding=binding)


def parse_measured(spec, max_cells):
    try:
        keys(spec, ("schema", "axis", "measurement", "primary", "output_binding_edges_eV"), label="spectrum")
        if spec["schema"] != SCHEMA or spec["axis"] != AXIS:
            raise ValueError("unknown measured spectrum schema or axis")
        m = spec["measurement"]
        keys(m, tuple(f.name for f in fields(BindingEnergyAxis)), label="measurement")
        ref = m["reference"]
        forms = REFERENCE_FIELDS
        if isinstance(ref, dict):
            keys(ref, ("kind",), (*forms["photon_and_offset"], *forms["paired_anchor"]), label="reference")
            kind = ref["kind"]
            if isinstance(kind, str) and kind in forms:
                keys(ref, ("kind", *forms[kind]), (*forms["photon_and_offset"], *forms["paired_anchor"]),
                     label="reference")
        for name, shape in RECORD_FIELDS.items():
            if isinstance(shape, dict) and isinstance(m[name], dict):
                keys(m[name], shape, label=name)
        p = spec["primary"]
        if not isinstance(p, dict) or p.get("kind") not in ("binned", "delta"):
            raise ValueError("invalid primary kind")
        binned = p["kind"] == "binned"
        pos = "binding_edges_eV" if binned else "binding_energies_eV"
        keys(p, ("kind", pos, "intensities"), PRIMARY_METADATA, label="primary")
        edge_spec = spec["output_binding_edges_eV"]
        compact = _grid(edge_spec, binding=True) if isinstance(edge_spec, dict) else None
    except (ValueError, TypeError, OverflowError) as exc:
        raise schema_error(exc) from exc

    if isinstance(ref, dict):
        if isinstance(kind, str) and kind in forms and ref.keys() - {"kind", *forms[kind]}:
            raise invalid("reference", "exclusive")
        ref = EnergyAxisReference(**ref)
    axis = BindingEnergyAxis(**{**m, "reference": ref})
    validate_axis(axis)
    cls = BindingBinnedPrimary if binned else BindingDeltaPrimary
    primary = cls(**{k: v for k, v in p.items() if k != "kind"})
    return primary, edges(edge_spec, max_cells, binding=True, compact=compact), axis, m
