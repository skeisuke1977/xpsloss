"""Fermi-referenced binding energies and their once-rounded kinetic coordinates."""

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction as Q

import numpy as np

from ._fold import _up
from ._order_types import MultipleLossInputError
from ._single_loss_types import BinnedPrimary, DeltaPrimary, _primary, vector
from .model import finite_number


@dataclass(frozen=True, kw_only=True)
class EnergyAxisReference:
    kind: str
    kinetic_energy_offset_eV: float | None = None
    binding_anchor_eV: float | None = None
    kinetic_anchor_eV: float | None = None
    line: str | None = None
    citation: str | None = None


@dataclass(frozen=True, kw_only=True)
class BindingEnergyAxis:
    photon_energy_eV: float
    reference: EnergyAxisReference
    primary_reference_binding_eV: float
    primary_reference_choice: str
    source: Mapping | str
    calibration: Mapping | str
    charge_correction: Mapping | str
    core_level: str
    emission_angle_deg: float | str
    acceptance_half_angle_deg: float | str


@dataclass(frozen=True)
class BindingDeltaPrimary:
    binding_energies_eV: np.ndarray
    intensities: np.ndarray
    intensity_unit: str = "arb. unit"
    description: str | None = None
    included_broadening: str | None = None
    truncation: str | None = None
    intrinsic_satellites: str | None = None

    def __post_init__(self):
        _primary(self, "binding_energies_eV", False, binding=True)


@dataclass(frozen=True)
class BindingBinnedPrimary:
    binding_edges_eV: np.ndarray
    intensities: np.ndarray
    intensity_unit: str = "arb. unit"
    description: str | None = None
    included_broadening: str | None = None
    truncation: str | None = None
    intrinsic_satellites: str | None = None

    def __post_init__(self):
        _primary(self, "binding_edges_eV", True, binding=True)


PRIMARY_METADATA = ("intensity_unit", "description", "included_broadening", "truncation", "intrinsic_satellites")
REFERENCE_FIELDS = {
    "photon_and_offset": ("kinetic_energy_offset_eV",),
    "paired_anchor": ("binding_anchor_eV", "kinetic_anchor_eV", "line", "citation"),
}
RECORD_FIELDS = {
    "source": {"name": str, "monochromatic": bool, "satellite_treatment": str},
    "calibration": {"line": str, "reference_binding_eV": float, "measured_binding_eV": float, "applied_shift_eV": float},
    "charge_correction": {"applied": bool, "method": str, "shift_eV": float},
    "core_level": str,
    "emission_angle_deg": float,
    "acceptance_half_angle_deg": float,
}
NR = "not_recorded"


def invalid(field, detail):
    return MultipleLossInputError("invalid_energy_axis", field=field, detail=detail)


def _value(value, field, kind, missing=None):
    if value is None:
        raise invalid(field, "missing")
    if isinstance(value, str) and value == NR:
        if missing is None:
            raise invalid(field, "not_recorded_not_allowed")
        missing.append(field)
        return NR
    if kind is str:
        if not isinstance(value, str):
            raise invalid(field, "type")
        if not value.strip():
            raise invalid(field, "empty")
        try:
            value.encode("utf-8")
        except UnicodeError as exc:
            raise invalid(field, "type") from exc
        return value
    if kind is bool:
        if not isinstance(value, (bool, np.bool_)):
            raise invalid(field, "type")
        return bool(value)
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise invalid(field, "type")
    try:
        return finite_number(value, field)
    except (ValueError, OverflowError) as exc:
        raise invalid(field, "not_finite") from exc


def _record_value(name, value, schema, missing):
    if not isinstance(schema, dict) or isinstance(value, str) and value == NR:
        return _value(value, name, schema, missing)
    if not isinstance(value, Mapping):
        raise invalid(name, "missing" if value is None else "type")
    if set(schema) != set(value):
        raise invalid(name, "missing" if schema.keys() - value.keys() else "type")
    record = {}
    for key, kind in schema.items():
        field = f"{name}.{key}"
        if name == "charge_correction" and key != "applied" and record["applied"] is not True:
            if value[key] is not None:
                raise invalid(field, "exclusive")
            record[key] = None
        else:
            record[key] = _value(value[key], field, kind, missing)
    return record


def validate_axis(axis):
    """Validate before constructing a KE primary. Copy the complete measurement record."""
    if not isinstance(axis, BindingEnergyAxis):
        raise invalid("energy_axis", "type")
    photon = _value(axis.photon_energy_eV, "photon_energy_eV", float)
    if photon <= 0:
        raise invalid("photon_energy_eV", "out_of_range")
    ref = axis.reference
    if not isinstance(ref, EnergyAxisReference):
        raise invalid("reference", "type")
    kind = _value(ref.kind, "reference.kind", str)
    if kind not in REFERENCE_FIELDS:
        raise invalid("reference.kind", "unknown_kind")
    allowed = REFERENCE_FIELDS[kind]
    if any(v is not None for k, v in vars(ref).items() if k != "kind" and k not in allowed):
        raise invalid("reference", "exclusive")
    missing = []
    values = {}
    for k in allowed:
        label = k in ("line", "citation")
        values[k] = _value(getattr(ref, k), "reference." + k, str if label else float, missing if label else None)
    ref = EnergyAxisReference(kind=kind, **values)
    binding = _value(axis.primary_reference_binding_eV, "primary_reference_binding_eV", float)
    choice = _value(axis.primary_reference_choice, "primary_reference_choice", str, missing)
    record = {k: _record_value(k, getattr(axis, k), schema, missing) for k, schema in RECORD_FIELDS.items()}
    angle = record["acceptance_half_angle_deg"]
    if angle != NR and angle < 0:
        raise invalid("acceptance_half_angle_deg", "out_of_range")
    return dict(axis="binding_energy_eV", binding_reference="fermi_level", reference_kind=kind,
                photon_energy_eV=photon, kinetic_energy_offset_eV=ref.kinetic_energy_offset_eV,
                binding_anchor_eV=ref.binding_anchor_eV, kinetic_anchor_eV=ref.kinetic_anchor_eV,
                anchor_line=ref.line, anchor_citation=ref.citation, primary_reference_binding_eV=binding,
                primary_reference_choice=choice, record=record, not_recorded=tuple(missing))


def _rounded(exact, field):
    try:
        return float(exact)
    except OverflowError as exc:
        raise invalid(field, "not_finite") from exc


def _direction(edges):
    return "increasing" if edges[0] < edges[-1] else "decreasing"


def convert(axis, primary, edges, energy=None):
    from ._binding_support import BindingDisplay

    info = validate_axis(axis)
    if not isinstance(primary, (BindingBinnedPrimary, BindingDeltaPrimary)):
        raise invalid("primary", "requires_binding_primary")
    edges = vector(edges, "invalid_output_edges", increasing=True, binding=True)
    photon = Q(info["photon_energy_eV"])
    offset = info["kinetic_energy_offset_eV"]
    constant = photon - Q(offset) if offset is not None else Q(info["binding_anchor_eV"]) + Q(info["kinetic_anchor_eV"])
    ref_ke = constant - Q(info["primary_reference_binding_eV"])
    info.update(conversion_constant_exact=f"{constant.numerator}/{constant.denominator}",
                conversion_constant_eV=_rounded(constant, "conversion_constant_eV"),
                implied_offset_eV=None if offset is not None else _rounded(photon - constant, "implied_offset_eV"),
                primary_reference_kinetic_eV=_rounded(ref_ke, "primary_reference_kinetic_eV"),
                electron_energy_minus_reference_eV=None)
    if energy is not None:
        try:
            e0 = finite_number(energy, "electron_energy_eV")
        except (ValueError, OverflowError):
            pass  # M4a validates its own argument after conversion and injectivity.
        else:
            info["electron_energy_minus_reference_eV"] = _rounded(Q(e0) - ref_ke, "electron_energy_minus_reference_eV")

    binned = isinstance(primary, BindingBinnedPrimary)
    field = "binding_edges_eV" if binned else "binding_energies_eV"
    pos = getattr(primary, field)
    image, roles = {}, {}
    largest = Q(0)
    for role, coordinates in (("primary." + field, pos), ("output_edges_eV", edges)):
        for index, x in enumerate(coordinates):
            if x in image:
                continue
            exact = constant - Q(x)
            y = _rounded(exact, role)
            if y < 0:
                raise invalid(role, "negative_kinetic_energy")
            roles[x] = f"{role}[{index}]"
            image[x] = y
            largest = max(largest, abs(Q(y) - exact))
    ordered = sorted(image)
    for x, y in zip(ordered, ordered[1:]):
        if image[x] <= image[y]:
            raise MultipleLossInputError("axis_conversion_not_injective", binding_eV=(x, y),
                                        kinetic_eV=(image[x], image[y]), roles=(roles[x], roles[y]))

    ke = np.array([image[x] for x in pos])
    intensities = primary.intensities
    meta = {k: getattr(primary, k) for k in PRIMARY_METADATA}
    if binned:
        if ke[0] > ke[-1]:
            ke, intensities = ke[::-1], intensities[::-1]
        kinetic = BinnedPrimary(ke, intensities, **meta)
    else:
        order = np.argsort(ke, kind="stable")
        kinetic = DeltaPrimary(ke[order], intensities[order], **meta)
    increasing = edges[0] < edges[-1]
    kinetic_edges = np.array([image[x] for x in (edges[::-1] if increasing else edges)])
    kinetic_edges.setflags(write=False)
    info.update(rounding_max_eV=_up(largest), input_primary=primary, input_binding_edges_eV=edges,
                kinetic_primary=kinetic, kinetic_output_edges_eV=kinetic_edges,
                direction=dict(output_edges=_direction(edges), primary=_direction(pos) if binned else None))

    # Duplicate delta positions retain the first input index, even when that line is dark.
    indices = {}
    support = set()
    for i, x in enumerate(pos):
        indices.setdefault(x, i)
    for i in np.flatnonzero(primary.intensities):
        support.add(pos[i])
        if binned:
            support.add(pos[i+1])
    points = [(x, image[x], indices[x]) for x in ordered if x in support]
    return BindingDisplay(info, points, image)
