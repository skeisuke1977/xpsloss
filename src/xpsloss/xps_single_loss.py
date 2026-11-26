"""Finite-window integration of Eq. (19), Simonsen et al., Phys. Rev. B 56, 1612.

The depth weight and the satellite multiplier use the same externally supplied
length in angstrom. Detector bins contain integrated intensity on the kinetic
energy axis, including separate low and high overflow regions.
"""

from fractions import Fraction
from importlib import resources
import json

import numpy as np

from . import _loss_panels
from ._loss_roots import candidates
from ._single_loss_types import (
    BinnedPrimary, COLUMNS, COMPONENTS, DeltaPrimary, LossWindowAreaResult, QUANTITIES,
    SingleLossConvergenceError, SingleLossInputError, SingleLossSpectrumResult,
)
from .model import Material, finite_number
from .xps_depth_averaged import DepthAverageConvergenceError, effective_cross_section_depth_averaged
from .xps_fixed_depth import EPSILON_MOMENTA, FIXED_ABS, FIXED_REL

SINGLE_LOSS_MODEL_ID = "simonsen1997_single_loss_spectrum_external_v0"
AREA_DEFAULTS = dict(epsrel=1e-6, area_epsabs=1e-10, panel_limit=400, kernel_epsabs=FIXED_ABS,
                     kernel_epsrel=FIXED_REL, kernel_limit=200)
SPECTRUM_DEFAULTS = dict(**AREA_DEFAULTS, spectrum_epsabs_fraction=1e-9, max_cells=2_000_000)


def _call_kernel(material, electron_energy_eV, loss_eV, **kw):
    return effective_cross_section_depth_averaged(material, electron_energy_eV, loss_eV, **kw)


def _load_comparison_table():
    path = resources.files("xpsloss").joinpath("data").joinpath("m4a_completed_comparisons.json")
    return json.loads(path.read_text(encoding="utf-8"))


def _primary_values(primary):
    binned = isinstance(primary, BinnedPrimary)
    field = "edges_eV" if binned else "energies_eV"
    return {"kind": "binned" if binned else "delta", field: getattr(primary, field).tolist(),
            "intensities": primary.intensities.tolist()}


def _comparison(conditions, edges, basis):
    table = _load_comparison_table()
    groups = {"areas": [*("window_area."+c for c in COLUMNS), "single_loss_weight",
                        "inverse_window_area_interval_A.lower", "inverse_window_area_interval_A.upper",
                        *("cutoff_sensitivity."+c for c in COLUMNS)]}
    entry = "loss_window_area"
    if edges is not None:
        entry = "single_loss_spectrum"
        names = ["outside.low", *(f"bin[{i}]" for i in range(len(edges)-1)), "outside.high"]
        measures = [*QUANTITIES, *(f"cutoff_sensitivity.{s}.{c}" for s in ("plus", "minus") for c in COMPONENTS)]
        groups["spectrum"] = [r+"."+q for r in names for q in measures]
    kernel = {k: v for k, v in conditions.items() if k not in ("primary", "output_edges_eV")}
    partial, match = False, None
    for record in table["records"]:
        if record["entry"] != entry:
            continue
        if record["conditions"] == conditions:
            match = record["key_status"]
            break
        partial |= all(record["conditions"].get(k) == v for k, v in kernel.items())
    result = dict(key_status={}, physical_validity="not_established", resolution_basis=basis,
                  status="success", version=table["version"], registration_sha256=table["registration_sha256"],
                  reference_sha256=table["reference_sha256"])
    for group, keys in groups.items():
        provided = {} if match is None else match.get(group, {})
        statuses = {k: provided.get(k, dict(status="no_reference", reason="Required key has no completed comparison."))
                    for k in keys}
        result["key_status"][group] = statuses
        if match is None and partial:
            result[group] = "kernel_conditions_checked"
        else:
            result[group] = "matched" if all(s["status"] == "pass" for s in statuses.values()) else "not_matched"
    return result


class _DepthSampler:
    def __init__(self, material, energy, window, kw, batch_size, diag):
        self.material, self.energy, self.window = material, energy, window
        self.kw, self.batch_size, self.diag = kw, batch_size, diag
        self.seen = set()
        self.warnings = {}

    def _warnings(self, result, grid, stage):
        for message in result.warnings:
            record = self.warnings.get(message)
            if record is None:
                record = dict(source="xps_depth_averaged", code="inherited", message=message, occurrences=[])
                self.warnings[message] = record
                self.diag["warnings"].append(record)
            points = grid
            if message.startswith("Near the hard gap edge at loss_eV="):
                points = [float(message.split("loss_eV=", 1)[1].split(",", 1)[0])]
            occurrences = record["occurrences"]
            if not occurrences or occurrences[-1]["stage"] != stage:
                occurrences.append(dict(stage=stage, loss_eV=[]))
            occurrences[-1]["loss_eV"].extend(map(float, points))

    def __call__(self, points, stage):
        a, b = self.window
        internal = points[(points != a) & (points != b)]
        size = self.batch_size or max(1, len(internal))
        calls = self.diag["kernel_calls"]
        calls["points_sent_per_stage"].append(0)
        values, errors = np.empty((11, len(points))), np.empty((11, len(points)))
        for start in range(0, max(1, len(internal)), size):
            grid = np.r_[a, internal[start:start+size], b]
            self.seen.update(grid)
            calls["distinct_points"] = len(self.seen)
            calls["calls"] += 1
            calls["batch_sizes"].append(len(grid))
            calls["points_sent_per_stage"][-1] += len(grid)
            calls["points_sent_total"] += len(grid)
            try:
                result = _call_kernel(self.material, self.energy, grid, **self.kw)
            except DepthAverageConvergenceError as exc:
                raise SingleLossConvergenceError(
                    "kernel_not_converged", node_kind="sentinel" if exc.loss_eV in self.window else "interior",
                    stage=stage, index=exc.index, column=exc.column, loss_eV=exc.loss_eV,
                    message_details=dict(estimate=exc.estimate, budget=exc.budget),
                    inner_exception=type(exc).__name__, inner_reason=exc.reason, inner=str(exc),
                    inner_attributes=vars(exc).copy(), diagnostics=self.diag) from exc
            self._warnings(result, grid, stage)
            indices = np.searchsorted(points, grid)
            valid = indices < len(points)
            valid[valid] &= points[indices[valid]] == grid[valid]
            for j, c in enumerate(COLUMNS):
                values[j, indices[valid]] = getattr(result, c)[valid]
                errors[j, indices[valid]] = result.component_errors[c][valid]
        for record in self.warnings.values():
            last = record["occurrences"][-1]
            last["loss_eV"] = sorted(set(last["loss_eV"]))
        return values, errors


def _warnings(result, primary):
    diag = result.diagnostics
    def add(code, message):
        diag["warnings"].append(dict(source="xps_single_loss", code=code, message=message))
    add("window_only", "The integral includes only losses inside the specified window.")
    add("weight_not_probability", "The single-loss weight is an area ratio and can exceed one.")
    add("kernel_fixed_at_E0", "The loss kernel uses the specified electron energy for every primary energy.")
    a, b = result.loss_window_eV
    if primary is not None:
        add("primary_assumed_without_intrinsic_loss", "The primary is assumed to exclude intrinsic loss satellites.")
        if isinstance(primary, BinnedPrimary) and a < min(np.diff(primary.edges_eV)):
            add("window_lower_below_primary_bin_width", "The window starts below the smallest primary bin width.")
        support = diag["primary_support"]
        if support["min_eV"] is not None:
            if not support["min_eV"] <= result.electron_energy_eV <= support["max_eV"]:
                add("E0_outside_primary_support", "The kernel energy lies outside the positive primary support.")
            if support["min_eV"] < b:
                add("primary_support_below_window_upper", "The primary support extends below the upper loss.")
    if b/result.electron_energy_eV > .1:
        add("window_over_E0_gt_0p1", "The largest loss exceeds 0.1 times the kernel energy.")
    if result.hole_charge != 1:
        add("hole_charge_not_unity", "The specified hole charge differs from unity.")
    if diag["resolution_basis"] == "computed_candidates_only":
        add("resolution_basis_computed_candidates_only", "Resolution uses computed candidates with limited coverage.")
    comparison = diag["comparison"]
    if any(comparison[k] != "matched" for k in comparison["key_status"]):
        add("outside_matched_scope", "The completed comparison table does not cover all requested outputs.")
    result.warnings = tuple(r["message"] for r in diag["warnings"])


def _run(material, energy, primary, edges, loss_window_eV, length, mode, theta, charge, settings, batch_size=None):
    if not isinstance(material, Material):
        raise ValueError("material must be a Material")
    if not isinstance(mode, str) or mode not in EPSILON_MOMENTA:
        raise ValueError("epsilon_momentum must be 'k_parallel' or 'k_total'")
    if finite_number(theta, "theta_deg") != 0:
        raise ValueError("only normal emission (theta_deg=0) is supported")
    energy = finite_number(energy, "electron_energy_eV")
    length = finite_number(length, "weight_length_A")
    charge = finite_number(charge, "hole_charge")
    if energy <= 0 or length <= 0:
        raise ValueError("weight_length_A and electron_energy_eV must be positive")
    window = _loss_panels.checked_window(loss_window_eV, energy)
    _loss_panels.validate_settings(settings)
    active = [o for o in material.oscillators if o.strength_eV2]
    if active and window[0] <= material.gap_eV:
        raise SingleLossInputError("window_not_above_gap", gap_eV=material.gap_eV, loss_window_eV=window,
                                   message_details=dict(loss_window_eV=window))
    diag = _loss_panels._diagnostics(settings)
    features = {}
    if active:
        diag.update(candidates(material, energy, mode, window))
        features.update(singularities=[(r["z_re"], r["z_im"]) for r in diag["candidates"]],
                        min_half_width_eV=Fraction(min(o.width_eV for o in active))/2)
    else:
        diag["resolution_basis_reason"] = "zero_response"
    kw = dict(weight_length_A=length, epsilon_momentum=mode, theta_deg=theta, hole_charge=charge,
              epsabs=settings["kernel_epsabs"], epsrel=settings["kernel_epsrel"], limit=settings["kernel_limit"])
    sampler = _DepthSampler(material, energy, window, kw, batch_size, diag)
    result, _ = _loss_panels._integrate(sampler, window, features, length, charge, primary, edges, settings,
                                        diag, kernel_errstate=np.geterr(), energy=energy, mode=mode,
                                        gap=material.gap_eV, zero=not active)
    conditions = dict(material=dict(gap_eV=material.gap_eV,
                                    oscillators=[[o.energy_eV, o.strength_eV2, o.width_eV, o.dispersion]
                                                 for o in material.oscillators]),
                      electron_energy_eV=energy, weight_length_A=length, theta_deg=theta, epsilon_momentum=mode,
                      hole_charge=charge, loss_window_eV=list(window), effective_settings=settings)
    if primary is not None:
        edges = result.output_edges_eV
        conditions.update(primary=_primary_values(primary), output_edges_eV=edges.tolist())
    diag["comparison"] = _comparison(conditions, edges, diag["resolution_basis"])
    _warnings(result, primary)
    return result


def loss_window_area(material, electron_energy_eV, *, loss_window_eV, weight_length_A, epsilon_momentum,
                     theta_deg=0., hole_charge=1., epsrel=1e-6, area_epsabs=1e-10, panel_limit=400,
                     kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200) -> LossWindowAreaResult:
    """Integrate eleven depth-averaged cross sections over a finite positive window (angstrom^-1)."""
    settings = dict(epsrel=epsrel, area_epsabs=area_epsabs, panel_limit=panel_limit,
                    kernel_epsabs=kernel_epsabs, kernel_epsrel=kernel_epsrel, kernel_limit=kernel_limit)
    return _run(material, electron_energy_eV, None, None, loss_window_eV, weight_length_A,
                epsilon_momentum, theta_deg, hole_charge, settings)


def _single_loss_spectrum(material, electron_energy_eV, primary, output_edges_eV, *, loss_window_eV,
                          weight_length_A, epsilon_momentum, theta_deg=0., hole_charge=1.,
                          kernel_batch_size=None, **numerics):
    if not isinstance(primary, (DeltaPrimary, BinnedPrimary)):
        raise SingleLossInputError("invalid_primary")
    return _run(material, electron_energy_eV, primary, output_edges_eV, loss_window_eV, weight_length_A,
                epsilon_momentum, theta_deg, hole_charge, {**SPECTRUM_DEFAULTS, **numerics}, kernel_batch_size)


def single_loss_spectrum(material, electron_energy_eV, primary, output_edges_eV, *, loss_window_eV,
                         weight_length_A, epsilon_momentum, theta_deg=0., hole_charge=1., epsrel=1e-6,
                         area_epsabs=1e-10, spectrum_epsabs_fraction=1e-9, panel_limit=400, max_cells=2_000_000,
                         kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200) -> SingleLossSpectrumResult:
    """Fold a delta or binned primary with the single-loss kernel at fixed electron energy."""
    return _single_loss_spectrum(
        material, electron_energy_eV, primary, output_edges_eV, loss_window_eV=loss_window_eV,
        weight_length_A=weight_length_A, epsilon_momentum=epsilon_momentum, theta_deg=theta_deg,
        hole_charge=hole_charge, epsrel=epsrel, area_epsabs=area_epsabs, spectrum_epsabs_fraction=spectrum_epsabs_fraction,
        panel_limit=panel_limit, max_cells=max_cells, kernel_epsabs=kernel_epsabs,
        kernel_epsrel=kernel_epsrel, kernel_limit=kernel_limit)

