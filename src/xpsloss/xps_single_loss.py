"""Finite-window integration of Eq. (19), Simonsen et al., Phys. Rev. B 56, 1612.

The depth weight and the satellite multiplier are constant lengths in angstrom.
Their finite-window conventions can be selected separately. Detector bins contain integrated intensity on the kinetic
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
    MULTIPLIERS, WEIGHT_KINDS, SINGLE_LOSS_MODEL_ID as SINGLE_LOSS_MODEL_ID, SINGLE_LOSS_CONVENTIONS_MODEL_ID,
    SingleLossConvergenceError, SingleLossInputError, SingleLossSpectrumResult,
)
from .model import Material, finite_number
from .xps_depth_averaged import DepthAverageConvergenceError, effective_cross_section_depth_averaged
from .xps_fixed_depth import EPSILON_MOMENTA, FIXED_ABS, FIXED_REL

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
    if result.model_id == SINGLE_LOSS_CONVENTIONS_MODEL_ID:
        if result.weight_length_kind == "inverse_13a_window_area":
            add("W-13a", "This window-dependent length is not an IMFP. Comparison with trajectory weighting "
                "covers only the registered Al and Si conditions at normal emission and unit hole charge.")
        if result.extra_inverse_length_per_A > 0:
            add("W-c", "The extra attenuation is supplied by the caller and was checked only in the "
                "supplementary M3a-2 conditions.")
        if result.multiplier_A != result.weight_length_A:
            add("W-len", "Different weight and multiplier lengths depart from the constant-attenuation derivation.")
        if result.multiplier_convention == "inverse_window_area":
            add("W-S2", "The single-loss weight is one by convention; it does not predict the relative satellite area.")
        add("W-sens", "Cutoff sensitivity holds the reported weight length and multiplier fixed.")
    result.warnings = tuple(r["message"] for r in diag["warnings"])


def _number(value, argument, reason, *, allow_zero=False):
    try:
        value = finite_number(value, argument)
        if value < 0 or (value == 0 and not allow_zero):
            raise ValueError
    except (ValueError, OverflowError) as exc:
        raise SingleLossInputError(reason, message_details=dict(argument=argument)) from exc
    return value


def _seed_length(material, energy, window, features, seed_length, extra, kw, batch_size, diag):
    settings = diag["effective_settings"]
    sampler = _DepthSampler(material, energy, window, dict(kw, weight_length_A=seed_length), batch_size, diag)
    try:
        seed, _ = _loss_panels._integrate(sampler, window, features, seed_length, kw["hole_charge"],
            None, None, settings, diag, kernel_errstate=np.geterr(), energy=energy,
            mode=kw["epsilon_momentum"], gap=material.gap_eV, seed=True, extra=extra)
    except (SingleLossInputError, SingleLossConvergenceError) as exc:
        exc.phase = "weight_length_seed"
        raise
    return seed


def _run(material, energy, primary, edges, loss_window_eV, length, mode, theta, charge, settings, batch_size=None,
         *, weight_length_kind="external_constant", extra_inverse_length_per_A=0., weight_length_seed_A=None,
         multiplier="weight_length", multiplier_A=None, on_final_panels=None, on_validated=None):
    if not isinstance(material, Material):
        raise ValueError("material must be a Material")
    if not isinstance(mode, str) or mode not in EPSILON_MOMENTA:
        raise ValueError("epsilon_momentum must be 'k_parallel' or 'k_total'")
    if finite_number(theta, "theta_deg") != 0:
        raise ValueError("only normal emission (theta_deg=0) is supported")
    energy = finite_number(energy, "electron_energy_eV")
    if not isinstance(weight_length_kind, str) or weight_length_kind not in WEIGHT_KINDS:
        raise SingleLossInputError("invalid_weight_length_kind", message_details=dict(argument="weight_length_kind"))
    derived = weight_length_kind == "inverse_13a_window_area"
    if derived:
        if length is not None:
            raise SingleLossInputError("invalid_weight_length", message_details=dict(argument="weight_length_A"))
    else:
        length = finite_number(length, "weight_length_A")
    charge = finite_number(charge, "hole_charge")
    if energy <= 0 or (not derived and length <= 0):
        raise ValueError("weight_length_A and electron_energy_eV must be positive")
    window = _loss_panels.checked_window(loss_window_eV, energy)
    _loss_panels.validate_settings(settings)
    active = [o for o in material.oscillators if o.strength_eV2]
    if active and window[0] <= material.gap_eV:
        raise SingleLossInputError("window_not_above_gap", gap_eV=material.gap_eV, loss_window_eV=window,
                                   message_details=dict(loss_window_eV=window))
    extra = _number(extra_inverse_length_per_A, "extra_inverse_length_per_A", "invalid_extra_inverse_length",
                    allow_zero=True)
    if extra and not derived:
        raise SingleLossInputError("extra_inverse_length_without_derived_weight",
                                  message_details=dict(argument="extra_inverse_length_per_A"))
    seed_length = 20.
    if weight_length_seed_A is not None:
        seed_length = _number(weight_length_seed_A, "weight_length_seed_A", "invalid_weight_length_seed")
        if not derived:
            raise SingleLossInputError("weight_length_seed_without_derived_weight",
                                      message_details=dict(argument="weight_length_seed_A"))
    if not isinstance(multiplier, str) or multiplier not in MULTIPLIERS:
        raise SingleLossInputError("invalid_multiplier", message_details=dict(argument="multiplier"))
    if multiplier == "external":
        multiplier_A = _number(multiplier_A, "multiplier_A", "invalid_multiplier_length")
    elif multiplier_A is not None:
        raise SingleLossInputError("invalid_multiplier_length", message_details=dict(argument="multiplier_A"))
    legacy = not derived and (multiplier == "weight_length" or (multiplier == "external" and multiplier_A == length))
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
    prepared = None
    weight_interval = None
    if not legacy or on_validated is not None:
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                prepared = _loss_panels._prepare(window, features, primary, edges, settings, diag, energy,
                                                 material.gap_eV, None if active else window)
                if not legacy:
                    _loss_panels.preflight(prepared, settings, diag)
        except ArithmeticError as exc:
            if legacy:
                raise
            raise SingleLossConvergenceError("nonfinite_output", diagnostics=diag) from exc
    if derived and not active:
        if not extra:
            raise SingleLossInputError("weight_length_undefined", message_details=dict(argument="weight_length_kind"))
        seed_diag = _loss_panels._diagnostics({k: settings[k] for k in AREA_DEFAULTS})
        center, interval = _loss_panels._denominator(dict(value=0., estimate=0.), extra,
            settings["epsrel"], "weight_length", None, seed_diag)
        seed = dict(length_A=center, interval_A=interval, area13=0., error13=0.,
                    seed_skipped="structural_zero", consistency="not_applicable")
    if not active and multiplier == "inverse_window_area":
        raise SingleLossInputError("multiplier_undefined", message_details=dict(argument="multiplier"))
    if on_validated is not None:
        on_validated(window, prepared[0], not active)
    if derived:
        if active:
            seed_diag = _loss_panels._diagnostics({k: settings[k] for k in AREA_DEFAULTS})
            seed = _seed_length(material, energy, window, features, seed_length, extra, kw, batch_size, seed_diag)
        length, weight_interval = seed["length_A"], seed["interval_A"]
        kw["weight_length_A"] = length
        diag["weight_length"] = dict(**seed, seed_A=seed_length, extra_inverse_length_per_A=extra,
            **{k: seed_diag[k] for k in ("kernel_calls", "effective_settings", "denominator", "final_panels",
                                          "stages", "warnings", "error_budget")})
    sampler = _DepthSampler(material, energy, window, kw, batch_size, diag)
    try:
        result, panels = _loss_panels._integrate(sampler, window, features, length, charge, primary, edges, settings,
            diag, kernel_errstate=np.geterr(), energy=energy, mode=mode, gap=material.gap_eV, zero=not active,
            prepared=prepared, multiplier=multiplier, multiplier_A=multiplier_A, weight_interval=weight_interval)
    except (SingleLossInputError, SingleLossConvergenceError) as exc:
        if not legacy:
            exc.phase = "production"
        raise
    result.weight_length_kind = weight_length_kind
    result.weight_length_interval_A = weight_interval
    result.extra_inverse_length_per_A = extra
    if derived and active:
        a, e = result.window_area["term_13a"], result.window_area_error["term_13a"]
        delta = abs(Fraction(seed["area13"])-Fraction(a))
        allowance = Fraction(seed["error13"])+Fraction(e)
        diag["weight_length"]["consistency"] = dict(beta_area13=a, beta_error13=e, delta=float(delta),
                                                    allowance=float(allowance), passed=delta <= allowance)
        if delta > allowance:
            raise SingleLossConvergenceError("weight_length_inconsistent", phase="production",
                stage=diag["stages"], diagnostics=diag, message_details=dict(argument="weight_length_kind"))
    conditions = dict(material=dict(gap_eV=material.gap_eV,
                                    oscillators=[[o.energy_eV, o.strength_eV2, o.width_eV, o.dispersion]
                                                 for o in material.oscillators]),
                      electron_energy_eV=energy, weight_length_A=length, theta_deg=theta, epsilon_momentum=mode,
                      hole_charge=charge, loss_window_eV=list(window), effective_settings=settings)
    if not legacy:
        result.model_id = SINGLE_LOSS_CONVENTIONS_MODEL_ID
        conditions.update(weight_length_kind=weight_length_kind, multiplier=multiplier, multiplier_A=result.multiplier_A)
        if derived:
            conditions["extra_inverse_length_per_A"] = extra
        diag.setdefault("multiplier", dict(convention=multiplier, value_A=result.multiplier_A, interval_A=None,
                                           rho=0.))
    if primary is not None:
        edges = result.output_edges_eV
        conditions.update(primary=_primary_values(primary), output_edges_eV=edges.tolist())
    diag["comparison"] = _comparison(conditions, edges, diag["resolution_basis"])
    _warnings(result, primary)
    if on_final_panels is not None:
        on_final_panels(result, panels)
    return result


def loss_window_area(material, electron_energy_eV, *, loss_window_eV, weight_length_A, epsilon_momentum,
                     theta_deg=0., hole_charge=1., epsrel=1e-6, area_epsabs=1e-10, panel_limit=400,
                     kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200,
                     weight_length_kind="external_constant", extra_inverse_length_per_A=0.,
                     weight_length_seed_A=None, multiplier="weight_length", multiplier_A=None) -> LossWindowAreaResult:
    """Integrate eleven depth-averaged cross sections over a finite positive window (angstrom^-1).

    weight_length_kind selects the constant depth weight. The default,
    "external_constant", requires a positive finite weight_length_A in angstrom.
    With "inverse_13a_window_area", pass weight_length_A=None. The length is
    1/(A_13a + c), using the term_13a area for this window and momentum choice.
    extra_inverse_length_per_A is the finite nonnegative attenuation c in
    inverse angstrom (default 0). A positive c requires the derived weight.
    weight_length_seed_A is a positive finite seed in angstrom for that
    calculation only. Its default None uses 20 angstrom for the derived weight.

    multiplier selects the length m in Eq. (19). The default "weight_length"
    uses the depth weight. "external" requires a positive finite multiplier_A
    in angstrom. For both other choices, multiplier_A must be None (default).
    "inverse_window_area" uses the reciprocal of the final total window area,
    with no added c in the denominator, so m times that area is one by
    convention. Each multiplier can be used with either weight. An external
    multiplier equal to the external weight follows the default numerical path.
    """
    settings = dict(epsrel=epsrel, area_epsabs=area_epsabs, panel_limit=panel_limit,
                    kernel_epsabs=kernel_epsabs, kernel_epsrel=kernel_epsrel, kernel_limit=kernel_limit)
    return _run(material, electron_energy_eV, None, None, loss_window_eV, weight_length_A,
                epsilon_momentum, theta_deg, hole_charge, settings, weight_length_kind=weight_length_kind,
                extra_inverse_length_per_A=extra_inverse_length_per_A, weight_length_seed_A=weight_length_seed_A,
                multiplier=multiplier, multiplier_A=multiplier_A)


def _single_loss_spectrum(material, electron_energy_eV, primary, output_edges_eV, *, loss_window_eV,
                          weight_length_A, epsilon_momentum, theta_deg=0., hole_charge=1.,
                          kernel_batch_size=None, weight_length_kind="external_constant", extra_inverse_length_per_A=0.,
                          weight_length_seed_A=None, multiplier="weight_length", multiplier_A=None, **numerics):
    if not isinstance(primary, (DeltaPrimary, BinnedPrimary)):
        raise SingleLossInputError("invalid_primary")
    return _run(material, electron_energy_eV, primary, output_edges_eV, loss_window_eV, weight_length_A,
                epsilon_momentum, theta_deg, hole_charge, {**SPECTRUM_DEFAULTS, **numerics}, kernel_batch_size,
                weight_length_kind=weight_length_kind, extra_inverse_length_per_A=extra_inverse_length_per_A,
                weight_length_seed_A=weight_length_seed_A, multiplier=multiplier, multiplier_A=multiplier_A)


def single_loss_spectrum(material, electron_energy_eV, primary, output_edges_eV, *, loss_window_eV,
                         weight_length_A, epsilon_momentum, theta_deg=0., hole_charge=1., epsrel=1e-6,
                         area_epsabs=1e-10, spectrum_epsabs_fraction=1e-9, panel_limit=400, max_cells=2_000_000,
                         kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200,
                         weight_length_kind="external_constant", extra_inverse_length_per_A=0., weight_length_seed_A=None,
                         multiplier="weight_length", multiplier_A=None) -> SingleLossSpectrumResult:
    """Fold a delta or binned primary with the single-loss kernel at fixed electron energy.

    weight_length_kind selects the constant depth weight. The default,
    "external_constant", requires a positive finite weight_length_A in angstrom.
    With "inverse_13a_window_area", pass weight_length_A=None. The length is
    1/(A_13a + c), using the term_13a area for this window and momentum choice.
    extra_inverse_length_per_A is the finite nonnegative attenuation c in
    inverse angstrom (default 0). A positive c requires the derived weight.
    weight_length_seed_A is a positive finite seed in angstrom for that
    calculation only. Its default None uses 20 angstrom for the derived weight.

    multiplier selects the length m in Eq. (19). The default "weight_length"
    uses the depth weight. "external" requires a positive finite multiplier_A
    in angstrom. For both other choices, multiplier_A must be None (default).
    "inverse_window_area" uses the reciprocal of the final total window area,
    with no added c in the denominator, so m times that area is one by
    convention. Each multiplier can be used with either weight. An external
    multiplier equal to the external weight follows the default numerical path.
    """
    return _single_loss_spectrum(
        material, electron_energy_eV, primary, output_edges_eV, loss_window_eV=loss_window_eV,
        weight_length_A=weight_length_A, epsilon_momentum=epsilon_momentum, theta_deg=theta_deg,
        hole_charge=hole_charge, epsrel=epsrel, area_epsabs=area_epsabs, spectrum_epsabs_fraction=spectrum_epsabs_fraction,
        panel_limit=panel_limit, max_cells=max_cells, kernel_epsabs=kernel_epsabs,
        kernel_epsrel=kernel_epsrel, kernel_limit=kernel_limit, weight_length_kind=weight_length_kind,
        extra_inverse_length_per_A=extra_inverse_length_per_A, weight_length_seed_A=weight_length_seed_A,
        multiplier=multiplier, multiplier_A=multiplier_A)

