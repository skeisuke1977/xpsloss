"""Noninteractive CLI. stdout is JSON; errors are JSON on stderr."""

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy

from . import __version__
from .config import load_config
from .dielectric import (
    DIELECTRIC_MODEL_ID, DielectricConvergenceError, DielectricGapEdgeError,
    dielectric_diagnostics, dielectric_response,
)
from .model import (BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2, MODEL_ID, FIXED_DEPTH_UNSUPPORTED_PHYSICS,
                    UNSUPPORTED_PHYSICS, bulk_diimfp, finite_number, loss_function)
from .xps_fixed_depth import (
    EPSILON_MOMENTA, FIXED_DEPTH_MODEL_ID, TERM_NAMES, FixedDepthConvergenceError,
    effective_cross_section_fixed_depth,
)


class JSONParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


def capabilities():
    return {
        "package": "xpsloss", "version": __version__, "model_id": MODEL_ID,
        "status": "experimental",
        "supported": ["strict_json_validation", "drude_lindhard_elf", "bulk_diimfp",
                      "causal_dielectric_response", "xps_fixed_depth"],
        "dielectric_model_id": DIELECTRIC_MODEL_ID,
        "unsupported_physics": list(UNSUPPORTED_PHYSICS),
        "unsupported_physics_model_id": MODEL_ID,
        "xps_fixed_depth": {
            "model_id": FIXED_DEPTH_MODEL_ID,
            "scope": "fixed depth, normal emission from a gapless material",
            "supported": ["surface_losses", "vacuum_losses", "core_hole", "electron_hole_interference"],
            "unsupported_physics": list(FIXED_DEPTH_UNSUPPORTED_PHYSICS),
            "implements": "Based on Eqs. (13)-(17) of Simonsen, Yubero & Tougaard (1997) with the causal "
                          "Drude-Lindhard response; epsilon momentum chosen explicitly",
            "verification": {
                "internal_tests": "I1-I8 of docs/M2_FIRST_CASE.md",
                "independent_reviews": "Separate agents reviewed the equations and implementation. "
                                       "There has been no human peer review.",
                "published_figures": "A separate reviewer agent reproduced Figs. 3, 5 and 8 within "
                                     "figure-reading uncertainty using independent code with the 1997 paper's "
                                     "Eq. (20) dielectric function, outside this package's oscillator model.",
                "queels_xps_v2_1": (
                    "At 1000 eV and normal emission, differences for the tested dispersionless synthetic Al "
                    "input were 0.1-0.4% of peak total. Both momentum choices failed the preregistered "
                    "criteria for dispersive Al at 9.967 angstrom. In three synthetic follow-up runs "
                    "(5.0-99.8 eV, depth nearest 10 angstrom), all three components agreed within 0.3% "
                    "of peak reference total with k_total and a candidate using the undispersed oscillator "
                    "energy in the numerator of Re(1/eps). The published Eq. (32) of Yubero et al. (1996) "
                    "used here exceeded that threshold for each input. The comparison tested numerical "
                    "outputs only. See docs/QUEELS_DISPERSION_TEST.md."),
                "experiment": "not compared",
            },
        },
        "queels_xps_equivalence": False,
        "material_validation": (
            "Literature ELF inputs (Si, Au, Al) and their provenance are available in benchmarks/inputs "
            "in the repository and source distribution. The wheel omits these inputs. Agreement with "
            "measured spectra and general equivalence to QUEELS-XPS have not been established."),
    }


def emit(data, stream=sys.stdout):
    print(json.dumps(data, ensure_ascii=False, allow_nan=False, indent=2), file=stream)


def build_parser():
    p = JSONParser(prog="xpsloss", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True, parser_class=JSONParser)
    cap = sub.add_parser("capabilities")
    cap.add_argument("--json", action="store_true", help="output is always JSON")
    val = sub.add_parser("validate")
    val.add_argument("config")
    for name in ("bulk", "elf", "dielectric", "xps-fixed-depth"):
        child = sub.add_parser(name)
        child.add_argument("config")
        child.add_argument("--out", type=Path, required=True, help="new output directory; never overwritten")
        if name in ("elf", "dielectric"):
            child.add_argument("--q-inv-A", dest="q_inv_A", type=float, default=0.0,
                               help="wavevector in inverse angstrom; default 0 (optical limit)")
        if name == "xps-fixed-depth":
            child.add_argument("--depth-A", required=True, type=float, help="creation depth in angstrom")
            child.add_argument("--epsilon-momentum", required=True, choices=EPSILON_MOMENTA)
            child.add_argument("--hole-charge", type=float, default=1.0, help="charge in units of +e")
            child.add_argument("--theta-deg", type=float, default=0.0,
                               help="angle from the surface normal (only 0 degrees)")
    return p


def run(args):
    if args.command == "capabilities":
        emit(capabilities())
        return 0
    config = load_config(args.config)
    if args.command == "validate":
        emit({"status": "valid", "model_id": MODEL_ID, "input_sha256": config.input_sha256,
              "config": config.public_dict(), "scope": capabilities()})
        return 0
    if args.out.exists() or args.out.is_symlink():
        raise FileExistsError(f"output path already exists: {args.out}")
    start = time.perf_counter()
    grid = config.grid()
    model_id = MODEL_ID
    warnings = []
    if args.command == "bulk":
        values = np.array([bulk_diimfp(config.material, config.electron_energy_eV, float(t),
                           epsabs=config.epsabs, epsrel=config.epsrel, limit=config.limit) for t in grid])
        table = np.column_stack((grid, values))
        header = "loss_eV,diimfp_per_A_eV,quadrature_error_per_A_eV"
        output_name = "bulk.csv"
        units = {"loss_eV": "eV", "diimfp_per_A_eV": "angstrom^-1 eV^-1",
                 "quadrature_error_per_A_eV": "angstrom^-1 eV^-1"}
        diagnostics = {"max_absolute_quadrature_error_per_A_eV": float(values[:, 1].max()),
                       "quadrature_error_interpretation": "Numerical estimate only; not physical uncertainty."}
    elif args.command == "elf":
        q = finite_number(args.q_inv_A, "q_inv_A")
        if q < 0:
            raise ValueError("q_inv_A must be nonnegative")
        values = loss_function(config.material, grid, q)
        table = np.column_stack((grid, values))
        header, output_name = "loss_eV,elf_dimensionless", "elf.csv"
        units = {"loss_eV": "eV", "elf_dimensionless": "dimensionless", "q_inv_A": "angstrom^-1"}
        diagnostics = {"q_inv_A": q, "optical_limit": q == 0,
                       "sum_rule_enforced": False, "real_dielectric_part_computed": False}
    elif args.command == "xps-fixed-depth":
        model_id = FIXED_DEPTH_MODEL_ID
        result = effective_cross_section_fixed_depth(
            config.material, config.electron_energy_eV, grid, args.depth_A, theta_deg=args.theta_deg,
            hole_charge=args.hole_charge, epsilon_momentum=args.epsilon_momentum,
            epsabs=config.epsabs, epsrel=config.epsrel, limit=config.limit,
        )
        columns = ("loss_eV", "total", "no_hole", "hole", *TERM_NAMES, "quadrature_error")
        table = np.column_stack([getattr(result, name) for name in columns])
        header, output_name = ",".join(columns), "xps_fixed_depth.csv"
        units = dict.fromkeys(columns, "angstrom^-1 eV^-1")
        units.update(loss_eV="eV", depth_A="angstrom", path_length_A="angstrom",
                     theta_deg="degree from surface normal", hole_charge="elementary charge")
        diagnostics = {
            "depth_A": result.depth_A, "path_length_A": result.path_length_A,
            "theta_deg": result.theta_deg, "hole_charge": result.hole_charge,
            "epsilon_momentum": result.epsilon_momentum, "velocity": "nonrelativistic",
            "delta_terms": "pi*delta(omega_perp) parts of Eqs. (14)-(16) contribute nothing for loss > 0 "
                           "and are omitted; only loss > 0 is accepted",
            "normalization": "per unit path length in the solid, including vacuum terms",
            "component_definitions": {
                "no_hole": "term_13a + term_13b + term_15_T1 + term_15_T3 + term_16_e",
                "hole": "term_14 + term_15_T2 + term_16_h (unit positive charge, includes interference)",
                "total": "no_hole + hole_charge * hole",
            },
            "quadrature": "adaptive 16/32-point Gauss-Legendre panels in log(k_parallel)",
            "max_absolute_quadrature_error": float(result.quadrature_error.max()),
            "quadrature_error_interpretation": "Sum of estimated panel and term errors, "
                "with hole terms weighted by max(1, abs(hole_charge)).",
        }
    else:
        model_id = DIELECTRIC_MODEL_ID
        q = finite_number(args.q_inv_A, "q_inv_A")
        response = dielectric_response(
            config.material, grid, q, epsabs=config.epsabs, epsrel=config.epsrel, limit=config.limit,
        )
        physics = dielectric_diagnostics(
            config.material, q, loss_eV=grid,
            epsabs=config.epsabs, epsrel=config.epsrel, limit=config.limit,
        )
        inverse, epsilon = response.negative_inverse_epsilon, response.epsilon
        table = np.column_stack((grid, inverse.real, inverse.imag, epsilon.real, epsilon.imag,
                                 response.real_quadrature_error))
        columns = ("loss_eV", "negative_inverse_epsilon_real", "negative_inverse_epsilon_imag",
                   "epsilon_real", "epsilon_imag", "negative_inverse_epsilon_real_quadrature_error")
        header, output_name = ",".join(columns), "dielectric.csv"
        units = dict.fromkeys(columns, "dimensionless")
        units.update(loss_eV="eV", q_inv_A="angstrom^-1")
        nearest_edge = None
        if config.material.gap_eV > 0:
            nearest_edge = float(np.min(abs(grid - config.material.gap_eV)))
        diagnostics = {
            **asdict(physics),
            "q_inv_A": q, "optical_limit": q == 0, "sum_rule_enforced": False,
            "electron_energy_used": False,
            "nearest_gap_edge_distance_eV": nearest_edge,
            "fourier_time_dependence": "exp(-i*omega*t)",
            "high_energy_limit_negative_inverse_epsilon": -1.0,
            "high_energy_tail": "untruncated_oscillator_model",
            "real_part_method": "kramers_kronig_hard_gap" if config.material.gap_eV else "analytic",
            "max_absolute_real_quadrature_error": float(response.real_quadrature_error.max()),
            "quadrature_error_interpretation": "Estimate for Re[-1/epsilon] only.",
        }
        if physics.gap_edge_exclusion_eV:
            diagnostics["gap_edge_policy"] = "reject_grid"
            warnings.append("The hard gap gives a logarithmic edge singularity and can produce a real-axis pole in epsilon.")
        warnings.extend(physics.warnings)
    # Computation succeeds before creating output. A run.json marks a completed run.
    args.out.mkdir(parents=True, exist_ok=False)
    np.savetxt(args.out / output_name, table, delimiter=",", header=header,
               comments="", fmt="%.17g", encoding="utf-8")
    manifest = {
        "schema_version": 1, "status": "complete", "command": args.command,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - start,
        "input_sha256": config.input_sha256, "input_file": str(args.config),
        "package_version": __version__, "model_id": model_id,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "config": config.public_dict(), "point_count": len(grid), "units": units,
        "constants": {"bohr_radius_A": BOHR_RADIUS_A, "hbar2_over_2me_eV_A2": HBAR2_OVER_2ME_EV_A2},
        "diagnostics": diagnostics, "files": [output_name],
        "unsupported_physics": list(UNSUPPORTED_PHYSICS),
        "queels_xps_equivalence": False,
        "references": ["https://doi.org/10.3390/nano13020339",
                       "https://doi.org/10.1002/sia.7095",
                       "https://www.quases.com/products/queels-xps/"],
    }
    if args.command in ("dielectric", "xps-fixed-depth"):
        manifest["input_config_model_id"] = manifest["config"]["model"]
        manifest["warnings"] = warnings
    if args.command == "dielectric":
        manifest["references"].append("https://doi.org/10.1209/0295-5075/90/10003")
    elif args.command == "xps-fixed-depth":
        manifest["unsupported_physics"] = list(FIXED_DEPTH_UNSUPPORTED_PHYSICS)
        manifest["references"].append("https://doi.org/10.1103/PhysRevB.56.1612")
    (args.out / "run.json").write_text(json.dumps(manifest, ensure_ascii=False, allow_nan=False, indent=2) + "\n", encoding="utf-8")
    emit({"status": "complete", "model_id": model_id, "output_dir": str(args.out.resolve()),
          "files": [output_name, "run.json"], "point_count": len(grid),
          "queels_xps_equivalence": False})
    return 0


def main(argv=None):
    try:
        return run(build_parser().parse_args(argv))
    except DielectricGapEdgeError as exc:
        emit({
            "status": "error", "error_type": type(exc).__name__, "message": str(exc),
            "gap_eV": exc.gap_eV, "gap_edge_exclusion_eV": exc.exclusion_eV,
            "offending_loss_eV": exc.offending_loss_eV, "remedy": exc.remedy,
        }, sys.stderr)
        return 2
    except (DielectricConvergenceError, FixedDepthConvergenceError) as exc:
        emit({"status": "not_converged", "error_type": type(exc).__name__, "message": str(exc)}, sys.stderr)
        return 2
    except (ValueError, OSError, OverflowError, RuntimeError) as exc:
        emit({"status": "error", "error_type": type(exc).__name__, "message": str(exc)}, sys.stderr)
        return 2
