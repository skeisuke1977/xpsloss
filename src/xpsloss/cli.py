"""Noninteractive CLI. stdout is JSON; errors are JSON on stderr."""

import argparse
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
from .dielectric import DIELECTRIC_MODEL_ID, DielectricConvergenceError, dielectric_response
from .model import (BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2, MODEL_ID,
                    UNSUPPORTED_PHYSICS, bulk_diimfp, finite_number, loss_function)


class JSONParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


def capabilities():
    return {
        "package": "xpsloss-starter", "version": __version__, "model_id": MODEL_ID,
        "status": "experimental_starter",
        "supported": ["strict_json_validation", "drude_lindhard_elf", "bulk_diimfp",
                      "causal_dielectric_response"],
        "dielectric_model_id": DIELECTRIC_MODEL_ID,
        "unsupported_physics": list(UNSUPPORTED_PHYSICS),
        "queels_xps_equivalence": False,
        "material_validation": "No real material benchmark is bundled or claimed.",
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
    val.add_argument("config", type=Path)
    for name in ("bulk", "elf", "dielectric"):
        child = sub.add_parser(name)
        child.add_argument("config", type=Path)
        child.add_argument("--out", type=Path, required=True, help="new output directory; never overwritten")
        if name in ("elf", "dielectric"):
            child.add_argument("--q-inv-A", dest="q_inv_A", type=float, default=0.0,
                               help="wavevector in inverse angstrom; default 0 (optical limit)")
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
    else:
        model_id = DIELECTRIC_MODEL_ID
        q = finite_number(args.q_inv_A, "q_inv_A")
        response = dielectric_response(
            config.material, grid, q, epsabs=config.epsabs, epsrel=config.epsrel, limit=config.limit,
        )
        inverse, epsilon = response.negative_inverse_epsilon, response.epsilon
        table = np.column_stack((grid, inverse.real, inverse.imag, epsilon.real, epsilon.imag,
                                 response.real_quadrature_error))
        columns = ("loss_eV", "negative_inverse_epsilon_real", "negative_inverse_epsilon_imag",
                   "epsilon_real", "epsilon_imag", "negative_inverse_epsilon_real_quadrature_error")
        header, output_name = ",".join(columns), "dielectric.csv"
        units = dict.fromkeys(columns, "dimensionless")
        units.update(loss_eV="eV", q_inv_A="angstrom^-1")
        diagnostics = {
            "q_inv_A": q, "optical_limit": q == 0, "sum_rule_enforced": False,
            "fourier_time_dependence": "exp(-i*omega*t)",
            "high_energy_limit_negative_inverse_epsilon": -1.0,
            "high_energy_tail": "untruncated_oscillator_model",
            "real_part_method": "kramers_kronig_hard_gap" if config.material.gap_eV else "analytic",
            "max_absolute_real_quadrature_error": float(response.real_quadrature_error.max()),
            "quadrature_error_interpretation": "Estimate for Re[-1/epsilon] only.",
        }
        if config.material.gap_eV > 0 and any(osc.strength_eV2 for osc in config.material.oscillators):
            warnings.append("The hard gap gives a logarithmic edge singularity and can produce a real-axis pole in epsilon.")
    # Computation succeeds before creating output. A run.json marks a completed run.
    args.out.mkdir(parents=True, exist_ok=False)
    np.savetxt(args.out / output_name, table, delimiter=",", header=header,
               comments="", fmt="%.17g", encoding="utf-8")
    manifest = {
        "schema_version": 1, "status": "complete", "command": args.command,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - start,
        "input_sha256": config.input_sha256, "input_file": str(args.config.resolve()),
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
    if args.command == "dielectric":
        manifest["warnings"] = warnings
        manifest["references"].append("https://doi.org/10.1209/0295-5075/90/10003")
    (args.out / "run.json").write_text(json.dumps(manifest, ensure_ascii=False, allow_nan=False, indent=2) + "\n", encoding="utf-8")
    emit({"status": "complete", "model_id": model_id, "output_dir": str(args.out.resolve()),
          "files": [output_name, "run.json"], "point_count": len(grid),
          "queels_xps_equivalence": False})
    return 0


def main(argv=None):
    try:
        return run(build_parser().parse_args(argv))
    except DielectricConvergenceError as exc:
        emit({"status": "not_converged", "error_type": type(exc).__name__, "message": str(exc)}, sys.stderr)
        return 2
    except (ValueError, OSError, OverflowError, RuntimeError) as exc:
        emit({"status": "error", "error_type": type(exc).__name__, "message": str(exc)}, sys.stderr)
        return 2
