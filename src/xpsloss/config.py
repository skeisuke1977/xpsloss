"""Strict JSON input: unknown keys, nonfinite values and duplicates are errors."""

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from .model import MODEL_ID, Material, Oscillator, finite_number


def keys(obj, required, optional=(), *, label):
    if not isinstance(obj, dict):
        raise ValueError(f"{label} must be a JSON object")
    missing = set(required) - obj.keys()
    unknown = obj.keys() - set(required) - set(optional)
    if missing:
        raise ValueError(f"{label}: missing keys {sorted(missing)}")
    if unknown:
        raise ValueError(f"{label}: unknown keys {sorted(unknown)}")


@dataclass(frozen=True)
class Config:
    material: Material
    electron_energy_eV: float
    loss_min_eV: float
    loss_max_eV: float
    loss_step_eV: float
    point_count: int
    epsabs: float
    epsrel: float
    limit: int
    input_sha256: str

    def grid(self):
        return np.linspace(self.loss_min_eV, self.loss_max_eV, self.point_count)

    def public_dict(self):
        return {
            "model": MODEL_ID,
            "material": asdict(self.material),
            "electron_energy_eV": self.electron_energy_eV,
            "loss_grid_eV": {"min": self.loss_min_eV, "max": self.loss_max_eV,
                             "step": self.loss_step_eV},
            "integration": {"epsabs": self.epsabs, "epsrel": self.epsrel, "limit": self.limit},
        }


def parse_config(obj, input_sha256=""):
    keys(obj, ["material", "electron_energy_eV", "loss_grid_eV"],
         ["model", "integration"], label="config")
    if obj.get("model", MODEL_ID) != MODEL_ID:
        raise ValueError(f"unsupported model; only {MODEL_ID!r} is implemented")
    mat = obj["material"]
    keys(mat, ["name", "gap_eV", "oscillators"], label="material")
    if not isinstance(mat["oscillators"], list) or not mat["oscillators"]:
        raise ValueError("material.oscillators must be a nonempty array")
    oscillators = []
    for i, osc in enumerate(mat["oscillators"]):
        keys(osc, ["energy_eV", "strength_eV2", "width_eV", "dispersion"], label=f"oscillator {i}")
        oscillators.append(Oscillator(**osc))
    material = Material(mat["name"], mat["gap_eV"], tuple(oscillators))
    energy = finite_number(obj["electron_energy_eV"], "electron_energy_eV")
    grid = obj["loss_grid_eV"]
    keys(grid, ["min", "max", "step"], label="loss_grid_eV")
    lo, hi, step = (finite_number(grid[k], f"loss_grid_eV.{k}") for k in ("min", "max", "step"))
    if not 0 < lo <= hi < energy or step <= 0:
        raise ValueError("require 0 < grid.min <= grid.max < electron_energy_eV and step > 0")
    intervals = (hi - lo) / step
    if not math.isfinite(intervals) or intervals > 999_999:
        raise ValueError("loss grid must contain at most 1,000,000 points")
    count = round(intervals)
    if not math.isclose(intervals, count, rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError("(grid.max - grid.min) must be an integer multiple of grid.step")
    integ = obj.get("integration", {})
    keys(integ, [], ["epsabs", "epsrel", "limit"], label="integration")
    ea = finite_number(integ.get("epsabs", 1e-9), "integration.epsabs")
    er = finite_number(integ.get("epsrel", 1e-7), "integration.epsrel")
    limit = integ.get("limit", 200)
    if ea <= 0 or er <= 0:
        raise ValueError("integration epsabs and epsrel must be positive")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("integration.limit must be a positive integer")
    return Config(material, energy, lo, hi, step, count + 1, ea, er, limit, input_sha256)


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError(f"nonfinite JSON constant: {value}")


def load_config(path):
    raw = Path(path).read_bytes()
    obj = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_unique_object,
                     parse_constant=_reject_constant)
    return parse_config(obj, hashlib.sha256(raw).hexdigest())
