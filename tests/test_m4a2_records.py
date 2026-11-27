"""Cross-version snapshots and manifest-only replay for both public commands."""

from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from dataclasses import fields
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import io
import json
import math
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
from unittest.mock import patch
import zipfile

import numpy as np
import xpsloss
from xpsloss import cli
from xpsloss._single_loss_cli import json_value
import m4a2_states as m


def invoke(argv, folder):
    stdout, stderr = io.StringIO(), io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        code = cli.main(argv)
    out = Path(argv[argv.index("--out") + 1])
    files = {p.name: p.read_bytes() for p in out.iterdir()} if out.exists() else {}
    manifest = json.loads(files.pop("run.json")) if "run.json" in files else None
    message = json.loads(stdout.getvalue()) if stdout.getvalue() else None
    if manifest:
        assert datetime.fromisoformat(manifest["created_utc"]).utcoffset() == timezone.utc.utcoffset(None)
        assert math.isfinite(manifest["elapsed_seconds"]) and manifest["elapsed_seconds"] >= 0
        assert Path(message["output_dir"]).resolve() == out.resolve()
        message["output_dir"] = out.name
        for obj, field, digest in (
            (manifest, "input_file", "input_sha256"),
            (manifest["spectrum_input"], "file", "sha256"),
        ):
            if obj:
                path = Path(obj[field])
                assert path.parent == folder and hashlib.sha256(path.read_bytes()).hexdigest() == obj[digest]
                obj[field] = path.name
        del manifest["created_utc"], manifest["elapsed_seconds"]
    return dict(
        code=code,
        stdout=message,
        stderr=stderr.getvalue(),
        exists=out.exists(),
        manifest=manifest,
        files=files,
    )


def combination(n):
    cfg = m.TABLE["conventions"][str(n)]
    kw = dict(multiplier=cfg["multiplier"])
    if cfg["weight"] == "derived":
        kw.update(m.DERIVED, extra_inverse_length_per_A=cfg["c"], weight_length_seed_A=cfg.get("seed_A"))
    else:
        kw["weight_length_A"] = cfg["length_A"]
    if "multiplier_A" in cfg:
        kw["multiplier_A"] = cfg["multiplier_A"]
    return kw


def flags(kw, sources):
    args = []

    def flag(opt, value):
        if value is not None:
            args.extend(["--" + opt, str(value)])

    if kw["weight_length_A"] is None:
        args.append("--weight-length-inverse-13a-window-area")
        if kw.get("extra_inverse_length_per_A") or sources.get("c") is not None:
            flag("extra-inverse-length-per-A", kw["extra_inverse_length_per_A"])
        flag("extra-inverse-length-source", sources.get("c"))
        flag("weight-length-seed-A", kw.get("weight_length_seed_A"))
    else:
        flag("weight-length-A", kw["weight_length_A"])
        flag("weight-length-source", sources.get("weight"))
    if kw.get("multiplier") == "external":
        flag("multiplier-A", kw["multiplier_A"])
        flag("multiplier-source", sources.get("multiplier"))
    if kw.get("multiplier") == "inverse_window_area":
        args.append("--multiplier-inverse-window-area")
    return args


def arguments(folder, spec, kw, sources=None, partial=False):
    cfg = dict(
        material=m.old.material_dict("synthetic_wide"),
        electron_energy_eV=1000.0,
        loss_grid_eV=dict(min=1, max=100, step=1),
        integration=dict(epsabs=1e-12, epsrel=1e-9, limit=200),
    )
    config = folder / "input.json"
    config.write_text(json.dumps(cfg), encoding="ascii")
    argv = [
        "xps-single-loss" if spec else "xps-loss-window-area",
        str(config),
        "--out",
        str(folder / "output"),
        "--epsilon-momentum",
        "k_total",
        "--loss-window-eV",
        "5",
        "6",
    ]
    if spec:
        edges = m.old.make_edges("E_RB1").tolist()
        if partial:
            edges[0] -= 1
        obj = dict(
            schema="xpsloss-single-loss-spectrum-v0",
            axis="kinetic_energy_eV",
            primary=dict(kind="delta", energies_eV=[1000.0], intensities=[1.0]),
            output_edges_eV=edges,
        )
        p = folder / "spectrum.json"
        p.write_text(json.dumps(obj), encoding="ascii")
        argv.extend(["--spectrum", str(p)])
    return argv + flags(kw, sources or {})


def replay(manifest, folder):
    d = manifest["spectrum_definition"]
    w = d["weight"]
    kw = dict(weight_length_A=w["length_A"])
    if w.get("length_kind") == "inverse_13a_window_area":
        kw.update(
            m.DERIVED,
            extra_inverse_length_per_A=w["extra_inverse_length_per_A"],
            weight_length_seed_A=w["seed_A"],
        )
    mult = d.get("multiplier", dict(convention="weight_length"))
    kw["multiplier"] = mult["convention"]
    if kw["multiplier"] == "external":
        kw["multiplier_A"] = mult["value_A"]
    cfg = folder / "input.json"
    cfg.write_text(json.dumps(manifest["config"]), encoding="utf-8")
    argv = [
        manifest["command"],
        str(cfg),
        "--out",
        str(folder / "output"),
        "--epsilon-momentum",
        d["epsilon_momentum"],
        "--theta-deg",
        str(d["theta_deg"]),
        "--hole-charge",
        str(d["hole_charge"]),
        "--loss-window-eV",
        *map(str, d["loss_window_eV"]),
    ]
    # Preserve which of the effective settings the caller supplied on the CLI.
    for k, v in manifest["effective_settings"].items():
        if manifest["effective_settings_sources"][k] == "cli":
            argv += ["--" + k.replace("_", "-"), str(v)]
    if manifest["primary"] is not None:
        p = {
            k: v
            for k, v in manifest["primary"].items()
            if k not in ("support", "primary_abs_integral", "primary_integral")
        }
        spec = folder / "spectrum.json"
        spec.write_text(
            json.dumps(
                dict(
                    schema="xpsloss-single-loss-spectrum-v0",
                    axis=d["axis"],
                    primary=p,
                    output_edges_eV=manifest["output_edges_eV"],
                )
            ),
            encoding="utf-8",
        )
        argv += ["--spectrum", str(spec)]
    return invoke(
        argv
        + flags(
            kw,
            dict(
                weight=w.get("length_source"),
                c=w.get("extra_inverse_length_source"),
                multiplier=mult.get("source"),
            ),
        ),
        folder,
    )


def legacy_input_fixture(spec):
    """The replay input fields of the e4e024f short manifest (also available in an sdist)."""
    settings = dict(
        epsrel=1e-6,
        area_epsabs=1e-10,
        panel_limit=400,
        kernel_epsabs=1e-12,
        kernel_epsrel=1e-9,
        kernel_limit=200,
    )
    if spec:
        settings.update(spectrum_epsabs_fraction=1e-9, max_cells=2000000)
    return dict(
        schema_version=1,
        command="xps-single-loss" if spec else "xps-loss-window-area",
        config=dict(
            material=m.old.material_dict("synthetic_wide"),
            electron_energy_eV=1000.0,
            loss_grid_eV=dict(min=1, max=100, step=1),
            integration=dict(epsrel=1e-9, epsabs=1e-12, limit=200),
        ),
        spectrum_definition=dict(
            axis="kinetic_energy_eV",
            loss_window_eV=[5.0, 6.0],
            epsilon_momentum="k_total",
            theta_deg=0.0,
            hole_charge=1.0,
            weight=dict(length_A=20.0, length_kind="external_constant", length_source=None),
        ),
        effective_settings=settings,
        effective_settings_sources={k: "config" if k.startswith("kernel_") else "default" for k in settings},
        primary=dict(kind="delta", energies_eV=[1000.0], intensities=[1.0]) if spec else None,
        output_edges_eV=m.old.make_edges("E_RB1").tolist() if spec else None,
    )


def snapshot(row, current):
    args, kw = m.fixture(row["inputs"]["fixture"], "spectrum" in row["entry"])
    i = row["inputs"]
    if row["id"].startswith("V9-old-length") or row["id"] == "V9-order-old":
        kw["weight_length_A"] = old_length_value(row)
    if row["id"].startswith("V9-old-length") or row["id"] in ("V9-order-old", "V1-old-N25"):
        try:
            m.call(args, kw)
        except ValueError as exc:
            return dict(type=type(exc).__name__, message=str(exc), attributes=vars(exc))
        except m.api.SingleLossConvergenceError as exc:
            return dict(type=type(exc).__name__, message=str(exc), attributes=vars(exc))
        raise AssertionError("baseline failure did not occur")
    if current:
        if i.get("mode") == "equal_external":
            kw.update(multiplier="external", multiplier_A=kw["weight_length_A"])
        elif i.get("mode") == "explicit_weight_length":
            kw["multiplier"] = "weight_length"
        if "extra_inverse_length_per_A" in i:
            kw["extra_inverse_length_per_A"] = 0.0
    if row["entry"].startswith("cli"):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            return invoke(arguments(folder, "spectrum" in row["entry"], kw), folder)
    return vars(m.call(args, kw))


def child(destination, source, ids, current):
    assert Path(xpsloss.__file__).resolve().is_relative_to(Path(source).resolve())
    result = {}
    for rid in ids:
        row = m.ROWS[rid]
        if rid == "V10-caps":
            s = io.StringIO()
            with redirect_stdout(s):
                assert cli.main(["capabilities", "--json"]) == 0
            result[rid] = json.loads(s.getvalue())
        else:
            if current and os.environ.get("XPSLOSS_M4A2_CHILD_MUTANT"):
                from test_m4a2_mutations import fault

                with fault(os.environ["XPSLOSS_M4A2_CHILD_MUTANT"]):
                    result[rid] = snapshot(row, current)
            else:
                result[rid] = snapshot(row, current)
    Path(destination).write_bytes(pickle.dumps(result))


@lru_cache(None)
def pairs(ids):
    commit = m.TABLE["baseline"]
    assert (
        subprocess.check_output(["git", "rev-parse", "e4e024f^{commit}"], cwd=m.ROOT, text=True).strip()
        == commit
    )
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        raw = subprocess.check_output(["git", "archive", "--format=zip", commit, "src/xpsloss"], cwd=m.ROOT)
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            z.extractall(folder)
        out = []
        for current, src in ((False, folder / "src"), (True, m.ROOT / "src")):
            dest = folder / "data.bin"
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    str(dest),
                    str(src),
                    ",".join(ids),
                    str(int(current)),
                ],
                cwd=m.ROOT,
                env=dict(os.environ, PYTHONPATH=str(src)),
                check=True,
                capture_output=True,
            )
            out.append(pickle.loads(dest.read_bytes()))
    return out


def baseline(rid):
    short = (
        "V1-api_area-omit",
        "V1-api_spectrum-omit",
        "V1-cli_area-omit",
        "V1-cli_spectrum-omit",
        "V10-caps",
    )
    ids = (
        tuple(
            dict.fromkeys(
                [
                    *short,
                    *(
                        r["id"]
                        for r in m.profile()
                        if r["rule"] == "paired_baseline" and r["inputs"]["fixture"] in ("short",)
                    ),
                ]
            )
        )
        if rid in short or (m.ROWS[rid]["inputs"]["fixture"] == "short")
        else (rid,)
    )
    if rid not in ids or os.environ.get("XPSLOSS_M4A2_CHILD_MUTANT"):
        ids = (rid,)
    a, b = pairs(ids)
    return a[rid], b[rid]


def old_length_value(row):
    kw = {
        "none": None,
        "text_number": "2",
        "bool": True,
        "np_bool": np.bool_(True),
        "nan": math.nan,
        "inf": math.inf,
        "zero": 0,
        "negative": -1,
    }
    return kw.get(row["inputs"].get("value_token"), -1)


def old_length_message(row):
    if (m.ROOT / ".git").exists():
        ids = tuple(r for r in m.ROWS if r.startswith("V9-old-length") or r == "V9-order-old")
        before, _ = pairs(ids)
        failure = before[row["id"]]
        assert failure["type"] == "ValueError" and not failure["attributes"]
        return failure["message"]
    value = old_length_value(row)
    if value is None or isinstance(value, (str, bool, np.bool_)):
        return "weight_length_A must be a finite real number"
    if not math.isfinite(value):
        return "weight_length_A must be finite"
    return "weight_length_A and electron_energy_eV must be positive"


@lru_cache(None)
def base_version():
    """The __version__ of the e4e024f baseline, read from that commit."""
    source = subprocess.check_output(
        ["git", "show", m.TABLE["baseline"] + ":src/xpsloss/__init__.py"], cwd=m.ROOT, text=True
    )
    found = [line for line in source.splitlines() if line.startswith("__version__ = ")]
    assert len(found) == 1, found
    return found[0].split("=", 1)[1].strip().strip('"')


def own_versions(before, after):
    """Post-freeze change (owner approval 2026-11-27, docs/VALIDATION.md): each side reports its own version.

    The frozen pre-registration (section 3.1) compares the package version with the baseline, which held only
    while both sides were 0.3.0.dev0. A release changes this string without changing behaviour.
    """
    return before == base_version() and after == xpsloss.__version__


def legacy_items(log, before, after, prefix=""):
    if isinstance(before, dict):
        additions = set()
        if "window_area_budget" in before:
            additions = m.NEW
        elif prefix == "manifest.spectrum_definition":
            additions = {"multiplier"}
        elif prefix == "manifest.units":
            additions = {
                "multiplier_A",
                "multiplier_interval_A",
                "weight_length_interval_A",
                "seed_A",
                "area13",
                "error13",
                "extra_inverse_length_per_A",
            }
        log.check((prefix or "legacy") + ".keys", set(after) - additions == set(before) - additions)
        for k, v in before.items():
            path = prefix + "." + k if prefix else k
            if k == "diagnostics" or k in ("comparison", "effective_settings", "effective_settings_sources"):
                log.equal(path, after[k], v)
            else:
                legacy_items(log, v, after[k], path)
    elif prefix == "manifest.package_version":
        log.check(prefix, own_versions(before, after))
    else:
        log.equal(prefix, after, before)


def matched_table(spec, absent):
    t = deepcopy(m.old.comparison_table("test_table_C02"))
    rec = t["records"][0]
    cond = rec["conditions"]
    cond["loss_window_eV"] = [5.0, 6.0]
    if absent:
        cond["epsilon_momentum"] = "k_parallel"
    if not spec:
        rec["entry"] = "loss_window_area"
        del cond["primary"], cond["output_edges_eV"]
        cond["effective_settings"] = m.api.AREA_DEFAULTS.copy()
    return t


def convention(row, log):
    from test_m4a2_driver import Simulation

    n = row["inputs"]["combination"]
    kw = combination(n)
    spec = "spectrum" in row["entry"]
    match = row["inputs"].get("comparison_fixture")
    sources = dict(
        weight="weight origin" if n <= 4 else None,
        multiplier="multiplier origin" if n in (2, 3, 6) else None,
        c="extra rate origin" if n in (6, 7) else None,
    )
    if row["id"] == "V10-coincidence":
        sim = Simulation(row)
        args, base = m.fixture()
        with sim.patches():
            r = m.call(args, base | kw)
        log.equal("multiplier_A", r.multiplier_A, 20.0)
        log.equal("multiplier_interval_A", r.multiplier_interval_A, m.reciprocal(1 / 20, 0.0)[1])
        log.check("model_id", r.model_id == m.api.SINGLE_LOSS_CONVENTIONS_MODEL_ID)
        log.check("comparison", r.diagnostics["comparison"]["areas"] == "not_matched")
        codes = m.old.warning_codes(r.diagnostics)
        log.check("warnings", "W-S2" in codes and "W-len" not in codes)
        return
    args, base = m.fixture("short", spec)
    if match == "partial":
        args[-1][0] -= 1
    with patch.object(m.api, "_load_comparison_table", return_value=matched_table(spec, match == "absent")):
        r = m.call(args, base | kw)
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            cliout = invoke(arguments(folder, spec, kw, sources, match == "partial"), folder)
        assert cliout["code"] == 0, cliout["stderr"]
        manifest = cliout["manifest"]
        with tempfile.TemporaryDirectory() as tmp:
            again = replay(manifest, Path(tmp))
    log.check("exit", again["code"] == 0 and cliout["stderr"] == "")
    log.equal("replay.csv", again["files"], cliout["files"])
    for k in (
        "results",
        "diagnostics",
        "warnings",
        "units",
        "spectrum_definition",
        "effective_settings",
        "effective_settings_sources",
        "comparison",
        "model_id",
        "runtime",
        "constants",
        "primary",
        "output_edges_eV",
    ):
        log.equal("replay." + k, again["manifest"][k], manifest[k])
    log.equal(
        "api",
        manifest["results"],
        json_value({k: v for k, v in vars(r).items() if k not in ("warnings", "diagnostics", "primary")}),
    )
    log.equal("api.diagnostics", manifest["diagnostics"], json_value(r.diagnostics))
    mid = m.api.SINGLE_LOSS_MODEL_ID if n <= 2 else m.api.SINGLE_LOSS_CONVENTIONS_MODEL_ID
    for k, val in (
        ("stdout.model_id", cliout["stdout"]["model_id"]),
        ("manifest.model_id", manifest["model_id"]),
        ("model_id", r.model_id),
    ):
        log.equal(k, val, mid)
    w = manifest["spectrum_definition"]["weight"]
    mult = manifest["spectrum_definition"]["multiplier"]
    log.equal("weight.source", w["length_source"], sources["weight"])
    log.equal("weight.length", w["length_A"], r.weight_length_A)
    log.equal("weight.kind", w["length_kind"], "external_constant" if n <= 4 else "inverse_13a_window_area")
    log.equal(
        "multiplier",
        mult,
        dict(
            convention=kw["multiplier"],
            value_A=r.multiplier_A,
            interval_A=json_value(r.multiplier_interval_A),
            source=sources["multiplier"],
        ),
    )
    if n >= 5:
        log.equal("seed_A", w["seed_A"], 10.0 if n == 6 else 20.0)
        log.equal("c.source", w["extra_inverse_length_source"], sources["c"])
        log.check("seed.settings", set(w["effective_settings"]) == set(m.api.AREA_DEFAULTS))
        log.check("seed.c", w["extra_inverse_length_per_A"] == kw["extra_inverse_length_per_A"])
        center, interval = m.reciprocal(w["area13"], w["error13"], w["extra_inverse_length_per_A"])
        log.equal("weight_interval", r.weight_length_interval_A, interval)
        log.equal("weight_center", r.weight_length_A, center)
    else:
        log.check(
            "external_weight",
            r.weight_length_interval_A is None and set(w) == {"length_A", "length_kind", "length_source"},
        )
    interval = (
        r.weight_length_interval_A if n == 5 else r.inverse_window_area_interval_A if n in (4, 7) else None
    )
    log.equal("multiplier_interval_A", r.multiplier_interval_A, interval)
    expected = set()
    if n > 2:
        expected.add("W-sens")
        if n >= 5:
            expected.add("W-13a")
        if n in (6, 7):
            expected.add("W-c")
        if n in (4, 7):
            expected.add("W-S2")
        if r.multiplier_A != r.weight_length_A:
            expected.add("W-len")
    codes = set(m.old.warning_codes(r.diagnostics))
    log.check(
        "warnings",
        {c for c in codes if c.startswith("W-")} == expected and manifest["warnings"] == list(r.warnings),
    )
    comp = manifest["comparison"]
    label = (
        "not_matched"
        if n > 2 or match == "absent"
        else "kernel_conditions_checked"
        if match == "partial"
        else "matched"
    )
    for g, keys in comp["key_status"].items():
        log.equal("comparison." + g, comp[g], label)
        if n > 2:
            log.check(g + ".keys", all(x["status"] == "no_reference" for x in keys.values()))
    log.equal("physical_validity", comp["physical_validity"], "not_established")
    for k in ("weight_length_A", "seed_A", "multiplier_A", "weight_length_interval_A", "multiplier_interval_A"):
        log.equal("units." + k, manifest["units"][k], "angstrom")
    for k in ("area13", "error13", "extra_inverse_length_per_A"):
        log.equal("units." + k, manifest["units"][k], "angstrom^-1")
    log.check("settings", len(manifest["effective_settings"]) == (8 if spec else 6))


def usage(row, log):
    spec = row["entry"] == "cli_spectrum"
    case = row["inputs"]["flag_case"]
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        argv = arguments(folder, spec, dict(weight_length_A=20))
        del argv[-2:]
        ext = ["--weight-length-A", "20"]
        derived = ["--weight-length-inverse-13a-window-area"]
        derived_only = "extra attenuation and seed options require --weight-length-inverse-13a-window-area"
        external_m_only = "--multiplier-source requires --multiplier-A"
        flags, message = dict(
            no_weight=([], "one of the arguments --weight-length-A --weight-length-inverse-13a-window-area is required"),
            both_weights=(ext + derived, "argument --weight-length-inverse-13a-window-area: not allowed with argument --weight-length-A"),
            derived_source=(derived + ["--weight-length-source", "text"], "--weight-length-source requires --weight-length-A"),
            external_c_zero=(ext + ["--extra-inverse-length-per-A", "0"], derived_only),
            external_c_source=(ext + ["--extra-inverse-length-source", "text"], derived_only),
            external_seed=(ext + ["--weight-length-seed-A", "20"], derived_only),
            missing_c_with_source=(derived + ["--extra-inverse-length-source", "text"],
                                   "--extra-inverse-length-source requires --extra-inverse-length-per-A"),
            both_m=(ext + ["--multiplier-A", "7", "--multiplier-inverse-window-area"],
                    "argument --multiplier-inverse-window-area: not allowed with argument --multiplier-A"),
            missing_m_with_source=(ext + ["--multiplier-source", "text"], external_m_only),
            s2_with_source=(ext + ["--multiplier-inverse-window-area", "--multiplier-source", "text"], external_m_only),
        )[case]
        with patch.object(m.api, "_call_kernel", wraps=m.api._call_kernel) as spy:
            out = invoke(argv + flags, folder)
    d = json.loads(out["stderr"])
    log.check(
        "stderr",
        d["status"] == "error"
        and d["error_type"] == "ValueError"
        and d["message"] == message
        and out["stderr"].isascii(),
        actual=d["message"], expected=message,
    )
    log.check("calls", spy.call_count == 0)
    log.check("stdout", out["stdout"] is None)
    log.check("exit", out["code"] == 2)
    log.check("outdir", not out["exists"])


def records(row, log):
    rid, rule = row["id"], row["rule"]
    if rule == "usage_error":
        usage(row, log)
    elif rid.startswith("V10-stderr"):
        from test_m4a2_driver import Simulation

        sim = Simulation(row)
        with tempfile.TemporaryDirectory() as tmp, sim.patches():
            folder = Path(tmp)
            out = invoke(arguments(folder, True, combination(7)), folder)
        d = json.loads(out["stderr"])
        for k in ("phase", "stage", "reason"):
            log.equal(k, d[k], row["expect"][k])
        log.check("stderr", d["status"] == "not_converged" and d["error_type"] == row["expect"]["exception"])
        log.equal("inner_attributes", d["inner_attributes"], json_value(vars(sim.inner)))
        log.check("exit", out["code"] == 2 and out["stdout"] is None and not out["exists"])
    elif rule == "convention_record":
        convention(row, log)
    elif not (m.ROOT / ".git").exists() and rule == "paired_baseline":
        log.record("baseline", "legacy", "not_run", reason="sdist_without_git_history")
    elif rule == "paired_baseline":
        old, new = baseline(rid)
        if row["entry"].startswith("cli"):
            legacy_items(log, old, new)
            result = new["manifest"]["results"]
        else:
            log.check("fields", set(new) == set(old) | m.NEW)
            for k, v in old.items():
                log.equal(k, new[k], v)
            result = new
        log.check(
            "new_fields",
            result["model_id"] == m.api.SINGLE_LOSS_MODEL_ID
            and result["weight_length_kind"] == "external_constant"
            and result["multiplier_A"] == result["weight_length_A"]
            and result["multiplier_interval_A"] is result["weight_length_interval_A"] is None
            and result["extra_inverse_length_per_A"] == 0.0
            and result["multiplier_convention"]
            == ("external" if row["inputs"].get("mode") == "equal_external" else "weight_length"),
        )
    elif rule == "old_manifest":
        if (m.ROOT / ".git").exists():
            old = baseline("V1-" + row["entry"] + "-omit")[0]
            legacy = old["manifest"]
        else:
            legacy = legacy_input_fixture(row["entry"] == "cli_spectrum")
            old = snapshot(row, False)
        with tempfile.TemporaryDirectory() as tmp:
            out = replay(legacy, Path(tmp))
        log.equal("csv", out["files"], old["files"])
        legacy_items(log, old["manifest"]["results"], out["manifest"]["results"])
    elif rule == "old_constructor":
        name = row["inputs"]["cls"]
        r = m.call(*m.fixture("short", name == "SingleLossSpectrumResult"))
        order = [f.name for f in fields(r) if not f.kw_only]
        if (m.ROOT / ".git").exists():
            old = baseline("V1-api_" + ("spectrum" if hasattr(r, "primary") else "area") + "-omit")[0]
            log.equal("field_order", order, list(old))
        other = type(r)(*(getattr(r, k) for k in order))
        log.equal("old_args", vars(other), vars(r))
        log.check("new_fields", {f.name for f in fields(r) if f.kw_only} == m.NEW)
    elif rule == "capabilities":
        if (m.ROOT / ".git").exists():
            old, new = baseline(rid)
            log.check(
                "other_blocks",
                own_versions(old.get("version"), new.get("version"))
                and {k: v for k, v in new.items() if k not in ("xps_single_loss", "version")}
                == {k: v for k, v in old.items() if k not in ("xps_single_loss", "version")},
            )
        else:
            new = cli.capabilities()
        b = new["xps_single_loss"]
        log.check(
            "xps_single_loss",
            b["model_ids"] == [m.api.SINGLE_LOSS_MODEL_ID, m.api.SINGLE_LOSS_CONVENTIONS_MODEL_ID]
            and b["model_id"] == b["model_ids"][0]
            and b["multiplier_conventions"] == list(m.api.MULTIPLIERS)
            and b["weight_length_kinds"] == list(m.api.WEIGHT_KINDS),
        )
    else:
        raise AssertionError(rule)


if __name__ == "__main__":
    child(sys.argv[1], sys.argv[2], sys.argv[3].split(","), bool(int(sys.argv[4])))
