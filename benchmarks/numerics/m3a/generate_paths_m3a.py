"""Independent paths of M3a (docs/M3A_PREREG.md section 5, R2): path 1 (V8) and path 3 (V9).

Written by R2 without the depth-averaged implementation.  The rules are those of the
pre-registration, section 5; the code of the rules is in tests/m3a_paths.py (one copy,
also used by the tests for V2) and is restated here:

Path 1 (numerical depth quadrature of the public fixed-depth function only)
  * Sets P (V8): (material, E0, momentum, hole charge, lambda) with the losses of the
    registered runs: Al and Y0 at 1000 eV with the nine MS5 losses, Si with the eighteen
    MS11 losses, both momenta, Q_h in {0, 0.5, 1}, lambda in {5, 20, 80} A; Al at 4000 eV
    with GA[4], GA[77], GA[149], GA[498] (MS2), both momenta, Q_h = 1, 20 A.
  * One fixed-depth call per depth node and hole charge evaluates all losses of P
    (theta = 0, hole_charge = Q_h).  A failing call is replayed loss by loss; a failing loss
    is a failure of the path at that point (recorded, never dropped from the set P or from
    its scale).  Once a loss has failed at a node, the later nodes of the same chunk (fixed
    chunks of consecutive nodes) omit it from the call, which only saves the cost of
    repeating a known failure.
  * Depth rule t = a/lambda on [0, 40]: registered edges, geometric edges 40 2^-m down to
    1/(8 k_par+,max lambda), equal split for a phase increment b_eff lambda dt <= pi/4,
    Gauss-Legendre 10 nodes per piece, node weight alpha_i = w_i t_i e^-t_i.
  * Variations, one at a time: (A) 16 nodes; (B) pi/8 and geometric ratio sqrt 2; (C) r = 50
    (the base pieces plus [40, 50] cut by the phase rule).
  * Preliminary run (registered and geometric edges, no phase split, 10 nodes, fixed-depth
    accuracy (1e-9, 1e-7, 200), looser on failure) gives max_T |R_total| of each set.  The
    nodes of the final runs request epsabs_i = 0.9 f_i eta D_floor / alpha_i with
    f_i = 1/N (N nodes of that node set), eta = 1/9, D_floor = 1e-8 max_T |R_total,prelim|,
    epsrel = 1e-15 (absolute budgets) and limit 200.  The final D_j = 1e-7 |R_j| +
    1e-8 max_{T in P} |R_total| comes from the final base values, and alpha_i e_ij <=
    f_i eta D_j is checked afterwards for every node, column and loss (base nodes; the
    variant node sets are reported with their own N).
  * u_j = max(|Delta_A|, |Delta_B|, |Delta_C|) + sum_i alpha_i e_ij + rounding + tail,
    tail = |kappa_j| (1+r) e^-r + C_j e^-r / lambda with the envelope of design section 2.2
    from the public dielectric function (upper side of a 16/8-node quadrature in ln k).
    Gapless nodes publish one common quadrature_error, given to every column.
Path 3 (I3 type, frozen tests/xps_route_b_reference.py, k_parallel)
  * Al and Al0 together (same k_par nodes), 1000 eV, MS5 losses, lambda in {5, 20} A, at the
    base depth nodes of path 1 (the nodes of the whole nine-loss set P).  Unit-charge terms;
    total = no_hole + Q_h hole.  The frozen route treats each loss independently, so the tasks
    run loss by loss; the V9 scale max_9 |K_total,I3| needs the base values of all nine losses.
  * Variations, one at a time: (A), (B), (C) as above; (D) kpar_nodes(n_base=24); (E) the
    k_perp core and tail node arrays of the module twice as fine (Gauss-Legendre 40 instead
    of 20, Gauss-Laguerre 160 instead of 80, the 1/k tail 64 instead of 32).  (D) and (E)
    change the module in memory in a separate process; the file is never edited.
  * Each change must be at most 1/10 of the V9 width 1e-6 |K_j,I3| + 1e-6 max|K_total,I3|.

Usage (resumable; partial results in WORK):
  python benchmarks/numerics/m3a/generate_paths_m3a.py plan     --work WORK
  python benchmarks/numerics/m3a/generate_paths_m3a.py work     --work WORK [--kinds prelim,env,node,i3]
                                                                 [--sets SUBSTR,...] [--variants ...] [--max-tasks N]
  python benchmarks/numerics/m3a/generate_paths_m3a.py task     --work WORK TASK_ID
  python benchmarks/numerics/m3a/generate_paths_m3a.py status   --work WORK
  python benchmarks/numerics/m3a/generate_paths_m3a.py assemble --work WORK [--out PATH]
Several "work" processes may run at once: a task is claimed by creating its claim file.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tests"))
import m3a_paths as mp3  # noqa: E402

TABLE = HERE / "expected_states_m3a.json"
TABLE_SHA256 = "d87dad8f02a2fd211679229bd7c5205d908dbb2b3b241b30a27e94469f414ed5"
OUT = HERE / "paths_m3a.json"
SCHEMA = "xpsloss-m3a-paths-v1"
FINAL_SAFETY = 0.9
CHUNK_P1 = 1500
CHUNK_I3 = 60
I3_LAMBDAS = (5.0, 20.0)
I3_MATERIALS = ("Al", "Al0")
I3_LOSS_ORDER = (1.0, 15.6, 0.05, 5.0, 0.2, 11.0, 30.0, 60.0, 99.8)
I3_SCALE_LOSS = 0.05
PRELIM_FALLBACK = ((1e-9, 1e-7, 200), (1e-6, 1e-5, 200), (1e-3, 1e-3, 200))


# ----------------------------------------------------------------------------- inputs

def load_table():
    raw = TABLE.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != TABLE_SHA256:
        raise SystemExit(f"expected-state table SHA-256 {sha} differs from the registered {TABLE_SHA256}")
    return json.loads(raw.decode("utf-8"))


def _family_losses(table, family, material):
    out = []
    for run in table["runs"]:
        if run["family"] == family and run["material"] == material:
            for h in run["losses"]["values_hex"]:
                if h not in out:
                    out.append(h)
    return [float.fromhex(h) for h in out]


def path1_sets(table):
    """All V8 sets.  The id names material, E0, momentum, charge and length."""
    ga = 0.2 * np.arange(1, 500)
    t_gapless = _family_losses(table, "MS5", "Al")
    if _family_losses(table, "MS5", "Y0") != t_gapless:
        raise SystemExit("MS5 losses differ between Al and Y0")
    t_si = _family_losses(table, "MS11", "Si")
    t_4000 = [float(ga[i]) for i in (4, 77, 149, 498)]
    out = []
    for mat, losses, e0, lams, charges in (("Al", t_gapless, 1000.0, (5.0, 20.0, 80.0), (0.0, 0.5, 1.0)),
                                           ("Y0", t_gapless, 1000.0, (5.0, 20.0, 80.0), (0.0, 0.5, 1.0)),
                                           ("Si", t_si, 1000.0, (5.0, 20.0, 80.0), (0.0, 0.5, 1.0)),
                                           ("Al", t_4000, 4000.0, (20.0,), (1.0,))):
        for spec in ("k_parallel", "k_total"):
            for lam in lams:
                for q in charges:
                    sid = f"p1_{mat}_{int(e0)}_{spec}_Q{q:g}_L{lam:g}"
                    out.append({"id": sid, "material": mat, "E0": e0, "spec": spec, "charge": q, "lam": lam,
                                "losses": losses})
    return out


def required_points():
    """V8 required points: Al, 1000 eV, k_total, Q_h = 1: (80 A, 99.8 eV), (20 A, 1.0 eV), (20 A, 15.6 eV)."""
    return [("p1_Al_1000_k_total_Q1_L80", 99.8), ("p1_Al_1000_k_total_Q1_L20", 1.0),
            ("p1_Al_1000_k_total_Q1_L20", 15.6)]


def i3_losses(table):
    return _family_losses(table, "MS5", "Al")


# ----------------------------------------------------------------------------- task list

def prelim_rule(s):
    return mp3.depth_rule(s["losses"], s["E0"], s["lam"], order=10, phase=math.inf)


def rule_for(s, variant):
    return mp3.variant_rule(s["losses"], s["E0"], s["lam"], variant)


def node_indices(rule, variant):
    if variant == "C":
        return np.flatnonzero(mp3.extension_mask(rule))
    return np.arange(rule.n)


def plan(table):
    tasks = []
    sets = path1_sets(table)
    priority = {sid for sid, _ in required_points()}
    for s in sets:
        tasks.append({"id": f"prelim:{s['id']}", "kind": "prelim", "set": s["id"], "prio": 0})
    envs = {}
    for s in sets:
        envs.setdefault((s["material"], s["E0"], s["spec"]), s["losses"])
    for (mat, e0, spec), losses in envs.items():
        tasks.append({"id": f"env:{mat}_{int(e0)}_{spec}", "kind": "env", "material": mat, "E0": e0, "spec": spec,
                      "prio": 0})
    for mat in I3_MATERIALS:
        key = (mat, 1000.0, "k_parallel")
        if key not in envs:
            tasks.append({"id": f"env:{mat}_1000_k_parallel", "kind": "env", "material": mat, "E0": 1000.0,
                          "spec": "k_parallel", "prio": 0})
    for s in sets:
        for variant in ("base", "A", "B", "C"):
            rule = rule_for(s, variant)
            idx = node_indices(rule, variant)
            nchunks = max(1, math.ceil(idx.size / CHUNK_P1))
            prio = (1 if s["id"] in priority else 3) + (0 if variant == "base" else 1)
            for c in range(nchunks):
                tasks.append({"id": f"node:{s['id']}:{variant}:{c:04d}", "kind": "node", "set": s["id"],
                              "variant": variant, "chunk": c, "nchunks": nchunks, "prio": prio})
    losses = i3_losses(table)
    # Path 3 is computed loss by loss (the frozen route treats each loss independently, so a loss
    # computed alone has the same value as in the batch).  Order: lambda = 20 A (routine) first.
    order = [losses.index(t) for t in I3_LOSS_ORDER]
    for li, lam in enumerate(I3_LAMBDAS[::-1]):
        for rank, j in enumerate(order):
            for variant in mp3.VARIANTS_PATH3:
                rule = mp3.variant_rule(losses, 1000.0, lam, variant)
                idx = node_indices(rule, variant)
                nchunks = max(1, math.ceil(idx.size / CHUNK_I3))
                # lambda = 20 A: the base of the lowest loss first (it sets the V9 scale), then whole points.
                prio = 10 + 100 * li + (0 if (j == losses.index(I3_SCALE_LOSS) and variant == "base") else 1 + rank)
                for c in range(nchunks):
                    tasks.append({"id": f"i3:L{lam:g}:T{j}:{variant}:{c:04d}", "kind": "i3", "lam": lam,
                                  "variant": variant, "loss_index": j, "chunk": c, "nchunks": nchunks, "prio": prio})
    tasks.sort(key=lambda t: (t["prio"], t["id"]))
    return tasks


def task_path(work, task_id):
    kind, rest = task_id.split(":", 1)
    return Path(work) / kind / (rest.replace(":", "__") + (".npz" if kind in ("node", "i3") else ".json"))


# ----------------------------------------------------------------------------- task execution

def _sets_by_id(table):
    return {s["id"]: s for s in path1_sets(table)}


def _atomic_write_bytes(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _atomic_npz(path, **arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + f".tmp{os.getpid()}.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def run_prelim(table, s, work):
    mat = mp3.material_from_table(s["material"], table)
    rule = prelim_rule(s)
    losses = np.asarray(s["losses"])
    vals = np.full((rule.n, 11, losses.size), np.nan)
    used = []
    start = time.perf_counter()
    for i, t in enumerate(rule.t):
        remaining = np.arange(losses.size)
        for epsabs, epsrel, limit in PRELIM_FALLBACK:
            v, e, st = mp3.fixed_depth_node(mat, s["E0"], losses[remaining], s["lam"] * t, s["spec"], epsabs,
                                            epsrel, limit, hole_charge=s["charge"])
            ok = np.array([not x for x in st])
            vals[i][:, remaining[ok]] = v[:, ok]
            for j in remaining[ok]:
                used.append([int(i), int(j), epsabs])
            remaining = remaining[~ok]
            if remaining.size == 0:
                break
        if remaining.size:
            raise RuntimeError(f"{s['id']}: preliminary node {i} failed at every fallback accuracy")
    r_total = np.array([math.fsum(rule.alpha * vals[:, 10, j]) for j in range(losses.size)])
    loose = [u for u in used if u[2] != PRELIM_FALLBACK[0][0]]
    out = {"set": s["id"], "nodes": rule.n, "R_total": r_total.tolist(),
           "max_abs_R_total": float(np.max(np.abs(r_total))), "looser_points": loose,
           "seconds": time.perf_counter() - start}
    _atomic_write_bytes(task_path(work, f"prelim:{s['id']}"), (json.dumps(out, indent=1) + "\n").encode("ascii"))


def run_env(table, task, work):
    mat = mp3.material_from_table(task["material"], table)
    losses = sorted({t for s in path1_sets(table) if (s["material"], s["E0"], s["spec"]) ==
                     (task["material"], task["E0"], task["spec"]) for t in s["losses"]})
    if task["material"] in I3_MATERIALS and task["E0"] == 1000.0 and task["spec"] == "k_parallel":
        losses = sorted(set(losses) | set(i3_losses(table)))
    start = time.perf_counter()
    out = {"material": task["material"], "E0": task["E0"], "spec": task["spec"], "points": {}}
    for t in losses:
        try:
            env = mp3.envelope(mat, task["E0"], t, task["spec"])
        except Exception as exc:  # noqa: BLE001 - an envelope that cannot be formed holds the comparison
            env = None
            out.setdefault("errors", {})[float(t).hex()] = f"{type(exc).__name__}: {exc}"
        out["points"][float(t).hex()] = env
    out["seconds"] = time.perf_counter() - start
    _atomic_write_bytes(task_path(work, task["id"]), (json.dumps(out, indent=1) + "\n").encode("ascii"))


def run_node(table, task, work):
    s = _sets_by_id(table)[task["set"]]
    prelim = json.loads(task_path(work, f"prelim:{s['id']}").read_text("ascii"))
    d_floor = mp3.V8_ABS * prelim["max_abs_R_total"]
    mat = mp3.material_from_table(s["material"], table)
    rule = rule_for(s, task["variant"])
    idx = node_indices(rule, task["variant"])
    chunk = idx[task["chunk"] * CHUNK_P1:(task["chunk"] + 1) * CHUNK_P1]
    eps_all = mp3.node_epsabs(rule, d_floor, 1.0, safety=FINAL_SAFETY)
    losses = np.asarray(s["losses"])
    vals = np.full((chunk.size, 11, losses.size), np.nan)
    errs = np.full((chunk.size, 11, losses.size), np.nan)
    status = np.full((chunk.size, losses.size), "", dtype=object)
    seconds = np.zeros(chunk.size)
    failed_at = {}
    for k, i in enumerate(chunk):
        t0 = time.perf_counter()
        # A loss that failed at an earlier node of this chunk is already a failed point of the
        # path; later nodes of the chunk omit it from the call (deterministic for a fixed chunk).
        active = np.array([j for j in range(losses.size) if j not in failed_at], dtype=int)
        if active.size:
            v, e, st = mp3.fixed_depth_node(mat, s["E0"], losses[active], s["lam"] * rule.t[i], s["spec"],
                                            float(eps_all[i]), hole_charge=s["charge"])
            vals[k][:, active], errs[k][:, active] = v, e
            for jj, j in enumerate(active):
                status[k, j] = st[jj]
                if st[jj]:
                    failed_at[j] = float(rule.t[i])
        for j, t_fail in failed_at.items():
            if j not in active:
                status[k, j] = f"skipped: failed at t={t_fail!r} earlier in this chunk"
        seconds[k] = time.perf_counter() - t0
    _atomic_npz(task_path(work, task["id"]), index=chunk, t=rule.t[chunk], alpha=rule.alpha[chunk],
                epsabs=eps_all[chunk], values=vals, errors=errs, status=status.astype(str), seconds=seconds,
                d_floor=np.array(d_floor), n_alloc=np.array(rule.n))


def run_i3(table, task, work):
    ref = mp3.import_route_b()
    mp3.configure_route_b(ref, task["variant"])
    losses = np.asarray(i3_losses(table))
    mats = [mp3.material_from_table(m, table) for m in I3_MATERIALS]
    rule = mp3.variant_rule(losses, 1000.0, task["lam"], task["variant"])  # nodes of the full set P
    idx = node_indices(rule, task["variant"])
    chunk = idx[task["chunk"] * CHUNK_I3:(task["chunk"] + 1) * CHUNK_I3]
    j = task["loss_index"]
    vals = np.full((chunk.size, len(mats), 10, 1), np.nan)
    seconds = np.zeros(chunk.size)
    for k, i in enumerate(chunk):
        t0 = time.perf_counter()
        vals[k] = mp3.i3_node(ref, mats, 1000.0, losses[j:j + 1], task["lam"] * rule.t[i])
        seconds[k] = time.perf_counter() - t0
    _atomic_npz(task_path(work, task["id"]), index=chunk, t=rule.t[chunk], alpha=rule.alpha[chunk], values=vals,
                seconds=seconds)


def run_task(table, task, work):
    kind = task["kind"]
    if kind == "prelim":
        run_prelim(table, _sets_by_id(table)[task["set"]], work)
    elif kind == "env":
        run_env(table, task, work)
    elif kind == "node":
        run_node(table, task, work)
    elif kind == "i3":
        run_i3(table, task, work)
    else:
        raise ValueError(kind)


def _ready(task, work):
    if task["kind"] == "node":
        return task_path(work, f"prelim:{task['set']}").exists()
    return True


def work_loop(table, work, kinds, sets, variants, max_tasks):
    tasks = json.loads((Path(work) / "plan.json").read_text("ascii"))
    claims = Path(work) / "claims"
    claims.mkdir(parents=True, exist_ok=True)
    done = 0
    for task in tasks:
        if task["kind"] not in kinds:
            continue
        if sets and not any(x in task.get("set", task["id"]) for x in sets):
            continue
        if variants and task.get("variant") not in variants:
            continue
        if task_path(work, task["id"]).exists() or not _ready(task, work):
            continue
        claim = claims / (task["id"].replace(":", "__") + ".claim")
        try:
            fd = os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            continue
        os.write(fd, f"{os.getpid()} {time.time()}\n".encode("ascii"))
        os.close(fd)
        start = time.perf_counter()
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "task", "--work", str(work), task["id"]],
                              stdin=subprocess.DEVNULL, capture_output=True, text=True)
        took = time.perf_counter() - start
        state = "ok" if proc.returncode == 0 and task_path(work, task["id"]).exists() else f"FAILED rc={proc.returncode}"
        print(f"{time.strftime('%H:%M:%S')} {task['id']} {state} {took:.1f}s", flush=True)
        if proc.returncode != 0:
            print(proc.stderr[-3000:], flush=True)
            (Path(work) / "failed").mkdir(exist_ok=True)
            (Path(work) / "failed" / (task["id"].replace(":", "__") + ".txt")).write_text(proc.stderr, "utf-8")
        else:
            claim.unlink(missing_ok=True)
        done += 1
        if max_tasks and done >= max_tasks:
            break


def status(work):
    tasks = json.loads((Path(work) / "plan.json").read_text("ascii"))
    by = {}
    for task in tasks:
        key = (task["kind"], task.get("set", f"L{task.get('lam', '')}"), task.get("variant", ""))
        tot, ok = by.get(key, (0, 0))
        by[key] = (tot + 1, ok + int(task_path(work, task["id"]).exists()))
    for key in sorted(by):
        tot, ok = by[key]
        if ok < tot:
            print(f"{key}: {ok}/{tot}")
    total = sum(v[0] for v in by.values())
    finished = sum(v[1] for v in by.values())
    print(f"tasks done {finished}/{total}")


# ----------------------------------------------------------------------------- assembly

def _load_nodes(work, prefix):
    files = sorted(Path(work).glob(prefix))
    return [np.load(f, allow_pickle=False) for f in files]


def _gather(chunks, n):
    """Concatenate node chunks in node order; returns None if a node is missing."""
    if not chunks:
        return None
    idx = np.concatenate([c["index"] for c in chunks])
    if idx.size != n or np.unique(idx).size != n:
        return None
    order = np.argsort(idx)
    out = {}
    for key in chunks[0].files:
        if key in ("d_floor", "n_alloc"):
            out[key] = chunks[0][key]
            continue
        out[key] = np.concatenate([c[key] for c in chunks])[order]
    return out


def _sum_columns(alpha, values):
    """Path values (11, n_losses) and rounding bound from node values (N, 11, n_losses)."""
    return mp3.weighted_sum(alpha, values)


def assemble_path1(table, work):
    sets = path1_sets(table)
    out_sets, failures = [], []
    required = {(sid, float(t).hex()) for sid, t in required_points()}
    for s in sets:
        losses = np.asarray(s["losses"])
        n_l = losses.size
        rules = {v: rule_for(s, v) for v in ("base", "A", "B", "C")}
        prelim_file = task_path(work, f"prelim:{s['id']}")
        prelim = json.loads(prelim_file.read_text("ascii")) if prelim_file.exists() else None
        data = {}
        for v, rule in rules.items():
            n = node_indices(rule, v).size
            data[v] = _gather(_load_nodes(work, f"node/{s['id']}__{v}__*.npz"), n)
        env_file = task_path(work, f"env:{s['material']}_{int(s['E0'])}_{s['spec']}")
        env = json.loads(env_file.read_text("ascii")) if env_file.exists() else None
        entry = {"id": s["id"], "material": s["material"], "electron_energy_eV": s["E0"], "epsilon_momentum": s["spec"],
                 "hole_charge": s["charge"], "hole_charge_hex": float(s["charge"]).hex(), "weight_length_A": s["lam"],
                 "weight_length_hex": float(s["lam"]).hex(), "losses_hex": [float(t).hex() for t in losses],
                 "nodes": {v: int(node_indices(r, v).size) for v, r in rules.items()},
                 "preliminary": None if prelim is None else {"nodes": prelim["nodes"],
                                                             "max_abs_R_total": prelim["max_abs_R_total"],
                                                             "looser_points": prelim["looser_points"]},
                 "complete": all(d is not None for d in data.values()) and prelim is not None}
        base = data["base"]
        if base is None or prelim is None:
            entry["points"] = [{"loss_hex": float(t).hex(), "loss_eV": float(t), "status": "not_computed"}
                               for t in losses]
            out_sets.append(entry)
            failures += [{"set": s["id"], "loss_hex": float(t).hex(), "status": "not_computed"} for t in losses]
            continue
        st = base["status"]
        node_fail = np.array([any(st[:, j] != "") for j in range(n_l)])
        R, rho = _sum_columns(base["alpha"], base["values"])
        E, _ = _sum_columns(base["alpha"], base["errors"])
        E = E * (1 + 4 * mp3.EPS * base["alpha"].size)
        sums = {}
        variant_fail = {}
        for v in ("A", "B", "C"):
            d = data[v]
            if d is None:
                sums[v] = None
                variant_fail[v] = np.ones(n_l, bool)
                continue
            variant_fail[v] = np.array([any(d["status"][:, j] != "") for j in range(n_l)])
            if v == "C":
                ext, _ = _sum_columns(d["alpha"], d["values"])
                sums[v] = R + ext
            else:
                sums[v], _ = _sum_columns(d["alpha"], d["values"])
        ok_total = ~node_fail
        scale_vals = np.abs(R[10])
        scale_parts = {"final": [bool(x) for x in ok_total]}
        if ok_total.all():
            scale = float(np.max(scale_vals))
            scale_source = "final"
        else:
            # A failed point has no final value: its preliminary value enters the scale.
            pre = np.abs(np.asarray(prelim["R_total"]))
            scale = float(max(np.max(np.where(ok_total, scale_vals, 0.0)), np.max(np.where(ok_total, 0.0, pre))))
            scale_source = "final where computed, preliminary at failed points"
        D = mp3.V8_REL * np.abs(R) + mp3.V8_ABS * scale
        # Allocation check alpha_i e_ij <= f_i eta D_j on every base node, column and loss.
        share = mp3.ETA / base["alpha"].size * D
        with np.errstate(invalid="ignore", divide="ignore"):
            alloc = base["alpha"][:, None, None] * base["errors"] / share[None]
        alloc_max = np.nanmax(np.where(np.isfinite(alloc), alloc, np.nan), axis=0)
        variant_alloc = {}
        for v in ("A", "B", "C"):
            d = data[v]
            if d is None:
                continue
            n_alloc = rules[v].n
            sh = mp3.ETA / n_alloc * D
            with np.errstate(invalid="ignore", divide="ignore"):
                a = d["alpha"][:, None, None] * d["errors"] / sh[None]
            variant_alloc[v] = np.nanmax(np.where(np.isfinite(a), a, np.nan), axis=0)
        points = []
        for j, t in enumerate(losses):
            key_hex = float(t).hex()
            p = {"loss_hex": key_hex, "loss_eV": float(t), "required": (s["id"], key_hex) in required}
            reasons = []
            if node_fail[j]:
                bad = [(float(base["t"][i]), str(st[i, j])) for i in range(st.shape[0])
                       if st[i, j] and not str(st[i, j]).startswith("skipped")]
                reasons.append(f"base node failure at {len(bad)} nodes, first t={bad[0][0]!r}: {bad[0][1][:300]}")
                p["failed_nodes"] = {"base": len(bad), "first_t": bad[0][0], "first_message": bad[0][1][:500],
                                     "skipped_after_failure": int(sum(str(x).startswith("skipped") for x in st[:, j]))}
            for v in ("A", "B", "C"):
                if sums[v] is None:
                    reasons.append(f"variant {v} not computed")
                elif variant_fail[v][j]:
                    d = data[v]
                    bad = [i for i in range(d["status"].shape[0])
                           if d["status"][i, j] and not str(d["status"][i, j]).startswith("skipped")]
                    reasons.append(f"variant {v} node failure at {len(bad)} nodes: {str(d['status'][bad[0], j])[:300]}")
            pt_env = None if env is None else env["points"].get(key_hex)
            if pt_env is None:
                reasons.append("tail envelope not available (comparison held)")
            if reasons:
                p["status"] = "failed" if node_fail[j] or any("node failure" in r for r in reasons) else "held"
                p["reasons"] = reasons
                failures.append({"set": s["id"], "loss_hex": key_hex, "status": p["status"], "reasons": reasons})
                points.append(p)
                continue
            cols = {}
            for c, name in enumerate(mp3.OUTPUTS):
                dA = float(sums["A"][c, j] - R[c, j])
                dB = float(sums["B"][c, j] - R[c, j])
                dC = float(sums["C"][c, j] - R[c, j])
                tail = mp3.tail_bound(pt_env, name, s["lam"], charge=s["charge"])
                u = max(abs(dA), abs(dB), abs(dC)) + float(E[c, j]) + float(rho[c, j]) + tail
                u *= 1 + 8 * mp3.EPS
                cols[name] = {"value": float(R[c, j]), "D": float(D[c, j]), "u": u,
                              "delta_A": dA, "delta_B": dB, "delta_C": dC, "node_error": float(E[c, j]),
                              "rounding": float(rho[c, j]), "tail": tail,
                              "allocation_max_ratio": float(alloc_max[c, j])}
            worst = max(cols, key=lambda n: cols[n]["u"] / cols[n]["D"])
            alloc_ok = all(x["allocation_max_ratio"] <= 1.0 for x in cols.values())
            u_ok = all(x["u"] <= x["D"] / 3 for x in cols.values())
            p["status"] = "ok" if (alloc_ok and u_ok) else "path_criterion_failed"
            p["max_u_over_D"] = cols[worst]["u"] / cols[worst]["D"]
            p["max_u_over_D_column"] = worst
            p["allocation_ok"] = alloc_ok
            # Variant node sets checked with their own N (reported, not part of the criterion).
            p["allocation_max_ratio_variants"] = {v: float(np.max(a[:, j])) for v, a in variant_alloc.items()}
            p["u_le_D_over_3"] = u_ok
            p["columns"] = cols
            if p["status"] != "ok":
                failures.append({"set": s["id"], "loss_hex": key_hex, "status": p["status"],
                                 "allocation_ok": alloc_ok, "u_le_D_over_3": u_ok, "max_u_over_D": p["max_u_over_D"]})
            points.append(p)
        entry["scale_max_abs_R_total"] = scale
        entry["scale_source"] = scale_source
        entry["node_seconds"] = {v: float(np.sum(d["seconds"])) for v, d in data.items() if d is not None}
        entry["points"] = points
        out_sets.append(entry)
    return out_sets, failures


def assemble_path3(table, work):
    losses = np.asarray(i3_losses(table))
    out, failures = [], []
    for lam in I3_LAMBDAS:
        rules = {v: mp3.variant_rule(losses, 1000.0, lam, v) for v in mp3.VARIANTS_PATH3}
        # sums[v][j]: (n_mat, 10) path sums for loss j, or None while that loss is incomplete
        sums = {v: [None] * losses.size for v in mp3.VARIANTS_PATH3}
        seconds = {v: [0.0] * losses.size for v in mp3.VARIANTS_PATH3}
        for v, rule in rules.items():
            n = node_indices(rule, v).size
            for j in range(losses.size):
                d = _gather(_load_nodes(work, f"i3/L{lam:g}__T{j}__{v}__*.npz"), n)
                if d is None:
                    continue
                flat = d["values"].reshape(d["values"].shape[0], -1)
                tot = np.array([math.fsum(col) for col in (d["alpha"][:, None] * flat).T])
                sums[v][j] = tot.reshape(d["values"].shape[1:3])
                seconds[v][j] = float(np.sum(d["seconds"]))
        for m, mat in enumerate(I3_MATERIALS):
            env_file = task_path(work, f"env:{mat}_1000_k_parallel")
            env = json.loads(env_file.read_text("ascii")) if env_file.exists() else None
            for q in (0.0, 0.5, 1.0):
                def cols(arr):
                    a = arr[m]
                    return np.concatenate((a, [a[8] + q * a[9]]))
                base = [None if x is None else cols(x) for x in sums["base"]]
                have = [b for b in base if b is not None]
                complete_scale = all(b is not None for b in base)
                scale = float(max(abs(b[10]) for b in have)) if have else None
                entry = {"material": mat, "weight_length_A": lam, "weight_length_hex": float(lam).hex(),
                         "hole_charge": q, "hole_charge_hex": float(q).hex(), "epsilon_momentum": "k_parallel",
                         "electron_energy_eV": 1000.0,
                         "nodes": {v: int(node_indices(r, v).size) for v, r in rules.items()},
                         "scale_max_abs_K_total": scale if complete_scale else None,
                         "scale_lower_bound": scale, "scale_complete": complete_scale, "points": []}
                for j, t in enumerate(losses):
                    key_hex = float(t).hex()
                    pt = {"loss_hex": key_hex, "loss_eV": float(t),
                          "node_seconds": {v: seconds[v][j] for v in mp3.VARIANTS_PATH3}}
                    missing = [v for v in mp3.VARIANTS_PATH3 if sums[v][j] is None]
                    if base[j] is not None:
                        R = base[j]
                        pt["base_values"] = {name: float(R[c]) for c, name in enumerate(mp3.OUTPUTS)}
                        pt["changes"] = {}
                        for v in ("A", "B", "C", "D", "E"):
                            if sums[v][j] is None:
                                continue
                            arr = cols(sums[v][j])
                            if v == "C":
                                arr = R + arr  # C adds the nodes on [40, 50] to the base nodes
                            pt["changes"][v] = {name: float(arr[c] - R[c]) for c, name in enumerate(mp3.OUTPUTS)}
                    if base[j] is None:
                        pt["status"] = "not_computed"
                        pt["missing_variants"] = missing
                        failures.append({"path": 3, "material": mat, "weight_length_A": lam, "hole_charge": q,
                                         "loss_hex": key_hex, "status": pt["status"], "missing_variants": missing})
                        entry["points"].append(pt)
                        continue
                    # Without all nine base values the scale is a lower bound, so the width is smaller
                    # than the registered one: the path criterion is then only checked conservatively.
                    width = mp3.V9_REL * np.abs(R) + mp3.V9_ABS * scale
                    pt_env = None if env is None else env["points"].get(key_hex)
                    pt["columns"] = {}
                    worst = 0.0
                    for c, name in enumerate(mp3.OUTPUTS):
                        col = {"value": float(R[c]), "width": float(width[c]),
                               "changes": {v: pt["changes"][v][name] for v in pt["changes"]}}
                        col["tail"] = None if pt_env is None else mp3.tail_bound(pt_env, name, lam, charge=q)
                        if col["changes"]:
                            col["max_change_over_width"] = max(abs(x) for x in col["changes"].values()) / float(width[c])
                            worst = max(worst, col["max_change_over_width"])
                        pt["columns"][name] = col
                    tail_ok = pt_env is not None and all(
                        pt["columns"][n]["tail"] <= 0.1 * pt["columns"][n]["width"] for n in mp3.OUTPUTS)
                    pt["max_change_over_width"] = worst
                    if missing:
                        pt["status"] = "incomplete"
                        pt["missing_variants"] = missing
                    elif not (worst <= 0.1 and tail_ok):
                        pt["status"] = "path_criterion_failed"
                    elif complete_scale:
                        pt["status"] = "ok"
                    else:
                        pt["status"] = "criterion_met_with_partial_scale"
                        pt["reason"] = "the V9 scale needs the base values of all nine losses (lower bound used)"
                    if pt["status"] != "ok":
                        failures.append({"path": 3, "material": mat, "weight_length_A": lam, "hole_charge": q,
                                         "loss_hex": key_hex, "status": pt["status"], "max_change_over_width": worst,
                                         "missing_variants": missing})
                    entry["points"].append(pt)
                out.append(entry)
    return out, failures


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assemble(table, work, out_path):
    p1, f1 = assemble_path1(table, work)
    p3, f3 = assemble_path3(table, work)
    import scipy
    import xpsloss
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True,
                                stdin=subprocess.DEVNULL).stdout.strip()
    except OSError:
        commit = None
    req = []
    for sid, t in required_points():
        entry = next(e for e in p1 if e["id"] == sid)
        pt = next(p for p in entry["points"] if p["loss_hex"] == float(t).hex())
        req.append({"set": sid, "loss_hex": float(t).hex(), "loss_eV": t, "status": pt["status"],
                    "max_u_over_D": pt.get("max_u_over_D")})
    doc = {
        "schema": SCHEMA, "prereg": "docs/M3A_PREREG.md", "expected_states_sha256": TABLE_SHA256,
        "generator": "benchmarks/numerics/m3a/generate_paths_m3a.py", "generator_sha256": _sha(__file__),
        "helper": "tests/m3a_paths.py", "helper_sha256": _sha(mp3.__file__),
        "route_b_reference": "tests/xps_route_b_reference.py",
        "route_b_reference_sha256": _sha(ROOT / "tests" / "xps_route_b_reference.py"),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                        "xpsloss": xpsloss.__version__, "platform": platform.platform(), "git_head": commit},
        "rules": {"registered_edges": list(mp3.REGISTERED_EDGES), "r": mp3.BASE_R, "decay": mp3.DECAY,
                  "eta": mp3.ETA, "f_i": "1/N", "final_request_safety": FINAL_SAFETY,
                  "node_epsrel": mp3.NODE_EPSREL, "node_limit": mp3.NODE_LIMIT,
                  "preliminary_fixed_depth_settings": [list(x) for x in PRELIM_FALLBACK],
                  "variants_path1": {k: {kk: (vv if vv is None or math.isfinite(vv) else str(vv)) for kk, vv in v.items()}
                                     for k, v in mp3.VARIANTS_PATH1.items()},
                  "variants_path3": {"A": "16 nodes", "B": "pi/8, ratio sqrt 2", "C": "r = 50",
                                     "D": "kpar_nodes n_base 24", "E": "k_perp GL 40, Laguerre 160, 1/k tail 64"},
                  "V8": {"D": "1e-7 |R_j| + 1e-8 max_P |R_total|", "path": "u_j <= D_j/3",
                         "comparison": "|Delta_j| + u_j <= D_j"},
                  "V9": {"width": "1e-6 |K_j,I3| + 1e-6 max_9 |K_total,I3|", "path": "each change <= width/10"}},
        "required_points": req,
        "path1": p1,
        "path3": p3,
        "not_ok": {"path1": f1, "path3": f3},
    }
    text = json.dumps(doc, indent=1, sort_keys=True, ensure_ascii=True, allow_nan=False) + "\n"
    data = text.encode("ascii")
    Path(out_path).write_bytes(data)
    print(json.dumps({"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data), "path1_not_ok": len(f1),
                      "path3_not_ok": len(f3), "required": req}, indent=1))


# ----------------------------------------------------------------------------- command line

def main(argv=None):
    parser = argparse.ArgumentParser(description="M3a paths 1 and 3 (R2)")
    parser.add_argument("command", choices=("plan", "work", "task", "status", "assemble"))
    parser.add_argument("task_id", nargs="?")
    parser.add_argument("--work", required=True)
    parser.add_argument("--kinds", default="prelim,env,node,i3")
    parser.add_argument("--sets", default="")
    parser.add_argument("--variants", default="")
    parser.add_argument("--max-tasks", type=int, default=0)
    parser.add_argument("--out", default=str(OUT))
    args = parser.parse_args(argv)
    work = Path(args.work)
    if args.command == "status":
        status(work)
        return
    table = load_table()
    if args.command == "plan":
        work.mkdir(parents=True, exist_ok=True)
        tasks = plan(table)
        (work / "plan.json").write_text(json.dumps(tasks, indent=0) + "\n", "ascii")
        kinds = {}
        for t in tasks:
            kinds[t["kind"]] = kinds.get(t["kind"], 0) + 1
        print(json.dumps({"tasks": len(tasks), "by_kind": kinds}))
    elif args.command == "task":
        tasks = {t["id"]: t for t in json.loads((work / "plan.json").read_text("ascii"))}
        run_task(table, tasks[args.task_id], work)
    elif args.command == "work":
        work_loop(table, work, set(args.kinds.split(",")), [x for x in args.sets.split(",") if x],
                  [x for x in args.variants.split(",") if x], args.max_tasks)
    elif args.command == "assemble":
        assemble(table, work, args.out)


if __name__ == "__main__":
    main()
