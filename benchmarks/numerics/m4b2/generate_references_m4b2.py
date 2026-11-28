"""R-C: the reference of M4b-2 (prereg docs/M4B2_PREREG.md section 6), written without the package.

Standard library only (``fractions`` for the exact arithmetic, ``float`` for the single
rounding).  ``xpsloss`` is not imported and no conversion function of the product is called.
The input is the frozen table ``expected_states_m4b2.json``; its SHA-256 is checked before
it is read.  The output ``references_m4b2.json`` is ASCII with LF line ends.

    python benchmarks/numerics/m4b2/generate_references_m4b2.py            write the file
    python benchmarks/numerics/m4b2/generate_references_m4b2.py --check    regenerate, compare bytes

What is fixed here, for the rows of prereg 6:

* the conversion: C as an irreducible fraction, RN(C), RN(C - B) once per distinct input
  coordinate (then handed to every output that uses the coordinate), ``rounding_max_eV``,
  the implied offset, the kinetic energy of the representative position, E0 minus it;
* the inputs of the physical control (ascending KE edges, the KE primary, the order);
* the map of the BE display bins to the KE bins and the zero loss of every region as an
  exact rational of the rounded-KE problem (a binned primary is spread over its KE widths);
* the both-sided support of prereg 2.3: S_B, S_K, S, the width of every bin rounded upward,
  complete or reached, the mismatch counts and the first examples, the warning, the state
  without a tolerance;
* the first colliding pair and its roles for the rows that must be rejected (prereg 2.2);
* the two edge hashes of the CLI rows.

The section ``keys`` of a row holds one entry per ledger key of the frozen table whose oracle
is R-C; every such entry is exact (a rational or a binary64), so the reference uncertainty is
zero and ``no_reference`` is empty.  The self checks and the known wrong answers of prereg 6
are evaluated on every run; a failed one stops the generation.
"""
from fractions import Fraction as Fr
import hashlib
import json
import math
from pathlib import Path
import struct
import sys

HERE = Path(__file__).resolve().parent
TABLE_PATH = HERE/"expected_states_m4b2.json"
OUT_PATH = HERE/"references_m4b2.json"
TABLE_SHA256 = "9c97907dab69a00e110414d6862d0d061bc80212000b8f3180f905bcb28108ea"
SCHEMA = "xpsloss-m4b2-references-rc-v1"
EPS = Fr(1, 2**52)

GEOMETRY_ROWS = ("BE1", "BE1r", "BE2", "BE3", "BE3r", "BE4", "BE5a", "BE5", "BE6", "BE7", "BE8", "BE9", "BR1", "BZ0",
                 "BT1a", "BT1b", "BT1c", "BDS", "BT2c", "BT2a",
                 "BS1", "BS4", "BS1t", "BS1d", "BS2", "BS3", "BS3n2", "BS5", "BS6",
                 "BW1", "BW2", "BW3", "BW4",
                 "BC01s", "BC01", "BC02", "BC07", "BC12")
COLLISION_ROWS = ("BN01", "BN23", "BN18", "BC06")
SUPPORT_ROWS = ("BT2b", "BS1x", "BS1dl")


class ReferenceError_(Exception):
    """A self check of the reference failed; nothing is written."""


def need(condition, message):
    if not condition:
        raise ReferenceError_(message)


# ------------------------------------------------------------------ exact numbers
def q(x):
    """The binary64 x read as an exact rational."""
    return Fr(float(x))


def rn(value):
    """One rounding of an exact rational to the nearest binary64 (ties to even)."""
    return float(Fr(value))


def up64(value):
    """The smallest binary64 that is not below the exact rational."""
    x = rn(value)
    return math.nextafter(x, math.inf) if Fr(x) < value else x


def fx(x):
    """A binary64 in the file: its hexadecimal form."""
    x = float(x)
    return {"hex": "0x0.0p+0" if x == 0.0 and math.copysign(1.0, x) > 0 else x.hex()}


def fq(value):
    """An exact rational in the file: 'numerator/positive denominator' (integers as p/1)."""
    value = Fr(value)
    return f"{value.numerator}/{value.denominator}"


def read(obj):
    """A number of the table: the hexadecimal form decides, the decimal form is checked."""
    if isinstance(obj, dict):
        x = float.fromhex(obj["hex"])
        need(x == float(obj["value"]), f"hex and decimal differ: {obj}")
        return x
    return float(obj)


def read_list(obj):
    xs = [float.fromhex(h) for h in obj["hex"]]
    need(xs == [float(v) for v in obj["values"]], "hex and decimal lists differ")
    return xs


def relation(x, y):
    return "<" if x < y else ("=" if x == y else ">")


# ------------------------------------------------------------------ the table
def load_table():
    raw = TABLE_PATH.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    need(digest == TABLE_SHA256, f"the table is not the frozen one: {digest}")
    return json.loads(raw.decode("utf-8"))


def axis_constant(axis):
    """C of prereg 2.1 from the binary64 inputs, never rounded before use."""
    ref = axis["reference"]
    hv = q(read(axis["photon_energy_eV"]))
    if ref["kind"] == "photon_and_offset":
        return hv - q(read(ref["kinetic_energy_offset_eV"])), hv, None
    need(ref["kind"] == "paired_anchor", ref["kind"])
    c = q(read(ref["binding_anchor_eV"])) + q(read(ref["kinetic_anchor_eV"]))
    return c, hv, hv - c


def option(argv, name, count):
    k = argv.index(name)
    return argv[k + 1:k + 1 + count]


def expand_grid(compact):
    """min, max, step -> the edges (exact interval count, then one rounding per point)."""
    lo, hi, step = (Fr(str(compact[k])) for k in ("min", "max", "step"))
    need(lo < hi and step > 0, "compact edges")
    count = (hi - lo)/step
    need(count.denominator == 1, "the interval count is not an integer")
    n = int(count)
    return [rn(lo + (hi - lo)*i/n) for i in range(n + 1)]


def row_inputs(table, row):
    """The numbers of a row that enter the geometry: axis, primary, edges, window, order, E0."""
    i = row["inputs"]
    if row["entry"] in ("cli",):
        s = i["spectrum"]
        measurement = s["measurement"]
        primary_name = s["primary"]
        edges_spec = s["output_binding_edges_eV"]
        window = tuple(float(x) for x in option(i["argv"], "--loss-window-eV", 2))
        max_order = int(option(i["argv"], "--max-order", 1)[0])
        e0 = float(i["config_electron_energy_eV"])
    else:
        measurement = i["energy_axis"]
        primary_name = i["primary"]
        edges_spec = i["output_edges_eV"]
        if "kernel" in i:
            window = tuple(read(x) for x in table["kernels"][i["kernel"]]["window_eV"])
            e0 = None
        else:
            window = tuple(read(x) for x in i["loss_window_eV"])
            e0 = float(i["electron_energy_eV"])
        max_order = i["max_order"]
    if isinstance(edges_spec, dict):
        registered = table["binding_output_edges"][edges_spec["compact_form_of"]]
        edges = expand_grid(registered["compact"])
        need([x.hex() for x in edges] == [float.fromhex(h).hex() for h in registered["hex"]],
             "the expanded compact edges are not the registered binary64 values")
    else:
        edges = read_list(table["binding_output_edges"][edges_spec])
    p = table["binding_primaries"][primary_name]
    kind = p["kind"]
    positions = read_list(p["binding_energies_eV" if kind == "delta" else "binding_edges_eV"])
    intensities = read_list(p["intensities"])
    return dict(axis=table["axes"][measurement["axis"]], axis_name=measurement["axis"],
                reference_binding=read(measurement["primary_reference_binding_eV"]),
                kind=kind, positions=positions, intensities=intensities, edges=edges,
                window=window, max_order=max_order, e0=e0, entry=row["entry"])


# ------------------------------------------------------------------ conversion and injectivity
def coordinates(inp):
    """Distinct BE coordinates with their roles: primary in input order, then the output edges."""
    name = "primary.binding_energies_eV" if inp["kind"] == "delta" else "primary.binding_edges_eV"
    seen = {}
    for k, b in enumerate(inp["positions"]):
        seen.setdefault(b, []).append(f"{name}[{k}]")
    for k, b in enumerate(inp["edges"]):
        seen.setdefault(b, []).append(f"output_edges_eV[{k}]")
    return seen


def convert(inp, c):
    """RN(C - B) once per distinct coordinate; the exact rounding error of each."""
    kinetic, error = {}, {}
    for b in coordinates(inp):
        exact = c - q(b)
        kinetic[b] = rn(exact)
        error[b] = abs(Fr(kinetic[b]) - exact)
    return kinetic, error


def first_collision(inp, kinetic):
    """The first adjacent pair, in ascending BE, whose rounded KE is not strictly decreasing."""
    roles = coordinates(inp)
    ordered = sorted(roles)
    for x, y in zip(ordered, ordered[1:]):
        if not kinetic[x] > kinetic[y]:
            return dict(binding=(x, y), kinetic=(kinetic[x], kinetic[y]), roles=(roles[x][0], roles[y][0]))
    return None


# ------------------------------------------------------------------ zero loss
def ke_regions(n_bins):
    return ["outside.low"] + [f"bin[{k}]" for k in range(n_bins)] + ["outside.high"]


def delta_region(e, edges):
    """KE rule of prereg 2.4: bins [e_k, e_k+1), the highest KE bin closed on the right too."""
    if e < edges[0]:
        return "outside.low"
    if e > edges[-1]:
        return "outside.high"
    if e == edges[-1]:
        return f"bin[{len(edges) - 2}]"
    k = max(j for j in range(len(edges) - 1) if edges[j] <= e)
    return f"bin[{k}]"


def overlap_shares(lo, hi, cuts):
    """Shares of the interval [lo, hi] (exact rationals, lo < hi) in the regions cut by ``cuts``."""
    bounds = [None] + list(cuts) + [None]
    out = []
    for left, right in zip(bounds, bounds[1:]):
        a = lo if left is None else max(lo, left)
        b = hi if right is None else min(hi, right)
        out.append(max(Fr(0), b - a)/(hi - lo))
    return out


def zero_loss_kinetic(inp, ke_edges, ke_positions, ke_intensities):
    """Exact zero loss of every KE region of the rounded-KE problem."""
    names = ke_regions(len(ke_edges) - 1)
    total = dict.fromkeys(names, Fr(0))
    if inp["kind"] == "delta":
        for e, s in zip(ke_positions, ke_intensities):
            total[delta_region(e, ke_edges)] += q(s)
    else:
        cuts = [q(e) for e in ke_edges]
        for k, s in enumerate(ke_intensities):
            shares = overlap_shares(q(ke_positions[k]), q(ke_positions[k + 1]), cuts)
            for name, share in zip(names, shares):
                total[name] += q(s)*share
    return total


def zero_loss_binding_overlap(inp):
    """Diagnostic only: a binned primary spread over its exact BE widths (BE region names)."""
    edges = sorted(inp["edges"])
    n = len(edges) - 1
    increasing = inp["edges"][0] < inp["edges"][-1]
    names = (["outside.low_binding"] + [f"bin[{k if increasing else n - 1 - k}]" for k in range(n)]
             + ["outside.high_binding"])
    total = dict.fromkeys(names, Fr(0))
    cuts = [q(e) for e in edges]
    for k, s in enumerate(inp["intensities"]):
        lo, hi = sorted((q(inp["positions"][k]), q(inp["positions"][k + 1])))
        for name, share in zip(names, overlap_shares(lo, hi, cuts)):
            total[name] += q(s)*share
    return total


# ------------------------------------------------------------------ one row
def support_points(inp):
    """P of prereg 2.3: positive deltas, or both edges of the bins of positive intensity (input index kept)."""
    index = {}
    if inp["kind"] == "delta":
        for k, (b, s) in enumerate(zip(inp["positions"], inp["intensities"])):
            if s > 0:
                index.setdefault(b, k)
    else:
        for k, s in enumerate(inp["intensities"]):
            if s > 0:
                index.setdefault(inp["positions"][k], k)
                index.setdefault(inp["positions"][k + 1], k + 1)
    return index


def analyse(inp, constant=None, kinetic_of=None):
    """Everything R-C fixes for a row.  ``constant`` and ``kinetic_of`` exist for the known wrong answers."""
    c, hv, implied = axis_constant(inp["axis"])
    if constant is not None:
        c = constant
    kinetic, error = convert(inp, c)
    if kinetic_of is not None:
        kinetic = {b: kinetic_of(c, b) for b in kinetic}
        error = {b: abs(Fr(kinetic[b]) - (c - q(b))) for b in kinetic}
    out = dict(c=c, hv=hv, implied=implied, kinetic=kinetic, error=error,
               negative=any(e < 0 for e in kinetic.values()), collision=first_collision(inp, kinetic))
    if out["negative"] or out["collision"]:
        return out
    edges_in = inp["edges"]
    n = len(edges_in) - 1
    increasing = edges_in[0] < edges_in[-1]
    need(all((x < y) == increasing for x, y in zip(edges_in, edges_in[1:])), "the output edges are not monotonic")
    ke_edges = sorted(kinetic[b] for b in edges_in)
    bin_map = [n - 1 - k if increasing else k for k in range(n)]
    for k in range(n):       # the KE bin of a display bin holds the images of its two BE edges
        pair = sorted((kinetic[edges_in[k]], kinetic[edges_in[k + 1]]))
        need(pair == ke_edges[bin_map[k]:bin_map[k] + 2], "bin map")
    if inp["kind"] == "delta":
        order = sorted(range(len(inp["positions"])), key=lambda k: kinetic[inp["positions"][k]])   # stable
        ke_positions = [kinetic[inp["positions"][k]] for k in order]
        ke_intensities = [inp["intensities"][k] for k in order]
        primary_direction = None
    else:
        p_increasing = inp["positions"][0] < inp["positions"][-1]
        primary_direction = "increasing" if p_increasing else "decreasing"
        ke_positions = sorted(kinetic[b] for b in inp["positions"])
        ke_intensities = inp["intensities"][::-1] if p_increasing else list(inp["intensities"])
    zero_ke = zero_loss_kinetic(inp, ke_edges, ke_positions, ke_intensities)
    zero = {"outside.low_binding": zero_ke["outside.high"], "outside.high_binding": zero_ke["outside.low"]}
    for k in range(n):
        zero[f"bin[{k}]"] = zero_ke[f"bin[{bin_map[k]}]"]

    # ---- both-sided support (prereg 2.3)
    a, b = q(inp["window"][0]), q(inp["window"][1])
    points = support_points(inp)
    high_edge = [max(edges_in[k], edges_in[k + 1]) for k in range(n)]
    if points:
        b_min = min(points)
        e_max = max(kinetic[p] for p in points)
        need(kinetic[b_min] == e_max, "the lowest BE point of the support has the highest KE")
        u_b = [q(d) - q(b_min) for d in high_edge]
        u_k = [q(e_max) - q(kinetic[d]) for d in high_edge]
        s_b = max(0, math.floor(max(u_b)/a))
        s_k = max(0, math.floor(max(u_k)/a))
    else:
        u_b = u_k = [None]*n
        s_b = s_k = 0
    support = max(s_b, s_k)
    requested = inp["max_order"]
    computed = support if requested == "display_support" else int(requested)
    limit = (computed + 1)*a
    if points:
        width = [up64(max(x, y)) for x, y in zip(u_b, u_k)]
        width_kinetic = [up64(y) for y in u_k]
        complete = [x < limit and y < limit for x, y in zip(u_b, u_k)]
        kinetic_reaches = [not y < limit for y in u_k]
    else:
        width = width_kinetic = [None]*n
        complete = [True]*n
        kinetic_reaches = [False]*n
    examples, computed_count = [], 0
    for order in range(1, computed + 1):
        for bound, c_ in (("lower", order*a), ("upper", order*b)):
            for p in sorted(points):
                for d in sorted(edges_in):
                    rel_b = relation(q(d) - q(p), c_)
                    rel_k = relation(q(kinetic[p]) - q(kinetic[d]), c_)
                    if rel_b != rel_k:
                        computed_count += 1
                        examples.append(dict(comparison="computed_order", order=order, bound=bound, primary=p,
                                             primary_index=points[p], edge=d, edge_index=edges_in.index(d),
                                             binding_relation=rel_b, kinetic_relation=rel_k))
    next_count, next_bins = 0, []
    if points:
        for d, k in sorted((high_edge[k], k) for k in range(n)):
            rel_b, rel_k = relation(u_b[k], limit), relation(u_k[k], limit)
            if rel_b != rel_k:
                next_count += 1
                examples.append(dict(comparison="next_order", order=computed + 1, bound="lower", primary=b_min,
                                     primary_index=points[b_min], edge=d, edge_index=edges_in.index(d),
                                     binding_relation=rel_b, kinetic_relation=rel_k, bin_index=k))
            if (rel_b == "<") != (rel_k == "<"):
                next_bins.append(k)
    next_bins.sort()
    warning = s_b != s_k or computed_count != 0 or next_count != 0 or bool(next_bins)
    excess = max((Fr(x) - Fr(y) for x, y in zip(width, width_kinetic)), default=Fr(0)) if points else Fr(0)
    need(excess >= 0, "the width of a bin is below its kinetic width")
    out.update(
        increasing=increasing, primary_direction=primary_direction, ke_edges=ke_edges, bin_map=bin_map,
        ke_positions=ke_positions, ke_intensities=ke_intensities, zero=zero, points=points,
        s_b=s_b, s_k=s_k, support=support, requested=requested, computed=computed, width=width,
        width_kinetic=width_kinetic, complete=complete, kinetic_reaches=kinetic_reaches,
        computed_count=computed_count, next_count=next_count, next_bins=next_bins, examples=examples,
        warning=warning, excess=excess, state="complete_on_display" if computed >= support else "truncated",
        reference_kinetic=rn(c - q(inp["reference_binding"])),
        e0_minus_reference=None if inp["e0"] is None else rn(q(inp["e0"]) - (c - q(inp["reference_binding"]))),
        rounding_max=up64(max(error.values())), exact=max(error.values()) == 0)
    if inp["kind"] == "binned":
        out["zero_binding_overlap"] = zero_loss_binding_overlap(inp)
    return out


def example_record(g, e):
    rec = dict(comparison=e["comparison"], order=e["order"], bound=e["bound"],
               primary_binding_eV=fx(e["primary"]), primary_kinetic_eV=fx(g["kinetic"][e["primary"]]),
               primary_index=e["primary_index"],
               edge_binding_eV=fx(e["edge"]), edge_kinetic_eV=fx(g["kinetic"][e["edge"]]), edge_index=e["edge_index"],
               binding_relation=e["binding_relation"], kinetic_relation=e["kinetic_relation"])
    if "bin_index" in e:
        rec["bin_index"] = e["bin_index"]
    return rec


def edge_hash(values):
    return hashlib.sha256(b"".join(struct.pack("<d", float(x)) for x in values)).hexdigest()


def geometry_record(inp, g):
    """The section of a geometry row: ``keys`` (one per R-C ledger key), the control inputs, the rest."""
    cli = inp["entry"] == "cli"
    pre = "run.json: " if cli else ""
    n = len(inp["edges"]) - 1
    keys = {
        pre + "energy_axis.conversion_constant_eV": fx(rn(g["c"])),
        pre + "energy_axis.conversion_constant_exact": fq(g["c"]),
        pre + "energy_axis.implied_offset_eV": None if g["implied"] is None else fx(rn(g["implied"])),
        pre + "energy_axis.rounding_max_eV": fx(g["rounding_max"]),
    }
    for k, e in enumerate(g["ke_edges"]):
        keys[pre + f"energy_axis.kinetic_output_edges_eV[{k}]"] = fx(e)
    zero_prefix = "CSV " if cli else ""
    for name in ["outside.low_binding"] + [f"bin[{k}]" for k in range(n)] + ["outside.high_binding"]:
        keys[zero_prefix + f"zero_loss.{name}"] = fq(g["zero"][name])
    if cli:
        keys[pre + "energy_axis.electron_energy_minus_reference_eV"] = fx(g["e0_minus_reference"])
        keys[pre + "energy_axis.support_order_binding"] = g["s_b"]
        keys[pre + "energy_axis.support_order_kinetic"] = g["s_k"]
    else:
        keys["energy_axis.primary_reference_kinetic_eV"] = fx(g["reference_kinetic"])
        if g["e0_minus_reference"] is not None:
            keys["energy_axis.electron_energy_minus_reference_eV"] = fx(g["e0_minus_reference"])
        for k, e in enumerate(g["ke_positions"]):
            keys[f"energy_axis.kinetic_primary.coordinates[{k}]"] = fx(e)
        keys["truncation.support_order"] = g["support"]
        keys["truncation.computed"] = g["computed"]
        keys["truncation.requested"] = g["requested"]
        if inp.get("tolerance") is None:
            keys["truncation.state"] = g["state"]
        keys["energy_axis.support_order_binding"] = g["s_b"]
        keys["energy_axis.support_order_kinetic"] = g["s_k"]
        for k in range(n):
            keys[f"truncation.bound_class[bin[{k}]]"] = "complete" if g["complete"][k] else "reached"
            keys[f"energy_axis.region_loss_upper_eV[bin[{k}]]"] = None if g["width"][k] is None else fx(g["width"][k])
        keys["energy_axis.support_comparison.computed_relation_mismatch_count"] = g["computed_count"]
        keys["energy_axis.support_comparison.next_relation_mismatch_count"] = g["next_count"]
        keys["energy_axis.support_comparison.next_reachability_mismatch_bins"] = g["next_bins"]
        keys["energy_axis.support_comparison.first_mismatches"] = [example_record(g, e) for e in g["examples"][:8]]
        keys["warning axis_conversion_support_disagreement present or absent"] = g["warning"]
    roles = coordinates(inp)
    record = dict(
        kind="geometry", keys=keys,
        control=dict(
            kinetic_output_edges_eV=[fx(e) for e in g["ke_edges"]],
            kinetic_primary=dict(kind=inp["kind"], coordinates=[fx(e) for e in g["ke_positions"]],
                                 intensities=[fx(s) for s in g["ke_intensities"]]),
            max_order="display_support" if g["requested"] == "display_support" and g["support"] == 0
            else g["computed"]),
        map=dict(
            direction=dict(output_edges="increasing" if g["increasing"] else "decreasing",
                           primary=g["primary_direction"]),
            kinetic_bin_of_binding_bin=g["bin_map"],
            outside={"low_binding": "high", "high_binding": "low"}),
        support=dict(
            support_order_binding=g["s_b"], support_order_kinetic=g["s_k"], support_order=g["support"],
            computed=g["computed"], state_without_tolerance=g["state"],
            region_loss_upper_eV=[None if w is None else fx(w) for w in g["width"]],
            kinetic_width_upper_eV=[None if w is None else fx(w) for w in g["width_kinetic"]],
            width_equals_kinetic=[x == y for x, y in zip(g["width"], g["width_kinetic"])],
            width_excess_over_kinetic_eV=fx(rn(g["excess"])), width_excess_over_kinetic_exact=fq(g["excess"]),
            bound_class=["complete" if x else "reached" for x in g["complete"]],
            kinetic_side_reaches_next_order=g["kinetic_reaches"],
            computed_relation_mismatch_count=g["computed_count"], next_relation_mismatch_count=g["next_count"],
            next_reachability_mismatch_bins=g["next_bins"],
            mismatch_examples_total=len(g["examples"]), warning=g["warning"]),
        conversion=dict(
            exact=g["exact"], primary_reference_kinetic_eV=fx(g["reference_kinetic"]),
            coordinates=[dict(binding_eV=fx(b_), kinetic_eV=fx(g["kinetic"][b_]), rounding_error_eV=fq(g["error"][b_]),
                              roles=roles[b_]) for b_ in sorted(roles)]))
    if "zero_binding_overlap" in g:
        differs = any(g["zero_binding_overlap"][k] != g["zero"][k] for k in g["zero"])
        record["zero_loss_exact_binding_overlap_diagnostic_only"] = dict(
            differs_from_kinetic_widths=differs, values={k: fq(v) for k, v in sorted(g["zero_binding_overlap"].items())})
    if cli:
        record["hashes"] = dict(binding_output_edges_sha256_float64_le=edge_hash(inp["edges"]),
                                kinetic_output_edges_sha256_float64_le=edge_hash(g["ke_edges"]))
        record["support"]["first_mismatches"] = [example_record(g, e) for e in g["examples"][:8]]
    return record


def tolerance_of(row):
    return row["settings"].get("neglect_tolerance_fraction")


# ------------------------------------------------------------------ known wrong answers
def within_membership(wrong, right):
    """The tolerance of prereg 7.2 for the zero loss: exact zero, else |y - r| <= 4 eps |r|."""
    return wrong == right if right == 0 else abs(Fr(wrong) - Fr(right)) <= 4*EPS*abs(Fr(right))


def wrong_answers(table, rows, analysed):
    out = []

    def add(ident, row, key, kind, wrong, right, note):
        if kind == "exact_conversion":
            outside = wrong is None or right is None or float(wrong).hex() != float(right).hex()
            w, r = (None if wrong is None else fx(wrong)), (None if right is None else fx(right))
        elif kind == "bin_membership":
            outside = not within_membership(wrong, right)
            w, r = fq(wrong), fq(right)
        else:
            outside = wrong != right
            w, r = wrong, right
        need(outside, f"known wrong answer {ident} {row} {key} is inside the tolerance")
        out.append(dict(id=ident, row=row, key=key, kind=kind, wrong=w, right=r, outside_tolerance=True, note=note))

    def inputs(rid):
        return row_inputs(table, rows[rid])

    # X34 (BE8): rounding C first changes all six KE edges
    inp, g = inputs("BE8"), analysed["BE8"]
    twice = analyse(inp, kinetic_of=lambda c, b: rn(Fr(rn(c)) - q(b)))
    need(not twice["collision"] and not twice["negative"], "BE8 rounded twice")
    need(len(g["ke_edges"]) == 6, "BE8 has six edges")
    for k in range(6):
        add("X34", "BE8", f"energy_axis.kinetic_output_edges_eV[{k}]", "exact_conversion",
            twice["ke_edges"][k], g["ke_edges"][k], "RN(RN(C) - B) instead of RN(C - B)")
    c_rounded = Fr(rn(g["c"]))
    second = up64(max(abs(Fr(rn(c_rounded - q(b))) - (c_rounded - q(b))) for b in g["kinetic"]))
    add("X34", "BE8", "energy_axis.rounding_max_eV", "exact_conversion", second, g["rounding_max"],
        "the error of the second rounding only, measured from RN(C) - B")
    add("X34", "BE8", "energy_axis.rounding_max_eV", "exact_conversion", twice["rounding_max"], g["rounding_max"],
        "the twice-rounded coordinates measured from the exact C - B")

    # X25: three wrong conversions
    inp, g = inputs("BE1"), analysed["BE1"]
    c = g["c"]
    minus = {b: rn(q(b) - c) for b in g["kinetic"]}
    need(all(v < 0 for v in minus.values()), "X25 branch 1: every KE of BE1 is negative")
    out.append(dict(id="X25(1)", row="BE1", key="terminal_state", kind="state", wrong="invalid_energy_axis/negative_kinetic_energy",
                    right="success", outside_tolerance=True, note="E = B - C: all coordinates negative"))
    plus = sorted(rn(c + q(b)) for b in inp["edges"])
    for k, (w, r) in enumerate(zip(plus, g["ke_edges"])):
        add("X25(2)", "BE1", f"energy_axis.kinetic_output_edges_eV[{k}]", "exact_conversion", w, r, "E = C + B")
    inp, g = inputs("BE7"), analysed["BE7"]
    dropped = analyse(inp, constant=g["hv"])
    add("X25(3)", "BE7", "energy_axis.conversion_constant_eV", "exact_conversion", rn(dropped["c"]), rn(g["c"]),
        "C = h nu (the offset dropped)")
    for k, (w, r) in enumerate(zip(dropped["ke_edges"], g["ke_edges"])):
        add("X25(3)", "BE7", f"energy_axis.kinetic_output_edges_eV[{k}]", "exact_conversion", w, r, "C = h nu")
    add("X25(3)", "BE7", "energy_axis.electron_energy_minus_reference_eV", "exact_conversion",
        dropped["e0_minus_reference"], g["e0_minus_reference"], "C = h nu")

    # X36 and X41 (paired anchor)
    inp, g = inputs("BE2"), analysed["BE2"]
    ref = inp["axis"]["reference"]
    sign = analyse(inp, constant=q(read(ref["kinetic_anchor_eV"])) - q(read(ref["binding_anchor_eV"])))
    add("X36", "BE2", "energy_axis.conversion_constant_eV", "exact_conversion", rn(sign["c"]), rn(g["c"]),
        "C = E_anchor - B_anchor")
    need(rn(sign["c"]) == 1314.25, "X36: the wrong constant is 1314.25")
    for k, (w, r) in enumerate(zip(sign["ke_edges"], g["ke_edges"])):
        add("X36", "BE2", f"energy_axis.kinetic_output_edges_eV[{k}]", "exact_conversion", w, r, "C = E_anchor - B_anchor")
    for rid in ("BE2", "BC02"):
        inp, g = inputs(rid), analysed[rid]
        pre = "run.json: " if rid == "BC02" else ""
        used = analyse(inp, constant=g["hv"])
        add("X41", rid, pre + "energy_axis.conversion_constant_eV", "exact_conversion", rn(used["c"]), rn(g["c"]),
            "C = h nu (the implied offset used as a correction)")
        for k, (w, r) in enumerate(zip(used["ke_edges"], g["ke_edges"])):
            add("X41", rid, pre + f"energy_axis.kinetic_output_edges_eV[{k}]", "exact_conversion", w, r, "C = h nu")

    # X32 (BS1): the support judged on the rounded KE differences only
    g = analysed["BS1"]
    need((g["s_b"], g["s_k"], g["support"], g["computed"]) == (2, 1, 2, 1), "BS1 is (2, 1, 2, 1)")
    limit = (g["computed"] + 1)*q(inputs("BS1")["window"][0])
    kinetic_only_complete = Fr(g["width_kinetic"][0]) < limit
    need(kinetic_only_complete and not g["complete"][0], "BS1: complete on the KE side only")
    add("X32", "BS1", "truncation.support_order", "truncation", g["s_k"], g["support"], "KE side only")
    add("X32", "BS1", "truncation.state", "truncation", "complete_on_display", g["state"], "KE side only (S_K = 1 = N)")
    add("X32", "BS1", "truncation.bound_class[bin[0]]", "truncation", "complete", "reached", "KE side only")
    add("X32", "BS1", "energy_axis.support_order_binding", "truncation", g["s_k"], g["s_b"], "KE side only")

    # BR1: the overlap on the exact BE widths is not the registered value
    g = analysed["BR1"]
    need(g["zero"]["bin[0]"] == Fr(1, 2) and g["zero"]["bin[1]"] == Fr(1, 2), "BR1: 1/2 and 1/2")
    need(g["zero_binding_overlap"]["bin[0]"] == Fr(3, 8) and g["zero_binding_overlap"]["bin[1]"] == Fr(5, 8),
         "BR1: 3/8 and 5/8 on the exact BE widths")
    for name in ("bin[0]", "bin[1]"):
        add("BR1-binding-overlap", "BR1", f"zero_loss.{name}", "bin_membership", g["zero_binding_overlap"][name],
            g["zero"][name], "exact BE overlap instead of the rounded-KE widths")

    # X27 and X43 (BE4), X26 branch 2 (BE1): membership
    g = analysed["BE4"]
    need([g["zero"][f"bin[{k}]"] for k in range(3)] == [Fr(1, 4), Fr(1, 2), Fr(1, 4)], "BE4: 1/4, 1/2, 1/4")
    add("X27", "BE4", "zero_loss.bin[1]", "bin_membership", Fr(0), g["zero"]["bin[1]"], "interior delta to the higher BE bin")
    add("X27", "BE4", "zero_loss.bin[2]", "bin_membership", Fr(3, 4), g["zero"]["bin[2]"], "interior delta to the higher BE bin")
    add("X43(1)", "BE4", "zero_loss.bin[0]", "bin_membership", Fr(0), g["zero"]["bin[0]"], "delta on the lowest edge to the outside")
    add("X43(1)", "BE4", "zero_loss.outside.low_binding", "bin_membership", Fr(1, 4), g["zero"]["outside.low_binding"],
        "delta on the lowest edge to the outside")
    add("X43(2)", "BE4", "zero_loss.bin[2]", "bin_membership", Fr(0), g["zero"]["bin[2]"], "delta on the highest edge to the outside")
    add("X43(2)", "BE4", "zero_loss.outside.high_binding", "bin_membership", Fr(1, 4), g["zero"]["outside.high_binding"],
        "delta on the highest edge to the outside")
    g = analysed["BE1"]
    need(g["zero"]["bin[0]"] == 1 and g["zero"]["outside.low_binding"] == Fr(1, 2), "BE1: bin[0] 1, low_binding 1/2")
    add("X27", "BE1", "zero_loss.bin[0]", "bin_membership", Fr(0), g["zero"]["bin[0]"], "interior delta to the higher BE bin")
    add("X26(2)", "BE1", "zero_loss.outside.low_binding", "bin_membership", g["zero"]["outside.high_binding"],
        g["zero"]["outside.low_binding"], "outside names not exchanged")
    add("X26(1)", "BE1", "zero_loss.bin[0]", "bin_membership", g["zero"]["bin[6]"], g["zero"]["bin[0]"], "arrays not reordered")
    return out


def self_checks(table, rows, analysed, collisions):
    done = []

    def check(name, condition):
        need(condition, f"self check failed: {name}")
        done.append(name)

    # AX_D and AX_DA give the same constant and the same KE
    c_d = axis_constant(table["axes"]["AX_D"])[0]
    c_da, hv_da, implied = axis_constant(table["axes"]["AX_DA"])
    check("AX_D and AX_DA: same C (4997/4) and implied offset 17/4", c_d == c_da == Fr(4997, 4) and implied == Fr(17, 4))
    a, b = analysed["BC01"], analysed["BC02"]
    check("AX_D (BC01) and AX_DA (BC02): the same ascending KE edges",
          [x.hex() for x in a["ke_edges"]] == [x.hex() for x in b["ke_edges"]])
    check("BE1 and BE1r: the same ascending KE edges and KE primary; bin maps reversed",
          analysed["BE1"]["ke_edges"] == analysed["BE1r"]["ke_edges"]
          and analysed["BE1"]["ke_positions"] == analysed["BE1r"]["ke_positions"]
          and analysed["BE1"]["ke_intensities"] == analysed["BE1r"]["ke_intensities"]
          and analysed["BE1"]["bin_map"] == analysed["BE1r"]["bin_map"][::-1])
    check("BE3 and BE3r: the same KE primary and zero loss",
          analysed["BE3"]["ke_positions"] == analysed["BE3r"]["ke_positions"]
          and analysed["BE3"]["ke_intensities"] == analysed["BE3r"]["ke_intensities"]
          and analysed["BE3"]["zero"] == analysed["BE3r"]["zero"])
    # design section 5 counterexample and its relatives
    check("design 5 counterexample (BS1): (S_B, S_K) = (2, 1)", (analysed["BS1"]["s_b"], analysed["BS1"]["s_k"]) == (2, 1))
    check("BS2: (1, 2)", (analysed["BS2"]["s_b"], analysed["BS2"]["s_k"]) == (1, 2))
    check("BS3: (8, 8), only bin[0] differs in reachability",
          (analysed["BS3"]["s_b"], analysed["BS3"]["s_k"]) == (8, 8) and analysed["BS3"]["next_bins"] == [0]
          and analysed["BS3"]["computed_count"] == 0)
    check("BS3n2: one computed-order mismatch (order 2, lower)",
          analysed["BS3n2"]["computed_count"] == 1 and analysed["BS3n2"]["examples"][0]["order"] == 2
          and analysed["BS3n2"]["examples"][0]["bound"] == "lower" and analysed["BS3n2"]["next_count"] == 0)
    check("BS5: no disagreement but a wider BE width", not analysed["BS5"]["warning"] and analysed["BS5"]["excess"] > 0)
    for exact_row, rounded_row in (("BS4", "BS1"), ("BS6", "BS5")):
        check(f"{exact_row}: exact conversion with the same upper width as {rounded_row}",
              analysed[exact_row]["exact"] and analysed[exact_row]["width"] == analysed[rounded_row]["width"]
              and analysed[exact_row]["excess"] == 0)
    check("BS4 bin[0] width equals BS3 bin[0] width", analysed["BS4"]["width"][0] == analysed["BS3"]["width"][0])
    check("BT1a, BT1b, BT1c: S = 3, 2, 3 and exact conversion",
          [analysed[r]["support"] for r in ("BT1a", "BT1b", "BT1c")] == [3, 2, 3]
          and all(analysed[r]["exact"] for r in ("BT1a", "BT1b", "BT1c")))
    check("BZ0: empty support", not analysed["BZ0"]["points"] and analysed["BZ0"]["support"] == 0)
    check("BE7: C - B_ref = 999 and E0 minus it is exactly 1",
          analysed["BE7"]["reference_kinetic"] == 999.0 and analysed["BE7"]["e0_minus_reference"] == 1.0)
    check("AX_R: C is not a binary64 and RN(C) = 1482", analysed["BE8"]["c"] != Fr(rn(analysed["BE8"]["c"]))
          and rn(analysed["BE8"]["c"]) == 1482.0)
    check("BE5a: delta and edge 73.1 are one coordinate; rounding_max > 0",
          analysed["BE5a"]["rounding_max"] > 0 and not analysed["BE5a"]["exact"])
    # design section 7 counterexample
    for rid in COLLISION_ROWS:
        col = collisions[rid]
        check(f"design 7 counterexample ({rid}): (1, 1 + 2^-45) -> (8191, 8191)",
              col["binding"] == (1.0, 1.0 + 2.0**-45) and col["kinetic"] == (8191.0, 8191.0)
              and col["roles"] == ("output_edges_eV[1]", "primary.binding_energies_eV[0]"))
    for rid in ("BN01", "BN23"):
        inp = row_inputs(table, rows[rid])
        c = axis_constant(inp["axis"])[0]
        check(f"{rid}: each array alone is injective (the union is not)",
              all(len({rn(c - q(b)) for b in arr}) == len(set(arr)) for arr in (inp["positions"], inp["edges"])))
    # exact conversion: the two zero-loss definitions agree
    for rid, g in sorted(analysed.items()):
        if "zero_binding_overlap" in g and g["exact"]:
            check(f"{rid}: exact conversion, the BE overlap equals the KE-width value",
                  all(g["zero_binding_overlap"][k] == g["zero"][k] for k in g["zero"]))
    for rid, g in sorted(analysed.items()):
        inp = row_inputs(table, rows[rid])
        need(sum(g["zero"].values()) == sum(q(s) for s in inp["intensities"]), f"{rid}: the zero loss sums to the primary")
    done.append("every geometry row: the zero loss of all regions sums to the primary integral")
    return done


# ------------------------------------------------------------------ build
def build():
    table = load_table()
    rows = {r["id"]: r for r in table["rows"]}
    out_rows, analysed, collisions, counts = {}, {}, {}, {}
    for rid in GEOMETRY_ROWS:
        inp = row_inputs(table, rows[rid])
        inp["tolerance"] = tolerance_of(rows[rid])
        g = analyse(inp)
        need(not g["negative"] and g["collision"] is None, f"{rid}: the conversion is not valid")
        analysed[rid] = g
        out_rows[rid] = geometry_record(inp, g)
    for rid in COLLISION_ROWS:
        inp = row_inputs(table, rows[rid])
        g = analyse(inp)
        need(g["collision"] is not None and not g["negative"], f"{rid}: no collision")
        collisions[rid] = g["collision"]
        col = g["collision"]
        out_rows[rid] = dict(kind="collision", collision=dict(
            binding_eV=[fx(x) for x in col["binding"]], kinetic_eV=[fx(x) for x in col["kinetic"]],
            roles=list(col["roles"]), conversion_constant_exact=fq(g["c"])), keys={})
    for rid in SUPPORT_ROWS:
        inp = row_inputs(table, rows[rid])
        g = analyse(inp)
        need(not g["negative"] and g["collision"] is None, f"{rid}: the conversion is not valid")
        pair = (g["s_b"], g["s_k"])
        if rid == "BT2b":
            keys = {f"R-C support orders ({pair[0]}, {pair[1]}) are the premise of the KE comparison": list(pair)}
        elif rid == "BS1x":
            keys = {f"required_order_estimate == R-C support order S = {g['support']}": g["support"]}
        else:
            keys = {"requested_order == R-C max(S_B, S_K) (the two orders are not read from the exception)": max(pair)}
        full = geometry_record(inp, g)
        out_rows[rid] = dict(kind="support_orders", keys=keys, control=full["control"], map=full["map"], support=dict(
            support_order_binding=g["s_b"], support_order_kinetic=g["s_k"], support_order=g["support"]))
    for rid, rec in out_rows.items():
        counts[rid] = len(rec["keys"])
    checks = self_checks(table, rows, analysed, collisions)
    wrong = wrong_answers(table, rows, analysed)
    return {
        "schema": SCHEMA,
        "prereg": "docs/M4B2_PREREG.md section 6",
        "table": {"file": "benchmarks/numerics/m4b2/expected_states_m4b2.json", "sha256": TABLE_SHA256},
        "method": "exact rationals (fractions.Fraction); one rounding float(Fraction) per coordinate; "
                  "the package xpsloss is not imported",
        "number_forms": {"binary64": "{'hex': float.hex()} (positive zero as 0x0.0p+0)",
                         "rational": "'numerator/positive denominator'"},
        "tolerance_of_reference": "every key is an exact rational or a binary64: u_r = 0",
        "no_reference": [],
        "row_groups": {"geometry": list(GEOMETRY_ROWS), "collision": list(COLLISION_ROWS),
                       "support_orders": list(SUPPORT_ROWS)},
        "key_counts": counts,
        "key_count_total": sum(counts.values()),
        "self_checks": checks,
        "known_wrong_answers": wrong,
        "rows": out_rows,
    }


def serialise(obj):
    text = json.dumps(obj, indent=1, sort_keys=True, ensure_ascii=True) + "\n"
    return text.encode("ascii")


def main(argv):
    if argv not in ([], ["--check"]):
        print("usage: generate_references_m4b2.py [--check]", file=sys.stderr)
        return 2
    data = serialise(build())
    if argv == ["--check"]:
        if not OUT_PATH.is_file() or OUT_PATH.read_bytes() != data:
            print("references_m4b2.json differs from the regenerated bytes", file=sys.stderr)
            return 1
        print(f"references_m4b2.json reproduced byte for byte; sha256 {hashlib.sha256(data).hexdigest()}")
        return 0
    OUT_PATH.write_bytes(data)
    print(f"wrote {OUT_PATH.name}; sha256 {hashlib.sha256(data).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
