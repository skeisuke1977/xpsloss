"""Display support on the original BE and rounded KE axes (M4b-2 §2.3)."""

from bisect import bisect_left, bisect_right
from fractions import Fraction as Q
from itertools import pairwise
from types import SimpleNamespace

from ._fold import _up
from ._single_loss_types import COMPONENTS


def _relation(value, threshold):
    return "<" if value < threshold else ">" if value > threshold else "="


def _partitions(axis, threshold):
    left = bisect_left(axis, threshold)
    return left, bisect_right(axis, threshold, lo=left)


def _relation_at(i, cuts):
    return "<" if i < cuts[0] else "=" if i < cuts[1] else ">"


class BindingDisplay:
    def __init__(self, info, points, image):
        self.info = info
        self.points = points  # (binding coordinate, kinetic coordinate, original primary index), increasing BE
        self.increasing = info["direction"]["output_edges"] == "increasing"
        edges = info["input_binding_edges_eV"]
        indices = range(len(edges)) if self.increasing else reversed(range(len(edges)))
        self.edges = [(edges[i], image[edges[i]], i) for i in indices]

    def support_order(self, lower):
        # Rows correspond to the high-BE edge of each bin, in increasing BE order.
        self.widths = []
        if self.points:
            b0, e0, _ = self.points[0]
            self.widths = [(Q(b) - Q(b0), Q(e0) - Q(e)) for b, e, _ in self.edges[1:]]
        sb, sk = (max(0, w // Q(lower)) for w in self.widths[-1]) if self.widths else (0, 0)
        widths = [max(pair) for pair in self.widths]
        upper = [_up(w) for w in widths] if widths else [None] * (len(self.edges) - 1)
        self.info.update(support_order_binding=sb, support_order_kinetic=sk,
                         region_loss_upper_eV=upper if self.increasing else upper[::-1])
        # truncation() subtracts each left edge from pmax, in ascending KE order.
        # Exact synthetic edges preserve strict comparisons. The exponential uses
        # the existing upward rounding inside bound(), after testing reachability.
        self.bound_geometry = SimpleNamespace(
            points=(0,), edges=tuple(-w for w in reversed(widths)) + (0,) if widths else (0,) * len(self.edges)
        )
        return max(sb, sk)

    def compare(self, window, n):
        examples = []
        reachability = []
        computed = next_count = 0
        if self.points:
            # Both edge lists increase. In each, the edges below and equal
            # to a threshold are contiguous; four cuts suffice for the count.
            b_edges = [Q(b) for b, _, _ in self.edges]
            neg_e_edges = [-Q(e) for _, e, _ in self.edges]
            for order in range(1, n + 1):
                for name, threshold in zip(("lower", "upper"), window):
                    loss = order * Q(threshold)
                    for point in self.points:
                        b, e, _ = point
                        bc = _partitions(b_edges, Q(b) + loss)
                        kc = _partitions(neg_e_edges, loss - Q(e))
                        if bc == kc:
                            continue
                        for lo, hi in pairwise(sorted({*bc, *kc})):
                            rb, rk = _relation_at(lo, bc), _relation_at(lo, kc)
                            if rb != rk:
                                computed += hi - lo
                                for i in range(lo, min(hi, lo + 8 - len(examples))):
                                    examples.append(self._example("computed_order", order, name, point, i, rb, rk))

            threshold = (n + 1) * Q(window[0])
            for i, (ub, uk) in enumerate(self.widths):
                rb, rk = _relation(ub, threshold), _relation(uk, threshold)
                if rb != rk:
                    next_count += 1
                    index = i if self.increasing else len(self.widths) - 1 - i
                    if (ub >= threshold) != (uk >= threshold):
                        reachability.append(index)
                    if len(examples) < 8:
                        example = self._example("next_order", n + 1, "lower", self.points[0], i + 1, rb, rk)
                        examples.append(dict(example, bin_index=index))
        if not self.increasing:
            reachability.reverse()
        self.info["support_comparison"] = dict(computed_relation_mismatch_count=computed,
                                               next_relation_mismatch_count=next_count,
                                               next_reachability_mismatch_bins=reachability,
                                               first_mismatches=[dict(sorted(e.items())) for e in examples])

    def _example(self, kind, order, name, point, i, rb, rk):
        b, e, index = point
        bd, ed, edge_index = self.edges[i]
        return dict(comparison=kind, order=order, bound=name, primary_binding_eV=b, primary_kinetic_eV=e,
                    primary_index=index, edge_binding_eV=bd, edge_kinetic_eV=ed, edge_index=edge_index,
                    binding_relation=rb, kinetic_relation=rk)

    def orient_truncation(self, trunc):
        if self.increasing:
            for k in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
                trunc[k].reverse()

    def orient_result(self, result):
        self.orient_truncation(result.truncation)
        if self.increasing:
            for key in ("zero_loss", "zero_loss_error", "spectrum", "spectrum_error"):
                setattr(result, key, getattr(result, key)[::-1].copy())
        sides = {"low_binding": "high", "high_binding": "low"}
        result.outside = {
            side: {**{k: result.single_loss.outside[ke][k] for k in ("zero_loss", "zero_loss_error")},
                   **result.outside[ke]} for side, ke in sides.items()
        }
        for order in result.orders:
            order.outside = {side: order.outside[ke] for side, ke in sides.items()}
            if self.increasing:
                for c in COMPONENTS:
                    for key in (c, c + "_error"):
                        setattr(order, key, getattr(order, key)[::-1].copy())
                    parts = getattr(order, c + "_error_parts")
                    for key in parts:
                        parts[key] = parts[key][::-1].copy()
        result.energy_axis = self.info
        warnings = result.diagnostics["warnings"]

        def warn(code, message):
            warnings.append(dict(source="xps_multiple_loss", code=code, message=message))

        warn("binding_axis_assumed_calibrated", "Input binding energies are assumed calibrated to the Fermi level.")
        if self.info["not_recorded"]:
            warn("energy_axis_field_not_recorded", "Measurement fields not recorded: " + ", ".join(self.info["not_recorded"]))
        angle = self.info["record"]["emission_angle_deg"]
        if angle != "not_recorded" and angle != result.single_loss.theta_deg:
            warn("recorded_emission_angle_differs_from_theta", "The recorded emission angle differs from the model angle.")
        counts = self.info["support_comparison"]
        if (self.info["support_order_binding"] != self.info["support_order_kinetic"]
                or counts["computed_relation_mismatch_count"] or counts["next_relation_mismatch_count"]):
            warn("axis_conversion_support_disagreement", "Binding and rounded kinetic coordinates give different support comparisons.")
        result.warnings = tuple(w["message"] for w in warnings)
