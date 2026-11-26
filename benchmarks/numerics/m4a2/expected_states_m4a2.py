"""Build the M4a-2 preregistration ledger from frozen inputs (standard library only)."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "expected_states_m4a2.json"
BASE = "e4e024fb3a4abbca078306c93e39c1f37d1cb845"
FILES = {
    "expected_states_m4a.json": "250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408",
    "references_rb_m4a.json": "942353c748b4959220ba4da15f078a8ebdbdd9833d3abdc2e80cb8acfd6e8d94",
}
AREA = "api_area"
SPEC = "api_spectrum"
SEED = "weight_length_seed"
PROD = "production"
IE = "SingleLossInputError"
CE = "SingleLossConvergenceError"
RB = [f"RB{i}" for i in range(1, 9)] + [f"RB{i}a" for i in range(1, 5)]
SEED_ROWS = RB + [f"CE{i}" for i in range(1, 7)]


def frozen(name):
    data = (HERE.parent / "m4a" / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != FILES[name]:
        raise ValueError(f"unexpected frozen input: {name}")
    return json.loads(data)


def failure(kind, reason, phase=None, stage=None, **details):
    status = "not_converged" if kind == CE else "error"
    return dict(status=status, exception=kind, reason=reason, phase=phase, stage=stage, **details)


def ledger():
    old = frozen("expected_states_m4a.json")
    references = frozen("references_rb_m4a.json")["rows"]
    old_rows = {row["id"]: row for row in old["rows"]}
    rows = []

    def add(row_id, entry, fixture, rule, keys, *, ci=True, expect=None, **params):
        rows.append(dict(id=row_id, entry=entry, inputs=dict(fixture=fixture, **params),
                         expect={"status": "success"} if expect is None else expect,
                         keys=keys.split() if isinstance(keys, str) else keys, rule=rule, ci=ci))

    for entry in (AREA, SPEC, "cli_area", "cli_spectrum"):
        modes = ["omit", "equal_external"]
        if entry.startswith("api"):
            modes.insert(1, "explicit_weight_length")
        for mode in modes:
            add(f"V1-{entry}-{mode}", entry, "short", "paired_baseline", "legacy new_fields", mode=mode)
    for rid in ("RB1", "RB1a", "RB3", "RB3a"):
        add("V1-"+rid, old_rows[rid]["entry"], rid, "paired_baseline", "legacy new_fields",
            ci=False, mode="equal_external")
    for rid in old["execution_order"]:
        row = old_rows[rid]
        expected = dict(m4a_state=rid, status=row["expect"]["status"])
        keys = "original_check_plan"
        params = {}
        if rid == "N25":
            expected.update(phase=None, stage=None)
            keys += " exception.phase exception.stage"
            params["modes"] = ["omit", "equal_external"]
        add("V1-old-"+rid, row["entry"], rid, "m4a_ledger", keys,
            ci=row["ci"], expect=expected, **params)
    for entry in (AREA, SPEC):
        for convention in ("external", "inverse_window_area", "weight_length"):
            add(f"V2-{entry}-{convention}", entry, "short", "beta_replay", "beta_numeric beta_diag",
                multiplier=convention, c=0.03125, external_m_A=7 if convention == "external" else None)
    add("V2-area-only", AREA, "short", "area_invariance", "area11 area_error11 area_budget11 area_parts11",
        weight_length_A=20, multiplier_A=7)
    add("V2-real", SPEC, "RB1", "beta_replay", "beta_numeric beta_diag", ci=False,
        multiplier="inverse_window_area", c=0)
    add("V2-real-weight_length", AREA, "RB1a", "beta_replay",
        "beta_numeric beta_diag returned_errors rho", ci=False, derived=True, multiplier="weight_length", c=0)

    rb_key_counts = Counter()
    for rid in RB:
        refkeys = sorted(references[rid]["keys"])
        # Every bin has five densities; outside regions have integrals only.
        densities = [k+"_per_eV" for k in refkeys if k.startswith("bin[") and k.count(".") == 1]
        rb_key_counts.update(references[rid]["keys"][k]["kind"] for k in refkeys)
        rb_key_counts["derived_density"] += len(densities)
        keys = refkeys + densities
        add("V3-"+rid, old_rows[rid]["entry"], rid, "external_map", keys,
            ci=rid == "RB1", multiplier="external", m_over_Lref=2)
        add("V4-"+rid, old_rows[rid]["entry"], rid, "s2_map",
            keys + ["multiplier_A", "multiplier_interval_A", "returned_errors", "rho", "denominator"],
            ci=rid == "RB1", multiplier="inverse_window_area")

    area_conditions = {}
    for rid in SEED_ROWS:
        if "window_area.term_13a" not in references[rid]["keys"]:
            raise ValueError(f"{rid}: 13a key missing")
        cond = references[rid]["conditions"]
        condition = (cond["gap_eV"], tuple(map(tuple, cond["oscillators_hex"])), tuple(cond["loss_window_eV"]),
                     cond["epsilon_momentum"], cond["electron_energy_eV"])
        area_conditions.setdefault(condition, []).append(rid)
        for c_name, c in (("z", 0), ("p", 0.03125)):
            for seed in (10, 20, 40):
                add(f"V5-{rid}-{c_name}-{seed}", "seed_test_entry", rid, "seed_reference",
                    "seed.area13 seed.error13 weight_length_A weight_length_interval_A",
                    ci=(rid, c_name, seed) == ("RB1a", "z", 20), c=c, seed_A=seed)
            add(f"V5-{rid}-{c_name}-intersection", "test_entry", rid, "seed_intersection",
                "common_interval", ci=False, c=c, seeds_A=[10, 20, 40])

    add("V3-layer", "integrate_window", "constant", "analytic_map", "regions5 densities Hpm H LbK_b w",
        external_m_A=4)
    add("V3-tails", "integrate_window", "tails", "direct_tails", "outside.low outside.high")
    add("V4-stages", "test_driver", "two_stages", "refold",
        "stage_m retained_panel.folded_values retained_panel.parts5 budget balance")
    add("V4-error", "test_finish", "asymmetric", "error_formula",
        "rho signed_error zero_value_error spectrum_error zero_loss_error density_error w_error "
        "satellite_errors outside_errors H_errors Hpm_errors upward_rounding")
    add("V6-wire", "_run+kernel_spy", "length_kernel", "wiring",
        "alpha.call_kwargs beta.call_kwargs seed.area13 window_area.term_13b kernel_calls LbK_b")
    for position in ("inside", "equal", "outside"):
        expected = failure(CE, "weight_length_inconsistent", PROD, 1) if position == "outside" else None
        add("V6-"+position, "test_driver", "consistency", "consistency",
            "consistency.delta consistency.allowance", expect=expected, position=position)
    for a in (1, 2):
        add(f"V7-a{a}", "integrate_window", "constant", "partial_sensitivity",
            "w H Hpm LbK_b balance", a=a, multiplier="inverse_window_area")
    add("V7-ext", "integrate_window", "constant", "partial_sensitivity", "H Hpm LbK_b w balance",
        a=2, external_m_A=4)
    for material, window, m, length, weight in (
        ("Al", [1, 50], 19.4, "27.0571259020694", "1.0861615888385654"),
        ("Si", [3, 50], 24.6, "25.549842388377478", "1.6938314906967444"),
    ):
        for seed in ("default", "historical"):
            add(f"V8-{material}-{seed}", AREA, "M3a2-"+material, "historical_replay",
                "weight_length_A single_loss_weight", ci=False, window=window, c=0,
                external_m_A=m, seed_A=None if seed == "default" else m,
                reference_length_A=length, reference_weight=weight)

    # New scalar arguments and the validation order at the public boundary.
    def reject(suffix, fixture, reason, kind=IE, **params):
        add("V9-"+suffix, params.pop("entry", AREA), fixture, "input_before_kernel",
            "exception calls", expect=failure(kind, reason, kernel_evaluations=0), **params)

    for kind in ("weight_length_kind", "multiplier"):
        for bad in ("unknown", "np_bool", "array2", "none"):
            reject(f"{kind}-{bad}", "short", "invalid_"+kind, argument=kind, value_token=bad)
    for arg, reason, fixture in (
        ("extra_inverse_length_per_A", "invalid_extra_inverse_length", "derived"),
        ("weight_length_seed_A", "invalid_weight_length_seed", "derived"),
        ("multiplier_A", "invalid_multiplier_length", "external_m"),
    ):
        bad_values = ["text_number", "bool", "np_bool", "nan", "inf", "negative", "overflow_int"]
        if arg != "extra_inverse_length_per_A":
            bad_values.append("zero")
        for bad in bad_values:
            reject(f"{arg}-{bad}", fixture, reason, argument=arg, value_token=bad)
    for bad in ("none", "text_number", "bool", "np_bool", "nan", "inf", "zero", "negative"):
        reject("old-length-"+bad, "short", None, "ValueError", argument="weight_length_A", value_token=bad)
    for suffix, fixture, arg, value, reason in (
        ("derived-length", "derived", "weight_length_A", 20, "invalid_weight_length"),
        ("external-c", "short", "extra_inverse_length_per_A", 1, "extra_inverse_length_without_derived_weight"),
        ("external-seed", "short", "weight_length_seed_A", 10, "weight_length_seed_without_derived_weight"),
        ("m-missing", "external_m", "multiplier_A", None, "invalid_multiplier_length"),
        ("m-extra", "short", "multiplier_A", 20, "invalid_multiplier_length"),
        ("s2-extra", "s2", "multiplier_A", 20, "invalid_multiplier_length"),
    ):
        reject(suffix, fixture, reason, argument=arg, value=value)
    reject("order-kind", "short", "invalid_weight_length_kind", weight_length_kind="bad", weight_length_A=-1)
    reject("order-old", "short", None, "ValueError", weight_length_A=-1, multiplier="bad")
    reject("order-window", "derived", "invalid_window", window=[6, 5], multiplier="bad")
    add("V9-api-zero-c", AREA, "short", "paired_baseline", "legacy", extra_inverse_length_per_A=0)
    add("V9-pre-N30", SPEC, "N30", "input_before_kernel", "exception calls seed_entries",
        expect=failure(IE, "cell_limit_exceeded", kernel_evaluations=0, seed_entries=0), derived=True)
    add("V9-pre-N35a", AREA, "N24", "input_before_kernel", "exception calls seed_entries required_panels",
        expect=failure(IE, "resolution_exceeds_panel_limit", kernel_evaluations=0, seed_entries=0, required_panels=401),
        derived=True)
    upper = float.fromhex("0x1.0000000000001p+1")
    add("V9-pre-N32a", AREA, "derived", "input_before_kernel", "exception calls seed_entries panel_eV panel_width_ulp",
        expect=failure(IE, "panel_width_unresolved", kernel_evaluations=0, seed_entries=0,
                       panel_eV=[2, upper], panel_width_ulp=1),
        window=[1, upper])
    for suffix, mode, c, conv, reason in (
        ("zero-w", "derived", 0, "weight_length", "weight_length_undefined"),
        ("zero-both", "derived", 0, "inverse_window_area", "weight_length_undefined"),
        ("zero-s2", "external", 0, "inverse_window_area", "multiplier_undefined"),
        ("zero-c-s2", "derived", 3, "inverse_window_area", "multiplier_undefined"),
    ):
        reject(suffix, "Z1", reason, entry=SPEC, weight=mode, c=c, multiplier=conv)
    for conv in ("weight_length", "external"):
        add("V9-zero-c-"+conv, SPEC, "Z1", "structural_zero", "seed_skipped kernel_calls interval w J0 J",
            c=3, multiplier=conv, external_m_A=7 if conv == "external" else None)
    add("V9-zero-overflow", AREA, "Z1", "structural_zero", "exception calls",
        expect=failure(CE, "weight_length_not_representable", kernel_evaluations=0),
        c="min_subnormal")
    add("V9-zero-primary", SPEC, "derived", "zero_primary", "area11 w w_error regions5 balance",
        primary="D1000zero", c=0.03125, multiplier="inverse_window_area")
    add("V9-signed", "_run+kernel_spy", "constant", "signed_s2",
        "window_area.hole single_loss_hole w weight_length_A multiplier_A",
        c=0.125, multiplier="inverse_window_area")

    # Both phases pass through the actual decision and refinement loop.
    for phase, prefix, tag in ((SEED, "weight_length", "a"), (PROD, "multiplier", "b")):
        for name, reason, kind, stage in (
            ("negative", "undefined", IE, 0), ("positive-limit", "unresolved", CE, 0),
            ("zero-limit", "unresolved", CE, 1), ("inverse-overflow", "not_representable", CE, 0),
            ("end-overflow", "not_representable", CE, 0),
        ):
            add(f"V9-{tag}-{name}", "test_driver", name, "denominator_failure",
                "exception denominator kernel_calls fold_calls worst_key",
                expect=failure(kind, prefix+"_"+reason, phase, stage), phase=phase)
        for name in ("release", "ties", "negative-hold"):
            add(f"V9-{tag}-{name}", "test_driver", name, "denominator_progress",
                "denominator split_panel fold_calls final_panels stage", phase=phase)
        for stage in (0, 1):
            add(f"V9-{tag}-up{stage}", "_run+kernel_spy", "upstream", "kernel_failure",
                "exception inner_attributes calls", expect=failure(CE, "kernel_not_converged", phase, stage),
                phase=phase, fail_stage=stage)
    add("V9-c-dominates", "test_driver", "small_area_plus_c", "denominator_progress",
        "denominator raw_budget weight_length_A weight_length_interval_A stage", phase=SEED)
    add("V9-exact-sum", "test_driver", "exact_sum", "denominator_progress",
        "denominator weight_length_A weight_length_interval_A", phase=SEED)
    for suffix, fixture, reason in (
        ("other-point", "budget_point", "kernel_error_exceeds_budget"),
        ("other-loss", "budget_loss", "loss_budget_not_met"),
        ("width", "refined_width", "panel_width_unresolved"),
        ("cells", "refined_cells", "cell_limit_reached"),
    ):
        # On geometry failure the new mesh has not reached the stage counter update.
        add("V9-"+suffix, "test_driver", fixture, "old_failure", "exception diagnostics.stages worst_or_geometry",
            expect=failure(CE, reason, PROD, diagnostics={"stages": 0}), phase=PROD)

    # Instrumented process-local CLI; every manifest must be replayable on its own.
    for n in range(1, 8):
        for entry in ("cli_area", "cli_spectrum"):
            add(f"V10-{n}-{entry}", entry, "short", "convention_record",
                "api stdout manifest csv replay warnings comparison units", combination=n)
    for entry in ("cli_area", "cli_spectrum"):
        add("V10-legacy-"+entry, entry, "short", "old_manifest", "legacy replay")
    for n in (1, 2):
        for match in ("partial", "absent"):
            add(f"V10-{n}-{match}", "cli_spectrum", "short", "convention_record",
                "api stdout manifest warnings comparison", combination=n, comparison_fixture=match)
    add("V10-coincidence", "test_finish", "s2_equals_weight", "convention_record",
        "model_id warnings comparison multiplier_interval_A", combination=4)
    add("V10-caps", "cli_capabilities", "none", "capabilities", "supported xps_single_loss other_blocks")
    for cls in ("LossWindowAreaResult", "SingleLossSpectrumResult"):
        add("V10-pos-"+cls, "constructor", "short", "old_constructor", "legacy new_fields", cls=cls)
    for name in ("no_weight", "both_weights", "derived_source", "external_c_zero", "external_c_source",
                 "external_seed", "missing_c_with_source", "both_m", "missing_m_with_source", "s2_with_source"):
        for entry in ("cli_area", "cli_spectrum"):
            add(f"V10-{entry}-{name}", entry, "short", "usage_error", "stderr stdout exit outdir calls",
                expect=failure("ValueError", None, exit_code=2, output_dir="absent"), flag_case=name)
    for phase, tag in ((SEED, "a"), (PROD, "b")):
        add(f"V10-stderr-{tag}", "cli_spectrum", "upstream", "kernel_failure", "stderr stdout exit outdir inner_attributes",
            expect=failure(CE, "kernel_not_converged", phase, 0, exit_code=2, output_dir="absent"),
            phase=phase, fail_stage=0)

    # A mutant counts only when its specified normal control passes and the witness fails.
    mutations = (
        ("scale_area", "V2-area-only", "window_area.total", "raw area depends on multiplier"),
        ("s2_13a", "V4-RB1", "multiplier_A", "wrong denominator after averaging"),
        ("weight_total", "V5-RB1a-z-20", "weight_length_A", "reciprocal must use the measured 13a column"),
        ("s2_add_c", "V2-api_spectrum-inverse_window_area", "cutoff_sensitivity.total", "control has identical kernel but no extra rate in S2"),
        ("weight_drop_c", "V9-c-dominates", "weight_length_A", "c dominates 13a"),
        ("interval_drop_c", "V9-c-dominates", "weight_length_interval_A", "reciprocal endpoints use A+c"),
        ("H_uses_Lw", "V7-ext", "cutoff_sensitivity.total", "fixed-m derivative"),
        ("w_uses_Lw", "V7-ext", "single_loss_weight", "area times external multiplier"),
        ("upper_uses_Lw", "V7-ext", "diagnostics.LbK_b.total", "upper endpoint uses m"),
        ("balance_uses_Lw", "V7-ext", "diagnostics.balance.single_loss_total", "target and roundoff use m"),
        ("return_seed", "V6-wire", "window_area.term_13b", "13b depends on sampled length"),
        ("sample_seed_in_beta", "V6-wire", "beta.call_kwargs.weight_length_A", "seed differs from production length"),
        ("skip_consistency", "V6-outside", "exception.reason", "first representable excess must stop"),
        ("round_inwards", "V9-zero-c-weight_length", "weight_length_interval_A", "one third is not representable"),
        ("first_m", "V4-stages", "retained_panel.folded_values", "resolved multiplier changes across stages"),
        ("omit_error", "V4-error", "signed_error", "nonzero interval and product term"),
        ("inflate_primary", "V4-error", "zero_loss_error", "primary is invariant"),
        ("stop_on_hold", "V9-b-ties", "split_panel", "all local errors zero but refinement required"),
        ("scale_spectrum", "V3-RB1", "bin[6].spectrum", "unit primary survives q=2"),
        ("false_matched", "V10-3-cli_spectrum", "comparison.spectrum", "old matching kernel is insufficient"),
        ("perturb_equal_m", "V1-api_spectrum-equal_external", "single_loss_total", "one ulp breaks legacy bits"),
        ("old_model_id", "V10-2-cli_area", "stdout.model_id", "equal external m uses legacy ID"),
        ("old_warning", "V1-api_area-omit", "warnings", "ordered legacy warnings change"),
        ("residual_low", "V3-tails", "outside.low.single_loss_total", "small direct tail survives large primary"),
        ("residual_high", "V3-tails", "outside.high.single_loss_total", "small direct tail survives large primary"),
        ("parallel_bulk", "V5-RB2a-z-20", "seed.area13", "dispersive parallel 13a differs from bulk"),
    )
    for index, (change, witness, key, why) in enumerate(mutations, 1):
        add(f"V11-M{index:02}", "mutation_harness", witness, "mutation", [key], ci=False,
            expect=dict(status="detected", control="pass", mutant="fail"),
            change=change, witness=witness, cause=why)

    ids = [row["id"] for row in rows]
    ci_ids = [row["id"] for row in rows if row["ci"]]
    reduction = ["V4-RB1", "V3-RB1"]
    if len(set(ids)) != len(ids) or not set(reduction) <= set(ci_ids):
        raise ValueError("duplicate row or invalid CI reduction")
    if set(references) != set(SEED_ROWS) or len(old_rows) != 100:
        raise ValueError("source row inventory differs from registration")
    if any(witness not in ids for _, witness, _, _ in mutations):
        raise ValueError("unknown mutation witness")
    conventions = {
        "1": dict(weight="external", length_A=20, multiplier="weight_length"),
        "2": dict(weight="external", length_A=20, multiplier="external", multiplier_A=20),
        "3": dict(weight="external", length_A=20, multiplier="external", multiplier_A=7),
        "4": dict(weight="external", length_A=20, multiplier="inverse_window_area"),
        "5": dict(weight="derived", c=0, multiplier="weight_length"),
        "6": dict(weight="derived", c=0.03125, seed_A=10, multiplier="external", multiplier_A=7),
        "7": dict(weight="derived", c=0.03125, multiplier="inverse_window_area"),
    }
    return dict(schema="xpsloss-m4a2-prereg-v1", baseline=BASE, frozen_inputs=FILES,
                conventions=conventions,
                reference_inventory=dict(rb_keys_per_map=dict(rb_key_counts),
                                         v5_distinct_area_conditions=len(area_conditions),
                                         v5_area_groups=list(area_conditions.values())),
                defaults=old["defaults"], columns=old["columns"], rules_document="docs/M4A2_PREREG.md",
                execution_order=ids, ci_reduction_order=reduction,
                limits_seconds=dict(ci_increment=60, full_unmutated=7200, mutations=7200),
                counts=dict(rows=len(ids), ci=len(ci_ids), inherited=len(old_rows), mutations=len(mutations)),
                rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify exact bytes without writing")
    args = parser.parse_args()
    table = ledger()
    content = json.dumps(table, ensure_ascii=True, allow_nan=False, indent=1) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != content.encode("utf-8"):
            parser.exit(1, "ledger differs\n")
    else:
        OUTPUT.write_text(content, encoding="utf-8", newline="\n")
    print(json.dumps(table["counts"]), hashlib.sha256(content.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
