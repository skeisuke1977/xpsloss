"""Each fault changes one production operation in memory, for one registered witness."""

from contextlib import contextmanager, ExitStack
import inspect
import math
import os
from unittest.mock import patch

from xpsloss import _single_loss_cli as cli_impl, bulk_diimfp
import m4a2_states as m

api, panels, fold = m.api, m.panels, m._fold


def replace(owner, name, *edits):
    original = getattr(owner, name)
    source = inspect.getsource(original)
    for old, new in edits:
        if source.count(old) != 1:
            raise AssertionError(f"mutation anchor {name} / {old!r} is not unique")
        source = source.replace(old, new)
    namespace = {}
    exec(compile(source, inspect.getfile(original), "exec"), vars(owner), namespace)
    return patch.object(owner, name, namespace[name])


def finish_change(change):
    original = panels._finish

    def finish(*args):
        r = original(*args)
        if change == "scale_area":
            r.window_area = {c: v * args[1] for c, v in r.window_area.items()}
        elif change == "scale_spectrum":
            r.spectrum *= 2.0
            for side in r.outside.values():
                side["spectrum"] *= 2.0
        else:
            side = "low" if change == "residual_low" else "high"
            other = "high" if side == "low" else "low"
            r.outside[side]["single_loss_total"] = r.primary_integral * r.single_loss_weight - math.fsum(
                [*r.single_loss_total, r.outside[other]["single_loss_total"]]
            )
        return r

    return patch.object(panels, "_finish", finish)


def kernel_change():
    original = api._call_kernel

    def kernel(mat, energy, loss, **kw):
        r = original(mat, energy, loss, **kw)
        if kw["epsilon_momentum"] == "k_parallel":
            vals, errs = zip(
                *(
                    bulk_diimfp(mat, energy, t, epsabs=kw["epsabs"], epsrel=kw["epsrel"], limit=kw["limit"])
                    for t in loss
                )
            )
            r.term_13a[:] = vals
            r.component_errors["term_13a"][:] = errs
        return r

    return patch.object(api, "_call_kernel", kernel)


@contextmanager
def fault(change):
    changes = {
        "s2_13a": (
            panels,
            "_integrate",
            ('if seed else "window_area.total"', 'if seed else "window_area.term_13a"'),
        ),
        "weight_total": (
            panels,
            "_integrate",
            (
                "if seed:\n                    outputs =",
                'denominator_row = outputs["window_area.total"]\n                if seed:\n                    outputs =',
            ),
            ("_denominator(outputs[denominator_key],", "_denominator(denominator_row,"),
        ),
        "s2_add_c": (
            api,
            "_run",
            ("prepared=prepared, multiplier=", "extra=extra, prepared=prepared, multiplier="),
        ),
        "weight_drop_c": (
            panels,
            "_denominator",
            ("return center, interval", "return (float(1/area) if extra and area else center), interval"),
        ),
        "interval_drop_c": (
            panels,
            "_denominator",
            (
                "return center, interval",
                "return center, (inverse_area(area, error)[1] if extra and area else interval)",
            ),
        ),
        "sample_seed_in_beta": (
            api,
            "_run",
            ('kw["weight_length_A"] = length', 'kw["weight_length_A"] = seed_length'),
        ),
        "skip_consistency": (api, "_run", ("if delta > allowance:", "if False:")),
        "first_m": (
            panels,
            "_integrate",
            ("    previous = {}", "    first_m = None\n    previous = {}"),
            (
                "folded = _fold._fold(panels, geometry, m, lower_minus)",
                "first_m = m if first_m is None else first_m\n                folded = _fold._fold(panels, geometry, first_m, lower_minus)",
            ),
        ),
        "omit_error": (
            fold,
            "propagate_multiplier",
            ("return _up(Q(error) + rho*(abs(Q(y)) + Q(error)))", "return error"),
        ),
        "inflate_primary": (
            fold,
            "propagate_multiplier",
            (
                "    for side in result.outside.values():",
                "    result.zero_loss_error = np.array([corrected(y,e) for y,e in zip(result.zero_loss, result.zero_loss_error)])\n    for side in result.outside.values():",
            ),
        ),
        "stop_on_hold": (
            panels,
            "_integrate",
            (
                "split[index] = True",
                'raise SingleLossConvergenceError(name+"_unresolved", diagnostics=diag, stage=stage)',
            ),
        ),
        "false_matched": (
            api,
            "_run",
            (
                'diag["comparison"] = _comparison(conditions, edges, diag["resolution_basis"])',
                'diag["comparison"] = _comparison({k:v for k,v in conditions.items() if k not in ("multiplier","multiplier_A","weight_length_kind")}, edges, diag["resolution_basis"])',
            ),
        ),
        "perturb_equal_m": (
            api,
            "_run",
            (
                "    return result",
                '    if legacy and multiplier == "external" and primary is not None:\n        vals = result.single_loss_total\n        vals[vals>0] = np.nextafter(vals[vals>0], np.inf)\n    return result',
            ),
        ),
        "old_model_id": (
            cli_impl,
            "run",
            (
                'return dict(status="complete", model_id=result.model_id',
                'return dict(status="complete", model_id=api.SINGLE_LOSS_CONVENTIONS_MODEL_ID',
            ),
        ),
        "old_warning": (
            api,
            "_warnings",
            (
                "    result.warnings =",
                '    add("W-sens","Fixed-length partial sensitivity.")\n    result.warnings =',
            ),
        ),
    }
    with ExitStack() as stack:
        if change in changes:
            stack.enter_context(replace(*changes[change]))
            if change == "s2_add_c":
                stack.enter_context(replace(panels, "_integrate", ("extra if seed else 0.", "extra")))
        elif change in ("scale_area", "scale_spectrum", "residual_low", "residual_high"):
            stack.enter_context(finish_change(change))
        elif change in ("H_uses_Lw", "w_uses_Lw", "upper_uses_Lw", "balance_uses_Lw"):
            target = {
                "H_uses_Lw": "h = scaled(Q(length)",
                "w_uses_Lw": "w = length*",
                "upper_uses_Lw": 'diag["LbK_b"] = scaled(Q(length)',
                "balance_uses_Lw": "expected = float(Q(mass)*Q(length)",
            }[change]
            # Supply the independently specified kernel weight from the actual driver.
            stack.enter_context(
                replace(
                    panels,
                    "_integrate",
                    (
                        "result = _finish(window, m,",
                        'diag["_mutant_weight"] = length\n                result = _finish(window, m,',
                    ),
                )
            )
            stack.enter_context(
                replace(panels, "_finish", (target, target.replace("length", 'diag["_mutant_weight"]')))
            )
        elif change == "return_seed":
            original = api._seed_length
            area = panels._areas
            values = {}

            def seed(*args):
                def observe(ps):
                    result = area(ps)
                    values.update({k.removeprefix("window_area."): v["value"] for k, v in result[0].items()})
                    return result

                with patch.object(panels, "_areas", observe):
                    return original(*args)

            stack.enter_context(patch.object(api, "_seed_length", seed))
            stack.enter_context(
                replace(
                    api,
                    "_run",
                    (
                        "    return result",
                        "    result.window_area.update(_mutant_seed_areas)\n    return result",
                    ),
                )
            )
            stack.enter_context(patch.object(api, "_mutant_seed_areas", values, create=True))
        elif change == "round_inwards":
            from xpsloss import _single_loss_types as types

            stack.enter_context(
                replace(
                    types,
                    "inverse_area",
                    ("if lo > lower:", "if lo < lower:"),
                    ("if hi < upper:", "if hi > upper:"),
                    ("nextafter(lo, -math.inf)", "nextafter(lo, math.inf)"),
                    ("nextafter(hi, math.inf)", "nextafter(hi, -math.inf)"),
                )
            )
            stack.enter_context(patch.object(panels, "inverse_area", types.inverse_area))
        elif change == "parallel_bulk":
            stack.enter_context(kernel_change())
        else:
            raise AssertionError(change)
        yield


def witness(row, change=None):
    from test_m4a2_records import pairs

    m.cached_seed.cache_clear()
    pairs.cache_clear()
    if change in ("perturb_equal_m", "old_warning"):
        with patch.dict(os.environ, XPSLOSS_M4A2_CHILD_MUTANT=change):
            return m.execute(row)
    if change:
        with fault(change):
            return m.execute(row)
    return m.execute(row)


def mutation(row, log):
    i = row["inputs"]
    target = m.ROWS[i["witness"]]
    (key,) = row["keys"]
    control = witness(target)
    altered = witness(target, i["change"])

    def verdict(run):
        # The V3 reference is also observed through the primary-only identity.
        # Credit only the registered reference when it is required.
        if i["change"] == "scale_spectrum":
            return run.items.get(("external", key), {}).get("status", "not_run")
        return run.status(key)

    a, b = verdict(control), verdict(altered)
    errors = [run.items.get(("execution", "unexpected_exception")) for run in (control, altered)]
    if any(errors):
        log.record(
            "mutation", key, "fail", control=a, mutant=b, reason="witness_execution_failed", errors=errors
        )
        return
    if a in ("not_run", "no_reference", "unresolved", "reference_path_failure"):
        extra = {}
        if i["change"] == "scale_spectrum":
            extra["primary_only_identity"] = [
                run.items["identity", key]["status"] for run in (control, altered)
            ]
        reason = next(
            (
                v["reason"]
                for (kind, name), v in control.items.items()
                if name == key and v["status"] == a and "reason" in v
            ),
            "control_without_decision",
        )
        log.record("mutation", key, a, control=a, mutant=b, reason=reason, **extra)
    else:
        log.check(
            key,
            a == "pass" and b == "fail",
            "mutation",
            control=a,
            mutant=b,
        )
