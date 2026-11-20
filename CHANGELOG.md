# Changelog

## [0.1.0] - unreleased

### Added

- M0 calculates oscillator energy-loss functions and bulk differential inverse inelastic mean free paths and records the input and numerical settings for each run.
- M1 provides causal dielectric response, Kramers–Kronig reconstruction for hard gaps, and traceable Al, Au, and Si inputs from the ELF data set.
- M2-1 calculates the fixed-depth XPS effective cross section for normal emission from gapless materials, including surface, vacuum, core-hole, and interference terms.

### Verification

- Tests cover analytical limits, sum rules, convergence, and an independent numerical integration of the fixed-depth terms. Agents independent of the implementer reviewed the equations and implementation. There has been no human peer review. A separate reviewer agent attempted to reproduce Figs. 3, 5, and 8 with independent code using Eq. (20) of Simonsen et al. (1997); per-point results are to be listed in the validation record before release (phase 4c). [Validation record](docs/VALIDATION.md)
- The two momentum choices of the fixed-depth model are checked differently: the independent route I3 covers `k_parallel` (dispersionless and dispersive); `k_total` is checked through the 13a bulk identity, the deep limit and a second implementation of the same closed forms, with no independent check of all finite-depth terms for dispersive inputs.
- Both momentum choices failed the preregistered QUEELS-XPS v2.1 comparison for dispersive Al. In three preregistered synthetic follow-up runs, all components of a candidate expression agreed with the v2.1 numerical output within the preset tolerance. The published expression used by xpsloss exceeded the tolerance for at least one component in each run. The comparison did not examine QUEELS-XPS internals. [Original comparison](docs/M2_FIRST_CASE.md) · [Follow-up](docs/QUEELS_DISPERSION_TEST.md)

### Changed (pre-release numerical robustness, phase 4a)

- Integrals over very narrow resonances no longer return wrong values as success. Bulk DIIMFP, the dielectric response (gapless oscillators, hard-gap and tabulated Kramers–Kronig, the gapped closed form), the optical sum rules and the fixed-depth terms either resolve such inputs within the stated tolerances or stop with a convergence error. The physical model is unchanged. [Validation record](docs/VALIDATION.md)
- New exception base class `NumericalConvergenceError` (a `RuntimeError`); `DielectricConvergenceError` and `FixedDepthConvergenceError` now derive from it, and `BulkConvergenceError` is new. Bulk DIIMFP raises `BulkConvergenceError` instead of `ValueError` when it cannot resolve an input, so `except ValueError` no longer catches numerical non-convergence there. Arithmetic overflow of the oscillator ELF for finite inputs raises `NumericalConvergenceError`.
- Integration settings: a non-positive tolerance, an invalid `limit`, or a `limit` too small for the initial panels raises `ValueError` (CLI `error`). The former special case for `limit == 1` is removed. Failure to reach the tolerance after the initial panels raises the path's convergence error (CLI `not_converged`). Both exit with code 2.
- CLI: on failure no output directory and no `run.json` are written. An existing output path is rejected before the computation with `FileExistsError` (`error`), and the final directory creation remains exclusive. Error messages give the offending loss value and index instead of the whole grid.

### Changed (pre-release, phase 4b)

- `xps-fixed-depth` writes an `applicability` block to `run.json`: the largest loss divided by E0, the electron energy, the momentum choice, the assumptions of Simonsen et al. (1997) (sudden pair creation at t = 0, static core hole with infinite lifetime, rectilinear trajectory at a constant nonrelativistic velocity from the initial kinetic energy, T ≪ E0, specular-reflection surface model), the ranges at whose ends all eight terms were checked against the independent route with `k_parallel` (max T/E0 ≤ 0.1, E0 1000–4000 eV; new test `tests/test_applicability_envelope.py`) with the basis of each bound, and whether the run lies inside them. Only these two quantities are compared. Outside the ranges the run completes and `warnings` says that the numerical checks do not cover the input and that physical applicability is not verified. `diagnostics.velocity` of the fixed-depth `run.json` is renamed `diagnostics.kinematics`. The success message on stdout of `xps-fixed-depth` and `dielectric` gives `warning_count`.
- `DielectricDiagnostics.static_status` (`stable` / `marginal` / `unstable`) is renamed `static_sign` (`negative` / `marginal` / `positive`), also in the `dielectric` `run.json` diagnostics. It is the sign of Re[−1/ε(q, 0)] of the oscillator representation relative to the tolerance, not a test of material stability. The rejection rules are unchanged. `schema_version` stays 1 (unreleased).
- `capabilities` states the checks per momentum choice and no longer claims that the figure reproduction is within reading uncertainty.

### Notes

- The QUEELS-XPS v2.1 User's Guide is cited for reference. Its Si input record and tutorial values are excluded from the release.
