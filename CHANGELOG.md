# Changelog

## [0.1.0] - unreleased

### Added

- M0 calculates oscillator energy-loss functions and bulk differential inverse inelastic mean free paths and records the input and numerical settings for each run.
- M1 provides causal dielectric response, Kramers–Kronig reconstruction for hard gaps, and traceable Al, Au, and Si inputs from the ELF data set.
- M2-1 calculates the fixed-depth XPS effective cross section for normal emission from gapless materials, including surface, vacuum, core-hole, and interference terms.

### Verification

- Tests cover analytical limits, sum rules, convergence, and an independent numerical integration of the fixed-depth terms. Agents independent of the implementer reviewed the equations and implementation. There has been no human peer review. A separate reviewer agent attempted to reproduce Figs. 3, 5, and 8 with independent code using Eq. (20) of Simonsen et al. (1997), outside the package's oscillator model. Of 84 comparison points, 36 lie within the figure-reading uncertainty, 18 outside it, and 30 depend on the assumed width. [Per-point results](benchmarks/paper_figures/RESULTS.md) · [Validation record](docs/VALIDATION.md)
- The two momentum choices of the fixed-depth model are checked differently: the independent route I3 covers `k_parallel` (dispersionless and dispersive); `k_total` is checked through the 13a bulk identity, the deep limit and a second implementation of the same closed forms, with no independent check of all finite-depth terms for dispersive inputs.
- Both momentum choices failed the preregistered QUEELS-XPS v2.1 comparison for dispersive Al. In three preregistered synthetic follow-up runs, all components of a candidate expression agreed with the v2.1 numerical output within the preset tolerance. The published expression used by xpsloss exceeded the tolerance for at least one component in each run. The comparison did not examine QUEELS-XPS internals. [Original comparison](docs/M2_FIRST_CASE.md) · [Follow-up](docs/QUEELS_DISPERSION_TEST.md)
- The 119 tests of the release candidate passed from the source tree on Windows with Python 3.10.22 (NumPy 1.24.0, SciPy 1.10.0) and 3.12.10 (NumPy 2.5.5, SciPy 1.18.1), and from the unpacked sdist with 3.12.10, with the default locale encoding (cp932) for redirected output. The previous 112-test suite also passed from the unpacked sdist on both, and a separate reviewer passed it from the sdist with Python 3.11.14 (NumPy 2.4.6, SciPy 1.17.1) on Windows. The earlier 75-test suite passed the four CI environments and the packaging job at `10a8dfd`. At commit `c7fdd22` the release candidate passed all five CI jobs (Linux with Python 3.10 and the minimum dependencies, Linux 3.12, macOS 3.13, Windows 3.12, and the packaging job with the unpacked sdist), including both benchmark checks. [Validation record](docs/VALIDATION.md)

### Changed (pre-release, phase 4e)

- The CLI writes JSON on stdout and stderr with non-ASCII characters escaped, so the output stays valid JSON whatever the encoding of the console, pipe or file. Before, on Windows without UTF-8 mode, output sent to a pipe or file used the locale encoding (e.g. cp932). A character outside it made `validate` (material name) and every calculation command (absolute output path, after the output directory and `run.json` had been written) end with `error`, and made error JSON on stderr containing such a character unreadable. `run.json` stays UTF-8 and keeps the characters.
- Input strings that cannot be encoded as UTF-8 (e.g. a lone surrogate in a material name) are rejected with `ValueError` (CLI `error`) when the input is read, before any output directory is created.
- `evaluate_prereg.py predict --check` also checks the bundled predictions against the registered SHA-256. `generate.py --check` compares the regenerated tables exactly for text, verdicts and counts, and with stated tolerances for numbers (absolute 2×10⁻⁶ eV for peak positions; relative 10⁻⁹ for the other CSV values; see the validation record).
- CI runs both checks in all four test environments and from the unpacked sdist.

### Added (pre-release, phase 4c)

- `benchmarks/queels_prereg/`: the derived candidate code of the preregistered QUEELS-XPS v2.1 follow-up (without the reference-software value table), the frozen predictions (byte-identical to the registered file) and a portable evaluator that applies the registered decision rule to the reader's own exports; `predict --check` verifies the predictions byte by byte. [README](benchmarks/queels_prereg/README.md)
- `benchmarks/paper_figures/`: independent Eq. (20) code for Figs. 3, 5 and 8 of Simonsen et al. (1997), cited figure readings, a uniform width rule based on Table I, refined computed peaks, and a generated per-point comparison; `generate.py --check` regenerates and compares the tables (with tolerances for numbers since phase 4e). [Results](benchmarks/paper_figures/RESULTS.md)
- The source distribution ships these folders and `benchmarks/numerics/`.

### Changed (pre-release numerical robustness, phase 4a)

- Phase 4a corrected the narrow-resonance failures found in the bulk, dielectric, optical sum-rule and fixed-depth calculations. In the recorded regression tests and scans, successful results met the comparison tolerances. Unresolved numerical cases raised convergence errors. An insufficient initial panel limit remains an input error. The physical model (equations, defaults, units and constants) is unchanged. [Validation record](docs/VALIDATION.md)
- New exception base class `NumericalConvergenceError` (a `RuntimeError`); `DielectricConvergenceError` and `FixedDepthConvergenceError` now derive from it, and `BulkConvergenceError` is new. Bulk DIIMFP raises `BulkConvergenceError` instead of `ValueError` when it cannot resolve an input, so `except ValueError` no longer catches numerical non-convergence there. Arithmetic overflow of the oscillator ELF for finite inputs raises `NumericalConvergenceError`.
- Integration settings: a non-positive tolerance, an invalid `limit`, or a `limit` too small for the initial panels raises `ValueError` (CLI `error`). The former special case for `limit == 1` is removed. Failure to reach the tolerance after the initial panels raises the path's convergence error (CLI `not_converged`). Both exit with code 2.
- Input validation and numerical failures occur before the CLI creates a new output directory. It rejects an existing output path before calculation with `FileExistsError` (`error`) and uses exclusive directory creation. An I/O failure during output can leave a partial directory. Numerical errors identify the affected loss and its array index when available.

### Changed (pre-release, phase 4b)

- `xps-fixed-depth` writes an `applicability` block to `run.json`: the largest loss divided by E0, the electron energy, the momentum choice, the assumptions of Simonsen et al. (1997) (sudden pair creation at t = 0, static core hole with infinite lifetime, rectilinear trajectory at a constant nonrelativistic velocity from the initial kinetic energy, T ≪ E0, specular-reflection surface model), the numerical bounds (E0 1000–4000 eV, max T/E0 ≤ 0.1) with the basis of each bound, and whether the run lies inside them. Only these two quantities are compared. At both energy endpoints, all eight terms were checked against the independent route with `k_parallel` at discrete losses up to T/E0 = 0.1 (`tests/test_applicability_envelope.py`). For a successful run outside either bound, `warnings` says that the numerical checks do not cover the input and that physical applicability is not verified. `diagnostics.velocity` of the fixed-depth `run.json` is renamed `diagnostics.kinematics`. The success message on stdout of `xps-fixed-depth` and `dielectric` gives `warning_count`.
- `DielectricDiagnostics.static_status` (`stable` / `marginal` / `unstable`) is renamed `static_sign` (`negative` / `marginal` / `positive`), also in the `dielectric` `run.json` diagnostics. It is the sign of Re[−1/ε(q, 0)] of the oscillator representation relative to the tolerance, not a test of material stability. The rejection rules are unchanged. `schema_version` stays 1 (unreleased).
- `capabilities` states the checks per momentum choice and no longer claims that the figure reproduction is within reading uncertainty.

### Changed (pre-release, phase 4d)

- Documentation updates. In user-visible strings, `capabilities` now states the fixed-depth scope as a positive depth, 0° from the surface normal and a gapless material, and describes how I5 was checked. The help texts and the error message for an unsupported input `model` were reworded. Calculations and JSON keys are unchanged.

### Notes

- The QUEELS-XPS v2.1 User's Guide is cited for reference. Its Si input record and tutorial values are excluded from the release.
