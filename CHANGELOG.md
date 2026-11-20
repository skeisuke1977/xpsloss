# Changelog

## [0.1.0] - unreleased

### Added

- M0 calculates oscillator energy-loss functions and bulk differential inverse inelastic mean free paths and records the input and numerical settings for each run.
- M1 provides causal dielectric response, Kramers–Kronig reconstruction for hard gaps, and traceable Al, Au, and Si inputs from the ELF data set.
- M2-1 calculates the fixed-depth XPS effective cross section for normal emission from gapless materials, including surface, vacuum, core-hole, and interference terms.

### Verification

- Tests cover analytical limits, sum rules, convergence, and an independent numerical integration of the fixed-depth terms. Agents independent of the implementer reviewed the equations and implementation. There has been no human peer review. A separate reviewer agent used independent code with Eq. (20) of Simonsen et al. (1997) to reproduce Figs. 3, 5, and 8 within figure-reading uncertainty. [Validation record](docs/VALIDATION.md)
- Both momentum choices failed the preregistered QUEELS-XPS v2.1 comparison for dispersive Al. In three preregistered synthetic follow-up runs, all components of a candidate expression agreed with the v2.1 numerical output within the preset tolerance. The published expression used by xpsloss exceeded the tolerance for at least one component in each run. The comparison did not examine QUEELS-XPS internals. [Original comparison](docs/M2_FIRST_CASE.md) · [Follow-up](docs/QUEELS_DISPERSION_TEST.md)

### Notes

- The QUEELS-XPS v2.1 User's Guide is cited for reference. Its Si input record and tutorial values are excluded from the release.
