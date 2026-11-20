# xpsloss

xpsloss is a Python library and command-line tool for researchers reproducing published models of electron energy loss in X-ray photoelectron spectroscopy (XPS). Version 0.1.0 covers oscillator energy-loss functions (ELFs), bulk differential inverse inelastic mean free paths in an infinite medium, causal dielectric response including hard gaps, and fixed-depth XPS cross sections. [日本語版](README.ja.md)

For emission at 0° from the surface normal, a positive creation depth, and a gapless material, xpsloss implements a fixed-depth effective cross section based on Eqs. (13)–(17) of Simonsen, Yubero and Tougaard (1997). The caller chooses the momentum argument of the dielectric function. `k_parallel` uses the approximation made in the closed-form derivation. At normal emission, `k_total` recovers the dispersive bulk limit. The two choices are checked differently. An independent numerical route (I3) covers `k_parallel` for dispersionless and dispersive inputs. `k_total` is checked through the bulk identity of term 13a, the deep limit, and a second implementation of the same closed forms; for dispersive inputs, no independent check covers all finite-depth terms. For dispersionless inputs the two choices coincide. Following the paper, the model assumes that the electron and the core hole are created suddenly at the same point at t = 0, that the hole is static with infinite lifetime, that the electron moves on a rectilinear trajectory at a constant velocity set by the initial kinetic energy (nonrelativistic), that losses are small compared with E0, and a specular-reflection model of the surface. Oblique emission, fixed-depth calculations for gapped materials, depth averaging, and full spectra are outside this release. [Equations and units (Japanese)](docs/EQUATIONS.md)

xpsloss is independent of QUEELS-XPS, QUASES and their authors. It is not an official, endorsed, or compatibility-guaranteed version of QUEELS-XPS and uses none of its executables, DLLs, source code, or figures.

## Quick start

From the root of the v0.1.0 source tree, use Python 3.10 or later in a virtual environment. The example material is synthetic.

```sh
python -m pip install .
xpsloss capabilities --json
xpsloss validate examples/synthetic_bulk.json
xpsloss bulk examples/synthetic_bulk.json --out runs/bulk
xpsloss dielectric examples/synthetic_bulk.json --q-inv-A 0.5 --out runs/dielectric
xpsloss xps-fixed-depth examples/synthetic_bulk.json --depth-A 10 --epsilon-momentum k_total --out runs/fixed-depth
```

Give each calculation a new output directory. The CSV contains values on the loss-energy grid. For fixed depth these include the total, the electron contribution with zero hole charge, and the signed unit-hole contribution including interference, in Å⁻¹ eV⁻¹ per path length in the solid. The accompanying `run.json` records the input path as supplied, its SHA-256, units, model, numerical settings, and runtime versions.

For `xps-fixed-depth`, `run.json` also has an `applicability` block. It records the model assumptions, the largest loss divided by E0, the checked bounds (E0 of 1000–4000 eV, max T/E0 of 0.1) with their basis, and whether the run lies inside them. At both energy bounds, all eight terms were checked against the independent route with `k_parallel` at discrete losses up to T/E0 = 0.1. Depth, material, line widths, dispersion, hole charge and the lower end of the loss grid are not part of this check. A successful run outside either bound reports a warning. Being inside the bounds establishes neither accuracy for arbitrary inputs nor physical validity.

Numerical non-convergence in the API raises `NumericalConvergenceError` or one of its subclasses. The CLI reports it as `not_converged` and invalid input as `error`. Both exit with code 2 before a new output directory is created. `dielectric` and `xps-fixed-depth` include `warning_count` in their success JSON. Their `run.json` diagnostics use `static_sign` and `kinematics`; the changelog lists the earlier names. JSON on stdout and stderr escapes non-ASCII characters, so it stays valid whatever the encoding of the console, pipe or file; `run.json` is UTF-8 and keeps them as written. [Changelog](CHANGELOG.md)

## Verification

Run `python -m unittest discover -s tests -v` from the source tree. The suite checks analytical limits, sum rules, convergence, fixed-depth terms against an independent numerical integration, and the comparison script using synthetic files. Agents independent of the implementer reviewed the equations and the implementation. There has been no human peer review. A separate reviewer agent attempted to reproduce Figs. 3, 5, and 8 with independent code using the dielectric function in Eq. (20) of Simonsen et al. (1997), outside this package's oscillator model. Of 84 comparison points, 36 lie within the figure-reading uncertainty, 18 lie outside it, and 30 depend on the assumed width (a uniform rule based on Table I of the paper; computed peaks refined beyond the grid). [Per-point results](benchmarks/paper_figures/RESULTS.md) [Validation record (Japanese)](docs/VALIDATION.md)

The 119 tests of the release candidate passed from the source tree on Windows with Python 3.12.10 (NumPy 2.5.5, SciPy 1.18.1) and Python 3.10.22 (NumPy 1.24.0, SciPy 1.10.0), and from the unpacked sdist with Python 3.12.10, without UTF-8 mode, so redirected output used the default locale encoding (cp932). The previous 112-test suite also passed from the unpacked sdist on both, and a separate reviewer ran it from the unpacked sdist on Windows with Python 3.11.14, NumPy 2.4.6 and SciPy 1.17.1. An earlier version of the suite (75 tests) passed in continuous integration on Linux with Python 3.10 and the minimum versions, on Linux and Windows with Python 3.12, and on macOS with Python 3.13, including a job that installs the wheel and tests the unpacked source distribution. The release candidate is to be checked again in continuous integration before release. [Validation record (Japanese)](docs/VALIDATION.md)

### Comparison with QUEELS-XPS v2.1

For dispersive Al from the ELF data set at 1000 eV, normal emission, and a depth of 9.967 Å, both momentum choices failed the prespecified pointwise and integrated-area criteria over 1.0–99.8 eV. The maximum absolute difference in `total`, relative to the maximum reference `total` over that interval, was 6.6% for `k_total` and 7.5% for `k_parallel`. [Original criteria and results (Japanese)](docs/M2_FIRST_CASE.md)

The follow-up hypothesis and decision rule were recorded in local Git before three new synthetic runs. Given the convergence of the reference components at low loss, the follow-up used 5.0–99.8 eV and the depth nearest 10 Å in each run. With `k_total` and a candidate expression using the undispersed oscillator energy in the numerator of the real part of the inverse dielectric function, the maximum absolute difference for each of the three components was within 0.3% of the maximum reference `total` for all three inputs. With the published expression, at least one component exceeded this threshold for each input. [Preregistered follow-up (Japanese)](docs/QUEELS_DISPERSION_TEST.md)

xpsloss uses the published Eq. (32) of Yubero et al. (1996). The comparison tested v2.1 numerical output under the stated conditions. Its internal implementation, the cause of the residual difference of order 10⁻³, and accuracy against measured spectra remain unverified. Reference-software output is excluded from distribution. [Validation record (Japanese)](docs/VALIDATION.md)

## License and citation

Keisuke Sato licenses the code under MIT and the project's own documentation under CC BY 4.0. Third-party quotations retain their owners' terms. Both grants exclude third-party data in `benchmarks/sources/**` and `benchmarks/inputs/*elfdb*`. Those files retain the original CC BY 4.0 license of N. Pauly, F. Yubero and S. Tougaard, *ELF dielectric functions for various materials*, v1.0 (2022). [MIT license](LICENSE) · [Notices and attribution](NOTICES.md) · [ELF data set](https://doi.org/10.5281/zenodo.6024064) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

Cite Keisuke Sato, *xpsloss*, version 0.1.0. The software DOI will be added at release. When using the ELF data set, cite its authors separately. [Software citation](CITATION.cff)

The Japanese documentation gives the [benchmark inputs and provenance](docs/BENCHMARK_INPUTS.md), [dielectric-response design](docs/DIELECTRIC_DESIGN.md), and [equation-to-implementation map](docs/M1_EQUATION_MAP.md).
