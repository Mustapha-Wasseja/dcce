# dcce 0.4.0

## New estimators

* **Interactive Fixed Effects (IFE)** estimator of Bai (2009), via
  `model = "ife"` in `dcce()`. Iterates between pooled OLS and
  principal-component extraction of factors/loadings until convergence.
  The number of factors can be specified or selected automatically via
  the BIC3 criterion of Bai & Ng (2002). Unlike CCE, IFE estimates
  factors and loadings directly rather than proxying them with
  cross-sectional averages.

* **Pooled CCE (CCEP)** estimator, via `model = "ccep"` in `dcce()`.
  Constrains slopes to be identical across units (precision-weighted
  pooled OLS with cross-sectional averages), as proposed by
  Pesaran (2006) alongside the existing CCE-MG. More efficient than
  CCE-MG when slopes are truly homogeneous.

## New diagnostics

* **Dumitrescu-Hurlin (2012) panel Granger causality test** via the
  new exported function `granger_test()`. Reports the W-bar statistic
  (cross-sectional average of unit Wald statistics), the Z-bar
  (large-T standardised), and the Z-bar tilde (small-sample adjusted)
  statistics.

* **IPS and LLC panel unit root tests** via `panel_ur_test()`. Implements
  Im, Pesaran & Shin (2003) IPS t-bar and Levin, Lin & Chu (2002) LLC
  common-root tests. These do not correct for CSD (use `cips_test()`
  for that); they are included as standard benchmarks.

* **Pedroni and Kao panel cointegration tests** via `panel_coint_test()`.
  Reports the Pedroni (1999, 2004) group-mean t and group-mean rho
  statistics, or the Kao (1999) ADF statistic.

## New features

* **Impulse response functions** via the new exported function `irf()`.
  Computes IRFs from fitted dynamic panel models (DCCE, CS-ARDL, PMG)
  using the MG ARDL lag polynomial. Optional cross-section bootstrap
  for confidence bands. Includes `print.dcce_irf()` and
  `plot.dcce_irf()` methods.

* **Half-panel jackknife** bias correction (Chudik & Pesaran 2015) via
  `bias_correction = "half_panel_jackknife"` in `dcce()`. Splits each
  unit's time series in half, fits on each half, and corrects the
  full-sample MG estimate: b_hpj = 2*b_full - 0.5*(b_half1 + b_half2).
  Targets the Nickell bias in dynamic CCE.

## Infrastructure

* **GitHub Actions CI** added (`.github/workflows/R-CMD-check.yaml`):
  runs `R CMD check --as-cran` on Ubuntu, Windows, and macOS on every
  push/PR to `main`.
* **pkgdown site** configuration added (`_pkgdown.yml`) with grouped
  reference index.
* **`cran-comments.md`** added for CRAN submission readiness.

## References

* Bai, J. (2009). Panel data models with interactive fixed effects.
  *Econometrica*, 77(4), 1229-1279.
* Dumitrescu, E.-I., & Hurlin, C. (2012). Testing for Granger
  non-causality in heterogeneous panels. *Economic Modelling*, 29(4),
  1450-1460.
* Im, K. S., Pesaran, M. H., & Shin, Y. (2003). Testing for unit
  roots in heterogeneous panels. *Journal of Econometrics*, 115(1),
  53-74.
* Kao, C. (1999). Spurious regression and residual-based tests for
  cointegration in panel data. *Journal of Econometrics*, 90(1), 1-44.
* Levin, A., Lin, C.-F., & Chu, C.-S. J. (2002). Unit root tests in
  panel data. *Journal of Econometrics*, 108(1), 1-24.
* Pedroni, P. (2004). Panel cointegration. *Econometric Theory*,
  20(3), 597-625.


# dcce 0.3.2

## Bug fixes

* **`structural_break_test(type = "unknown")` now uses the correct
  asymptotic distribution for the sup-Wald statistic.** The previous
  implementation applied a Bonferroni correction over the candidate
  breakdate grid, which was wildly over-conservative and routinely
  inflated borderline p-values all the way to 1. The sup-Wald statistic
  has a non-standard asymptotic distribution (supremum of a squared
  Brownian bridge) whose critical values were tabulated by Andrews
  (1993, *Econometrica* 61(4), Table I). `dcce` now ships those
  critical values as internal data and reports p-values by interpolating
  on the log-scale between the tabulated 1%, 5%, and 10% levels.
  Values above the 1% critical value are reported as `p <= 0.01`; below
  the 10% critical value they are reported as `p > 0.10`.
* `print.dcce_break()` for unknown-date tests now displays the full
  set of Andrews critical values next to the sup-Wald statistic, and
  formats the p-value as a human-readable bracket ("> 0.10", "<= 0.01",
  or the interpolated value) instead of a bare number. The known-date
  (Chow) path is unchanged and still uses the standard chi-square
  p-value.
* The `dcce_break` object now carries a `critical_values` element (a
  named vector with `cv10`, `cv05`, `cv01`) and a `trim` element so
  downstream code and custom print methods can access them.
* **Silenced Armadillo runtime warnings.** `src/unit_ols.cpp` now
  defines `ARMA_WARN_LEVEL 0` in addition to `ARMA_DONT_PRINT_ERRORS`,
  fully suppressing the `solve(): system is singular` messages that
  used to leak through to the console during `structural_break_test()`
  and near-rank-deficient CCE fits. The R-level `.unit_ols()` fallback
  still handles rank deficiency via `pinv()` and continues to log its
  own `cli_warn` diagnostics where appropriate.

## Internal

* New tests in `test-structural-break.R` covering the Andrews
  critical-value lookup, p-value interpolation, and a regression
  test reproducing the v0.3.1 Bonferroni blow-up (sup-Wald = 7.20,
  q = 3, pi0 = 0.15) to lock in the fix.


# dcce 0.3.1

## Bug fixes / ergonomics

* **New `dcce_bootstrap()` alias.** `broom` exports a function called
  `bootstrap()` (for resampling data frames) that has a completely
  different signature from `dcce::bootstrap()`. If `broom` is loaded
  after `dcce`, the broom version masks ours on the search path and
  calls of the form `bootstrap(fit, type = "crosssection", reps = 199)`
  fail with an "unused arguments" error. This release adds
  `dcce_bootstrap()`, an exported alias with identical semantics to
  `dcce::bootstrap()` but a conflict-free name. `bootstrap()` itself is
  unchanged and still available; users who previously wrote
  `dcce::bootstrap(fit, ...)` can keep doing so.
* `?bootstrap` now documents the conflict and three workarounds.
* README updated with the conflict note under "Verify your
  installation".


# dcce 0.3.0

## New features

* **`dcce_rolling()` for rolling-window estimation.** Fits a sequence
  of `dcce()` models on overlapping time windows and returns a
  `dcce_rolling` object containing the full list of fits and a tidy
  tibble of coefficients indexed by window end-date. Includes
  `print.dcce_rolling()` and a `plot.dcce_rolling()` method that draws
  one coefficient-path panel per regressor with a 95% confidence
  ribbon. Useful for detecting parameter drift and regime shifts in
  long panels.

* **`absorb` argument for high-dimensional fixed effects.** A new
  argument to `dcce()` that accepts a one-sided formula (or character
  vector) naming grouping factors to project out of \eqn{y} and
  \eqn{X} before the main unit loop runs. A single factor uses the
  within transformation; multiple factors use the alternating
  projections of Guimaraes & Portugal (2010) / Correia (2016). The
  unit fixed effects of CCE estimators are still kept via unit
  intercepts — `absorb` is for *additional* categorical effects on top
  of the cross-section (e.g. industry, region, sub-period dummies).

* **Spatial CCE via the new `spatial_weights` argument.** When supplied
  with an \eqn{N \times N} row-normalised weight matrix, `dcce()`
  replaces the global cross-sectional averages of classical CCE with
  **local**, unit-specific weighted averages
  \eqn{\bar y^W_{i,t} = \sum_j w_{ij} y_{j,t}}. The matrix is
  automatically row-normalised, its diagonal is zeroed, and rows/columns
  are aligned to the panel's unit identifiers if they are named. This
  enables spatial CCE estimation that respects the topology of
  cross-sectional dependence (geographical contiguity, trade links,
  input-output connections).

* **`structural_break_test()` for panel structural breaks.** New
  exported function implementing a Wald-type Chow test at a known
  break date, a sup-Wald test over a trimmed candidate window for an
  unknown break date, and a sequential Bai & Perron (1998) procedure
  for multiple breaks. Returns a `dcce_break` object with the test
  statistic, p-value, estimated break date(s), the full Wald profile
  across candidates, and fitted `dcce_fit` objects for the pre- and
  post-break regimes. Works with any base estimator supported by
  `dcce()`. References: Andrews (1993); Bai & Perron (1998); Ditzen,
  Karavias & Westerlund (2024). This closes the gap with Stata's
  `xtbreak` package.

## Internal

* Added `R/rolling.R`, `R/absorb.R`, `R/spatial_cce.R`, and
  `R/structural_break.R`.
* Added `ARMA_DONT_PRINT_ERRORS` to `src/unit_ols.cpp` to silence
  Armadillo diagnostics on near-singular `XtX` matrices. The R-level
  `.unit_ols()` fallback still handles rank deficiency via `pinv()`.
* 60+ new test assertions across four new test files.

## References

* Andrews, D. W. K. (1993). Tests for parameter instability and
  structural change with unknown change point. *Econometrica*,
  61(4), 821-856.
* Bai, J., & Perron, P. (1998). Estimating and testing linear models
  with multiple structural changes. *Econometrica*, 66(1), 47-78.
* Correia, S. (2016). Linear models with high-dimensional fixed
  effects: An efficient and feasible estimator. Working paper.
* Ditzen, J., Karavias, Y., & Westerlund, J. (2024). Multiple
  structural breaks in interactive effects panel data models.
  *Journal of Applied Econometrics*.
* Guimaraes, P., & Portugal, P. (2010). A simple feasible procedure
  to fit models with high-dimensional fixed effects. *Stata Journal*,
  10(4), 628-649.


# dcce 0.2.0

## New features

* **RcppArmadillo C++ inner OLS loop** for the unit-level estimator. The
  pure-R `.unit_ols()` is retained as an automatic fallback, and the new
  `fast = TRUE` argument to `dcce()` controls which path is used. Pure-R
  and C++ paths are numerically identical to well under the `1e-6`
  tolerance on all bundled tests. Compilation requires the standard R
  toolchain (Rtools on Windows, Xcode command-line tools on macOS).

* **Parallel unit-level estimation** via the new `n_cores` argument to
  `dcce()`. When `n_cores > 1L` on Unix/macOS, the unit loop is
  dispatched to `parallel::mclapply()`; on Windows the argument is
  silently ignored.

* **Augmented Mean Group (AMG) estimator** (Eberhardt & Teal 2010; Bond
  & Eberhardt 2013), available as `model = "amg"` in `dcce()`. The
  estimator extracts a Common Dynamic Process (CDP) from a pooled
  first-difference regression with time dummies, cumulates it within
  each unit, and adds the level-form CDP as a nuisance regressor in the
  unit-level OLS. The CDP is stored on the fit object as `fit$cdp`.

* **Westerlund (2007) panel cointegration tests** exposed via the new
  exported function `cointegration_test()`. Supports all four
  statistics from Westerlund (2007): Gt, Ga, Pt, Pa. Asymptotic
  p-values use the Westerlund (2007) Table 3 critical values; a
  cross-sectional bootstrap path (`n_bootstrap > 0`) is also provided
  for dependence-robust inference. Returns a `dcce_cointegration`
  object with a dedicated print method.

* **`dcce_workflow()` diagnostic pipeline** that runs the recommended
  pre-estimation sequence and returns a structured report plus a
  printable recommended `dcce()` call. Steps include the panel summary,
  CIPS unit root tests on each variable, a pooled CD test,
  (conditionally) the Westerlund cointegration test, the rank condition
  classifier, and IC-based CSA lag selection. Returns a
  `dcce_workflow` object.

* **`marginaleffects` compatibility.** S3 methods on the internal
  `marginaleffects` generics (`get_coef`, `get_vcov`, `get_predict`,
  `find_predictors`, `find_response`) are registered dynamically in
  `.onLoad()` when the `marginaleffects` package is available, enabling
  `avg_slopes()`, `avg_predictions()`, and `hypotheses()` on `dcce_fit`
  objects. `marginaleffects` is in `Suggests` only, so the integration
  is a zero-cost soft dependency.

* **`predict.dcce_fit()` now supports `newdata`.** When a new data
  frame is supplied, predictions are computed using the Mean Group
  coefficients on the structural regressors.

## Internal

* Added `src/unit_ols.cpp` with RcppArmadillo implementations of
  `unit_ols_cpp()` and `batch_ols_cpp()`.
* Added `src/Makevars` and `src/Makevars.win` linking to R's
  `LAPACK_LIBS`, `BLAS_LIBS`, and `FLIBS`.
* Added `Rcpp (>= 1.0.0)` to `Imports`; added `Rcpp` and
  `RcppArmadillo` to `LinkingTo`.
* Added `marginaleffects` and `parallel` to `Suggests`.
* Added `.run_unit_loop()` dispatcher in `R/est_mg.R`.
* Added `.onLoad()` hook in `R/zzz.R` for dynamic method registration.

## References

* Bond, S., & Eberhardt, M. (2013). Accounting for unobserved
  heterogeneity in panel time series models. *Economics Letters*.
* Eberhardt, M., & Teal, F. (2010). Productivity Analysis in Global
  Manufacturing Production. Economics Series Working Paper 515, Oxford.
* Westerlund, J. (2007). Testing for Error Correction in Panel Data.
  *Oxford Bulletin of Economics and Statistics*, 69(6), 709-748.


# dcce 0.1.0

* Initial release.
* Mean Group (MG), Common Correlated Effects (CCE), and Dynamic CCE (DCCE)
  estimators.
* Regularized CCE (rCCE) with Ahn-Horenstein ER/GR criteria.
* Long-run estimators: PMG, CS-DL, CS-ARDL.
* Cross-sectional dependence tests: CD, CDw, PEA, CD*.
* Exponent of cross-sectional dependence (BKP 2016, 2019).
* Information criteria for CSA selection.
* Rank condition classifier.
* Cross-section and wild bootstrap inference.
