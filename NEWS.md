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
