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
