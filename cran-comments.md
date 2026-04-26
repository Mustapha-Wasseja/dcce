## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.5.2, Rtools 4.5
* GitHub Actions: ubuntu-latest, windows-latest, macos-latest (R release)

## Notes

* This is a new submission.
* The package uses RcppArmadillo for an optional C++ inner loop that
  accelerates unit-level OLS. A pure-R fallback is always available
  and is used automatically if the compiled routines are not loaded.
* `marginaleffects` is in Suggests; S3 methods are registered
  dynamically in `.onLoad()` only when that package is installed.
* The exported function `D()` masks `stats::D()` (symbolic
  differentiation). This is intentional — `dcce::D()` is a panel-data
  differencing operator used inside `dcce()` formulas. Users who need
  symbolic differentiation can call `stats::D()` explicitly.
