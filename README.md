# dcce

**Dynamic Common Correlated Effects Estimation for Panel Data**

`dcce` is an R package implementing the family of **Common Correlated Effects (CCE)**
estimators for heterogeneous coefficient panel data models with cross-sectional
dependence. It is an R port of Jan Ditzen's `xtdcce2` Stata package and provides
the standard estimators of Pesaran (2006), Chudik & Pesaran (2015), and related
contributions, together with a comprehensive cross-sectional dependence (CD) test
suite.

---

## Features

### Estimators

| Estimator | Reference | Notes |
|-----------|-----------|-------|
| Mean Group (MG) | Pesaran & Smith (1995) | heterogeneous slopes |
| Common Correlated Effects (CCE-MG) | Pesaran (2006) | static, with CSAs |
| Dynamic CCE (DCCE) | Chudik & Pesaran (2015) | dynamic panel + CSA lags |
| Regularized CCE (rCCE) | Juodis (2022) | PCA-regularized CSA factor |
| CS-DL (long-run) | Chudik et al. (2016) | direct LR via level of x |
| CS-ARDL (short + long run) | Chudik et al. (2016) | full SR / adjustment / LR blocks via delta method |
| Pooled Mean Group (PMG) | Shin, Pesaran & Smith (1999) | inverse-variance pooled LR |

All three long-run estimators produce a three-block output: **short-run**
coefficients, the **adjustment** (speed of return to equilibrium), and
**long-run** elasticities with delta-method standard errors.

### Cross-sectional dependence tests

| Test | Reference | Description |
|------|-----------|-------------|
| CD | Pesaran (2015) | benchmark Pesaran CD |
| CDw | Juodis & Reese (2022) | Rademacher-weighted |
| **CDw+** | Baltagi, Feng & Kao (2012) | bias-adjusted LM with weighting |
| PEA | Fan, Liao & Yao (2015) | power-enhanced for sparse alternatives |
| CD\* | Pesaran & Xie (2021) | bias-corrected for strong factors |

### Other diagnostics

| Tool | Reference |
|------|-----------|
| **Pesaran CIPS panel unit root test** | Pesaran (2007) |
| **Swamy / Pesaran-Yamagata slope heterogeneity test** | Swamy (1970); Pesaran & Yamagata (2008) |
| **Hausman-style MG vs Pooled test** | — |
| Exponent of cross-sectional dependence | Bailey, Kapetanios & Pesaran (2016, 2019) |
| IC for CSA selection | Margaritella & Westerlund (2023) |
| Rank condition classifier | De Vos, Everaert & Sarafidis (2024) |
| Cross-section / wild bootstrap inference | — |

### Extensions

| Tool | Description |
|------|-------------|
| **`dcce_rolling()`** | Rolling-window estimation with coefficient path tibble and `plot` method |
| **`absorb` argument** | High-dimensional fixed-effect absorption via alternating projections |
| **`spatial_weights` argument** | Spatial CCE with user-supplied weight matrix |
| **`structural_break_test()`** | Chow / sup-Wald tests, breakdate estimation, sequential Bai-Perron (R port of Stata `xtbreak`) |

### S3 methods and ergonomics

- `broom`-compatible `tidy()` and `glance()` (tidy includes short-run, adjustment, and long-run rows for LR estimators)
- `confint()` with `type = c("mg", "lr", "adjustment")`
- `plot()` for unit-level coefficient histograms and residual diagnostics
- `update()` for refitting with modified arguments
- `coef(fit, type = "unit")` for unit-level coefficient extraction
- Native support for `L()`, `D()`, and `Lrange()` operators in formulas (xtdcce2-compatible syntax)
- Unbalanced panel handling

---

## Installation

The package is in active development and not yet on CRAN. Install the development
version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Mustapha-Wasseja/dcce")
```

To build the vignette during installation:

```r
remotes::install_github("Mustapha-Wasseja/dcce", build_vignettes = TRUE)
```

### System requirements

- **R** >= 4.1.0
- **Imports:** `stats`, `Matrix`, `collapse` (>= 2.0.0), `sandwich`, `generics`,
  `rlang` (>= 1.1.0), `cli` (>= 3.0.0), `tibble`
- **Suggests:** `broom`, `ggplot2`, `lifecycle`, `plm`, `testthat` (>= 3.0.0),
  `knitr`, `rmarkdown`

---

## Quick start

```r
library(dcce)

# Load the bundled Penn World Tables 8 dataset (93 countries, 1960-2007)
data(pwt8)

# Fit a Dynamic CCE growth regression with 3 lags of CSAs
fit <- dcce(
  data               = pwt8,
  unit_index         = "country",
  time_index         = "year",
  formula            = d_log_rgdpo ~ L(log_rgdpo, 1) + log_hc + log_ck + log_ngd,
  model              = "dcce",
  cross_section_vars = ~ log_rgdpo + log_hc + log_ck + log_ngd,
  cross_section_lags = 3
)

print(fit)

# Verify that DCCE has removed cross-sectional dependence
pcd_test(fit, test = "pesaran")

# Tidy output (broom compatible)
tidy(fit)
glance(fit)
```

For a complete walkthrough including motivation, theory, all estimators, and the
Ditzen (2018) replication, see the package vignette:

```r
vignette("dcce-introduction", package = "dcce")
```

---

## Verify your installation

The fastest way to confirm the package works on your system is to run a few
of the worked examples on the bundled datasets:

```r
library(dcce)

# Example 1: Mean Group on the simulated dataset
data(dcce_sim)
fit_mg <- dcce(
  data = dcce_sim, unit_index = "unit", time_index = "time",
  formula = y ~ L(y, 1) + x,
  model = "mg", cross_section_vars = NULL
)
coef(fit_mg)

# Example 2: DCCE with CD test on residuals
data(pwt8)
fit_dcce <- dcce(
  data = pwt8, unit_index = "country", time_index = "year",
  formula = d_log_rgdpo ~ L(log_rgdpo, 1) + log_hc + log_ck + log_ngd,
  model = "dcce",
  cross_section_vars = ~ log_rgdpo + log_hc + log_ck + log_ngd,
  cross_section_lags = 3
)
pcd_test(fit_dcce, test = "pesaran")  # Should be insignificant after DCCE

# Example 3: Bootstrap inference
set.seed(42)
boot <- bootstrap(fit_dcce, type = "crosssection", reps = 199)
print(boot)
```

### Optional: validate against `plm`

If you have `plm` installed, the package's static CCE estimator matches
`plm::pmg(..., model = "cmg")` to three decimal places on the Produc
dataset. This is checked automatically by the bundled
`tests/testthat/test-produc-validation.R` file.

---

## Usage notes

**Formula operators.** The package extends standard R formulas with three
panel-aware operators:

- `L(x, k)` — k-th lag of `x` within each unit
- `D(x, k)` — k-th difference of `x` within each unit
- `Lrange(x, k0, k1)` — lags `k0` through `k1` (used in CS-ARDL)

**Cross-section variables.** Use `cross_section_vars = ~ .` to include all
regressors plus the dependent variable as CSAs (the default), or provide an
explicit one-sided formula such as `~ log_rgdpo + log_hc`.

**CSA lags.** For dynamic models the Chudik-Pesaran rule
`p_T = floor(T^(1/3))` is the standard recommendation
(`cross_section_lags = 3` for `T ≈ 30-50`).

---

## References

- Bailey, N., Kapetanios, G., & Pesaran, M. H. (2016). Exponent of cross-sectional dependence: estimation and inference. *Journal of Applied Econometrics*, 31(6), 929–960.
- Chudik, A., & Pesaran, M. H. (2015). Common correlated effects estimation of heterogeneous dynamic panel data models with weakly exogenous regressors. *Journal of Econometrics*, 188(2), 393–420.
- Chudik, A., Mohaddes, K., Pesaran, M. H., & Raissi, M. (2016). Long-run effects in large heterogeneous panel data models with cross-sectionally correlated errors. In *Essays in Honor of Aman Ullah*, 36, 85–135.
- De Vos, I., Everaert, G., & Sarafidis, V. (2024). Rank condition for the CCE estimator with heterogeneous slopes. *Journal of Econometrics*, 240(2), 105703.
- Ditzen, J. (2018). Estimating dynamic common-correlated effects in Stata. *The Stata Journal*, 18(3), 585–617.
- Fan, J., Liao, Y., & Yao, J. (2015). Power enhancement in high-dimensional cross-sectional tests. *Econometrica*, 83(4), 1497–1541.
- Juodis, A., & Reese, S. (2022). A randomized CD test for error cross-sectional dependence in panel data models with possibly many regressors. *Journal of Business & Economic Statistics*, 40(4), 1500–1514.
- Margaritella, L., & Westerlund, J. (2023). Information criteria for factor model selection in panel data. *Econometric Reviews*, 42(7), 619–641.
- Pesaran, M. H. (2006). Estimation and inference in large heterogeneous panels with a multifactor error structure. *Econometrica*, 74(4), 967–1012.
- Pesaran, M. H. (2015). Testing weak cross-sectional dependence in large panels. *Econometric Reviews*, 34(6–10), 1089–1117.
- Pesaran, M. H., & Smith, R. (1995). Estimating long-run relationships from dynamic heterogeneous panels. *Journal of Econometrics*, 68(1), 79–113.
- Shin, Y., Pesaran, M. H., & Smith, R. (1999). An autoregressive distributed-lag modelling approach to cointegration analysis. In *Econometrics and Economic Theory in the 20th Century*, 371–413.

## License

GPL (>= 3)

## Citation

If you use `dcce` in published work, please cite the package and the relevant
methodological references above.

```r
citation("dcce")
```

## Issues

Bug reports and feature requests are welcome at
<https://github.com/Mustapha-Wasseja/dcce/issues>.
