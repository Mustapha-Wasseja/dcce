test_that("Rcpp and pure-R paths produce numerically identical coefficients", {
  data(dcce_sim)

  fit_r <- dcce(
    data               = dcce_sim,
    unit_index         = "unit",
    time_index         = "time",
    formula            = y ~ L(y, 1) + x,
    model              = "dcce",
    cross_section_vars = ~ .,
    cross_section_lags = 3,
    fast               = FALSE
  )
  fit_cpp <- dcce(
    data               = dcce_sim,
    unit_index         = "unit",
    time_index         = "time",
    formula            = y ~ L(y, 1) + x,
    model              = "dcce",
    cross_section_vars = ~ .,
    cross_section_lags = 3,
    fast               = TRUE
  )

  # Coefficients must match to 1e-6
  expect_equal(coef(fit_r), coef(fit_cpp), tolerance = 1e-6)
  # Standard errors must match too
  expect_equal(fit_r$se, fit_cpp$se, tolerance = 1e-6)
  # Residuals
  expect_equal(fit_r$residuals, fit_cpp$residuals, tolerance = 1e-6)
})


test_that("Rcpp and pure R agree on pwt8 DCCE", {
  data(pwt8)

  fit_r <- dcce(
    data               = pwt8,
    unit_index         = "country",
    time_index         = "year",
    formula            = d_log_rgdpo ~ L(log_rgdpo, 1) + log_hc + log_ck + log_ngd,
    model              = "dcce",
    cross_section_vars = ~ log_rgdpo + log_hc + log_ck + log_ngd,
    cross_section_lags = 3,
    fast               = FALSE
  )
  fit_cpp <- dcce(
    data               = pwt8,
    unit_index         = "country",
    time_index         = "year",
    formula            = d_log_rgdpo ~ L(log_rgdpo, 1) + log_hc + log_ck + log_ngd,
    model              = "dcce",
    cross_section_vars = ~ log_rgdpo + log_hc + log_ck + log_ngd,
    cross_section_lags = 3,
    fast               = TRUE
  )
  expect_equal(coef(fit_r), coef(fit_cpp), tolerance = 1e-6)
})


test_that("n_cores argument is accepted (silently ignored on Windows)", {
  data(dcce_sim)
  fit <- dcce(
    data               = dcce_sim,
    unit_index         = "unit",
    time_index         = "time",
    formula            = y ~ x,
    model              = "mg",
    cross_section_vars = NULL,
    n_cores            = 2L
  )
  expect_s3_class(fit, "dcce_mg_fit")
})


test_that(".unit_ols_cpp and .unit_ols return identical results on a single unit", {
  set.seed(42)
  T_val <- 40
  X <- cbind(1, rnorm(T_val), rnorm(T_val))
  colnames(X) <- c("(Intercept)", "x1", "x2")
  y <- 0.5 + 1.2 * X[, 2] - 0.4 * X[, 3] + rnorm(T_val, sd = 0.3)

  res_r   <- .unit_ols(y, X)
  if (.dcce_cpp_available()) {
    res_cpp <- .unit_ols_cpp(X, y)
    # C++ returns arma::vec wrapped as a column matrix; compare numerically.
    b_cpp <- as.numeric(res_cpp$b)
    names(b_cpp) <- colnames(X)
    expect_equal(res_r$b,  b_cpp, tolerance = 1e-8)
    expect_equal(res_r$r2, as.numeric(res_cpp$r2), tolerance = 1e-8)
  }
})
