# Tests for v0.4.0 features: IFE, CCEP, half-panel JK, IPS/LLC, Pedroni/Kao, IRF

# ── IFE (Bai 2009) ──────────────────────────────────────────────────────────

test_that("dcce(model='ife') returns an IFE fit", {
  set.seed(200)
  N <- 20; T_val <- 30
  f <- rnorm(T_val)
  df <- data.frame(id = rep(1:N, each = T_val), t = rep(1:T_val, N))
  df$x <- NA_real_; df$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    li <- rnorm(1)
    df$x[idx] <- 0.5 * f + rnorm(T_val)
    df$y[idx] <- 0.6 * df$x[idx] + li * f + rnorm(T_val, sd = 0.3)
  }
  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ x, model = "ife", n_factors = 1L)
  expect_s3_class(fit, "dcce_ife_fit")
  expect_s3_class(fit, "dcce_fit")
  expect_equal(fit$n_factors, 1L)
  expect_true(abs(coef(fit)["x"] - 0.6) < 0.2)
})


test_that("IFE print method works", {
  set.seed(201)
  df <- data.frame(id = rep(1:10, each = 20), t = rep(1:20, 10),
                   y = rnorm(200), x = rnorm(200))
  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ x, model = "ife", n_factors = 1L)
  expect_output(print(fit), "Interactive Fixed Effects")
})


# ── CCEP (Pooled CCE) ────────────────────────────────────────────────────────

test_that("dcce(model='ccep') returns a pooled CCE fit", {
  data(pwt8)
  fit_ccep <- dcce(data = pwt8, unit_index = "country", time_index = "year",
                   formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
                   model = "ccep", cross_section_vars = ~ .)
  expect_s3_class(fit_ccep, "dcce_ccep_fit")

  # CCEP produces a valid coefficient vector
  b_ccep <- coef(fit_ccep)
  expect_true(all(is.finite(b_ccep)))
  expect_true("log_ck" %in% names(b_ccep))
  # SE should be smaller than or comparable to CCE-MG SE (more efficient)
  fit_cce <- dcce(data = pwt8, unit_index = "country", time_index = "year",
                  formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
                  model = "cce", cross_section_vars = ~ .)
  expect_lte(fit_ccep$se["log_ck"], fit_cce$se["log_ck"] * 1.5)
})


test_that("CCEP tidy() output has correct structure", {
  data(pwt8)
  fit <- dcce(data = pwt8, unit_index = "country", time_index = "year",
              formula = d_log_rgdpo ~ log_hc + log_ck,
              model = "ccep", cross_section_vars = ~ .)
  td <- tidy(fit)
  expect_s3_class(td, "tbl_df")
  expect_true(all(c("term", "estimate", "std.error") %in% names(td)))
})


# ── Half-panel jackknife ─────────────────────────────────────────────────────

test_that("half_panel_jackknife bias correction runs without error", {
  data(dcce_sim)
  fit_plain <- dcce(data = dcce_sim, unit_index = "unit", time_index = "time",
                    formula = y ~ L(y, 1) + x, model = "dcce",
                    cross_section_vars = ~ ., cross_section_lags = 3,
                    bias_correction = "none")
  fit_hpj <- dcce(data = dcce_sim, unit_index = "unit", time_index = "time",
                  formula = y ~ L(y, 1) + x, model = "dcce",
                  cross_section_vars = ~ ., cross_section_lags = 3,
                  bias_correction = "half_panel_jackknife")
  expect_s3_class(fit_hpj, "dcce_dcce_fit")
  # HPJ should produce different coefficients from plain
  expect_false(isTRUE(all.equal(coef(fit_plain), coef(fit_hpj), tolerance = 1e-6)))
})


# ── IPS and LLC panel unit root ──────────────────────────────────────────────

test_that("panel_ur_test(test='ips') runs and returns dcce_unit_root", {
  set.seed(300)
  X <- matrix(cumsum(rnorm(20 * 30)), 20, 30)
  res <- panel_ur_test(X, test = "ips", lags = 1)
  expect_s3_class(res, "dcce_unit_root")
  expect_equal(res$test, "IPS")
  expect_true(is.finite(res$statistic))
})


test_that("panel_ur_test(test='llc') runs and returns dcce_unit_root", {
  set.seed(301)
  X <- matrix(rnorm(25 * 40), 25, 40)   # stationary
  res <- panel_ur_test(X, test = "llc", lags = 0)
  expect_s3_class(res, "dcce_unit_root")
  expect_equal(res$test, "LLC")
  # Stationary series should reject
  expect_lt(res$p_value, 0.10)
})


test_that("IPS does NOT reject for random walks", {
  set.seed(302)
  X <- matrix(NA, 20, 40)
  for (i in 1:20) X[i, ] <- cumsum(rnorm(40))
  res <- panel_ur_test(X, test = "ips", lags = 1)
  expect_gt(res$p_value, 0.05)
})


# ── Pedroni and Kao cointegration ────────────────────────────────────────────

test_that("panel_coint_test(test='pedroni') returns expected structure", {
  data(pwt8)
  res <- panel_coint_test(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc + log_ck, test = "pedroni", lags = 1
  )
  expect_s3_class(res, "dcce_cointegration_extra")
  expect_equal(res$test, "PEDRONI")
  expect_true(nrow(res$statistics) == 2L)
  expect_true(all(c("test", "statistic", "p_value") %in% names(res$statistics)))
})


test_that("panel_coint_test(test='kao') returns expected structure", {
  data(pwt8)
  res <- panel_coint_test(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc + log_ck, test = "kao", lags = 1
  )
  expect_s3_class(res, "dcce_cointegration_extra")
  expect_equal(res$test, "KAO")
  expect_true(nrow(res$statistics) == 1L)
})


test_that("panel_coint_test print method works", {
  data(pwt8)
  res <- panel_coint_test(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc, test = "pedroni", lags = 1
  )
  expect_output(print(res), "PEDRONI")
  expect_output(print(res), "no cointegration")
})


# ── Impulse Response Functions ───────────────────────────────────────────────

test_that("irf() returns a dcce_irf object with correct length", {
  data(dcce_sim)
  fit <- dcce(data = dcce_sim, unit_index = "unit", time_index = "time",
              formula = y ~ L(y, 1) + x,
              model = "dcce", cross_section_vars = ~ .,
              cross_section_lags = 3)
  ir <- irf(fit, impulse = "x", horizon = 10, boot_reps = 0)
  expect_s3_class(ir, "dcce_irf")
  expect_equal(length(ir$irf), 11L)
  expect_equal(ir$impulse, "x")
  expect_true(abs(ir$irf[1]) > 0)   # immediate response is non-zero
})


test_that("irf() with bootstrap produces confidence bands", {
  data(dcce_sim)
  fit <- dcce(data = dcce_sim, unit_index = "unit", time_index = "time",
              formula = y ~ L(y, 1) + x,
              model = "dcce", cross_section_vars = ~ .,
              cross_section_lags = 3)
  ir <- irf(fit, impulse = "x", horizon = 5, boot_reps = 49, seed = 42)
  expect_true(!is.null(ir$lower))
  expect_true(!is.null(ir$upper))
  expect_equal(length(ir$lower), 6L)
  # Lower should be below upper at each horizon
  expect_true(all(ir$lower <= ir$upper))
})


test_that("irf() errors on non-dynamic models", {
  data(pwt8)
  fit <- dcce(data = pwt8, unit_index = "country", time_index = "year",
              formula = d_log_rgdpo ~ log_hc + log_ck,
              model = "cce", cross_section_vars = ~ .)
  expect_error(irf(fit, impulse = "log_hc"), regexp = "dynamic")
})


test_that("irf() print method works", {
  data(dcce_sim)
  fit <- dcce(data = dcce_sim, unit_index = "unit", time_index = "time",
              formula = y ~ L(y, 1) + x,
              model = "dcce", cross_section_vars = ~ .,
              cross_section_lags = 3)
  ir <- irf(fit, impulse = "x", horizon = 5, boot_reps = 0)
  expect_output(print(ir), "Impulse Response")
})
