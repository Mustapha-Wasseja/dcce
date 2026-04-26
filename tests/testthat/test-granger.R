test_that("granger_test returns a dcce_granger object with expected fields", {
  data(pwt8)
  gc <- granger_test(
    data = pwt8, unit_index = "country", time_index = "year",
    y = "d_log_rgdpo", x = "log_ck", lags = 1L
  )
  expect_s3_class(gc, "dcce_granger")
  expect_true(is.finite(gc$W_bar))
  expect_true(is.finite(gc$Z_bar))
  expect_true(is.finite(gc$p_value_Z))
  expect_true(gc$p_value_Z >= 0 && gc$p_value_Z <= 1)
  expect_equal(gc$lags, 1L)
  expect_equal(gc$y, "d_log_rgdpo")
  expect_equal(gc$x, "log_ck")
})


test_that("granger_test detects causality in a constructed DGP", {
  set.seed(100)
  N <- 25; T_val <- 40
  df <- data.frame(id = rep(1:N, each = T_val), t = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    df$y[idx[1]] <- rnorm(1)
    for (tt in 2:T_val) {
      df$y[idx[tt]] <- 0.3 * df$y[idx[tt - 1]] + 0.5 * df$x[idx[tt - 1]] + rnorm(1, sd = 0.3)
    }
  }
  gc <- granger_test(data = df, unit_index = "id", time_index = "t",
                     y = "y", x = "x", lags = 1L)
  # x truly Granger-causes y -> should reject
  expect_lt(gc$p_value_Z, 0.01)
})


test_that("granger_test does NOT reject when there is no causality", {
  set.seed(101)
  N <- 20; T_val <- 35
  df <- data.frame(id = rep(1:N, each = T_val), t = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- rnorm(N * T_val)
  gc <- granger_test(data = df, unit_index = "id", time_index = "t",
                     y = "y", x = "x", lags = 1L)
  # No causality -> should NOT reject at 5%
  expect_gt(gc$p_value_Z, 0.05)
})


test_that("granger_test works with multiple lags", {
  data(pwt8)
  gc <- granger_test(
    data = pwt8, unit_index = "country", time_index = "year",
    y = "d_log_rgdpo", x = "log_hc", lags = 3L
  )
  expect_s3_class(gc, "dcce_granger")
  expect_equal(gc$lags, 3L)
})


test_that("granger_test print method works", {
  data(pwt8)
  gc <- granger_test(
    data = pwt8, unit_index = "country", time_index = "year",
    y = "d_log_rgdpo", x = "log_ck", lags = 1L
  )
  expect_output(print(gc), "Granger")
  expect_output(print(gc), "W-bar")
  expect_output(print(gc), "Z-bar")
})


test_that("granger_test errors on non-existent variables", {
  data(pwt8)
  expect_error(
    granger_test(data = pwt8, unit_index = "country", time_index = "year",
                 y = "not_here", x = "log_ck"),
    regexp = "not found"
  )
})


test_that("granger_test Z-bar tilde is finite for pwt8", {
  data(pwt8)
  gc <- granger_test(
    data = pwt8, unit_index = "country", time_index = "year",
    y = "d_log_rgdpo", x = "log_ck", lags = 1L
  )
  expect_true(is.finite(gc$Z_bar_tilde))
  expect_true(is.finite(gc$p_value_Zt))
})
