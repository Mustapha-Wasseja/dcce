test_that("CDw+ test returns a finite statistic and valid p-value", {
  set.seed(10)
  N <- 15; T_val <- 30
  f <- rnorm(T_val)
  em <- matrix(NA_real_, N, T_val)
  for (i in seq_len(N)) em[i, ] <- 0.8 * f + rnorm(T_val, sd = 0.3)

  res <- pcd_test(em, test = "cdwplus", n_reps = 49)
  expect_s3_class(res, "dcce_cd")
  stat <- res$statistics$statistic[res$statistics$test == "cdwplus"]
  p <- res$statistics$p_value[res$statistics$test == "cdwplus"]
  expect_true(is.finite(stat))
  expect_true(p >= 0 && p <= 1)
})


test_that("CIPS test runs and returns a valid statistic", {
  set.seed(11)
  N <- 20; T_val <- 40
  X <- matrix(NA, N, T_val)
  for (i in seq_len(N)) X[i, ] <- cumsum(rnorm(T_val))  # random walks

  res <- cips_test(X, lags = 1)
  expect_s3_class(res, "dcce_cips")
  expect_true(is.finite(res$statistic))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  # Random walks: should NOT reject the unit root null
  expect_gt(res$p_value, 0.05)
})


test_that("CIPS rejects for stationary series", {
  set.seed(12)
  N <- 25; T_val <- 60
  X <- matrix(rnorm(N * T_val), N, T_val)  # white noise (stationary)
  res <- cips_test(X, lags = 0)
  expect_s3_class(res, "dcce_cips")
  # Should reject unit root for stationary series
  expect_lt(res$p_value, 0.10)
})


test_that("Swamy test detects heterogeneity in a heterogeneous DGP", {
  set.seed(13)
  N <- 20; T_val <- 30
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- NA_real_

  # Heterogeneous slopes: beta_i ~ N(1, 0.5)
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    beta_i <- rnorm(1, 1, 0.5)
    df$y[idx] <- beta_i * df$x[idx] + rnorm(T_val, sd = 0.3)
  }

  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ x, model = "mg", cross_section_vars = NULL)

  sw <- swamy_test(fit)
  expect_s3_class(sw, "dcce_swamy")
  expect_true(is.finite(sw$S_stat))
  expect_true(is.finite(sw$delta_stat))
  # Should reject homogeneity
  expect_lt(sw$p_swamy, 0.05)
})


test_that("Hausman test runs on an MG fit", {
  set.seed(14)
  N <- 20; T_val <- 30
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- df$x + rnorm(N * T_val)

  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ x, model = "mg", cross_section_vars = NULL)

  h <- hausman_test(fit)
  expect_s3_class(h, "dcce_hausman")
  expect_true(is.finite(h$statistic))
  expect_true(h$p_value >= 0 && h$p_value <= 1)
})


test_that("confint.dcce_fit handles mg/lr/adjustment types", {
  set.seed(15)
  N <- 15; T_val <- 30
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- numeric(N * T_val)
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    for (tt in 2:T_val) {
      df$y[idx[tt]] <- 0.5 * df$y[idx[tt - 1]] + 0.3 * df$x[idx[tt]] + rnorm(1)
    }
  }

  fit_mg <- dcce(data = df, unit_index = "id", time_index = "t",
                 formula = y ~ x, model = "mg", cross_section_vars = NULL)
  ci_mg <- confint(fit_mg, type = "mg")
  expect_true(is.matrix(ci_mg))
  expect_equal(ncol(ci_mg), 2L)

  fit_csardl <- dcce(data = df, unit_index = "id", time_index = "t",
                     formula = y ~ L(y, 1) + x + L(x, 1),
                     model = "csardl",
                     cross_section_vars = ~ y + x,
                     cross_section_lags = 2)
  ci_lr <- confint(fit_csardl, type = "lr")
  expect_true(is.matrix(ci_lr))
  ci_adj <- confint(fit_csardl, type = "adjustment")
  expect_equal(nrow(ci_adj), 1L)
})


test_that("tidy() includes adjustment and lr rows for CS-ARDL fits", {
  set.seed(16)
  N <- 15; T_val <- 30
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- numeric(N * T_val)
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    for (tt in 2:T_val) {
      df$y[idx[tt]] <- 0.5 * df$y[idx[tt - 1]] + 0.3 * df$x[idx[tt]] + rnorm(1)
    }
  }

  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ L(y, 1) + x + L(x, 1),
              model = "csardl",
              cross_section_vars = ~ y + x,
              cross_section_lags = 2)
  td <- tidy(fit)
  expect_true("adjustment" %in% td$type)
  expect_true("lr" %in% td$type)
})
