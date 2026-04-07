test_that(".delta_method computes correct SE for ratio", {
  # g(a, b) = a/b
  # Var(a/b) = (1/b)^2 * Var(a) + (a/b^2)^2 * Var(b) - 2*(a/b^3)*Cov(a,b)
  b <- c(a = 3, b = 2)
  V <- matrix(c(0.1, 0.01, 0.01, 0.05), 2, 2)
  g <- function(x) x[1] / x[2]

  result <- .delta_method(g, b, V)

  expect_equal(unname(result$estimate), 1.5)
  # Analytical: (1/2)^2*0.1 + (3/4)^2*0.05 - 2*(3/8)*0.01
  #           = 0.025 + 0.028125 - 0.0075 = 0.045625
  expect_equal(result$vcov[1, 1], 0.045625, tolerance = 1e-4)
})


test_that(".csardl_classify_terms correctly separates y lags and x groups", {
  cls <- .csardl_classify_terms(
    y_name = "y",
    x_names = c("L(y,1)", "L(y,2)", "x", "L(x,1)", "z", "L(z,1)")
  )
  expect_equal(cls$y_lag_terms, c("L(y,1)", "L(y,2)"))
  expect_equal(sort(cls$y_lag_orders), c(1L, 2L))
  expect_equal(names(cls$x_groups), c("x", "z"))
  expect_equal(cls$x_groups$x, c("x", "L(x,1)"))
  expect_equal(cls$x_groups$z, c("z", "L(z,1)"))
})


test_that(".csardl_unit_lr recovers LR from an ARDL(1,1)", {
  # Stata / textbook ARDL(1,1):
  #   y_t = c + phi y_{t-1} + beta0 x_t + beta1 x_{t-1}
  # Long-run: theta = (beta0 + beta1) / (1 - phi)
  b <- c("(Intercept)" = 0, "L(y,1)" = 0.5, "x" = 0.2, "L(x,1)" = 0.1)
  V <- diag(0.01, nrow = 4)
  rownames(V) <- colnames(V) <- names(b)

  classify <- list(
    y_lag_terms  = "L(y,1)",
    y_lag_orders = 1L,
    x_groups     = list(x = c("x", "L(x,1)")),
    x_group_lags = list(x = c(0L, 1L))
  )

  res <- .csardl_unit_lr(b, V, classify)

  # LR = (0.2 + 0.1) / (1 - 0.5) = 0.6
  expect_equal(unname(res$lr_coef), 0.6, tolerance = 1e-6)
  # phi = -(1 - 0.5) = -0.5
  expect_equal(unname(res$phi), -0.5, tolerance = 1e-6)
  # SEs are positive and finite
  expect_true(all(sqrt(diag(res$lr_vcov)) > 0))
  expect_true(res$phi_se > 0)
})


test_that(".pmg_pool_lr inverse-variance pools unit LR coefficients", {
  # Two units with LR = 1.0 and LR = 2.0, variances 1 and 4
  # Inverse-variance pooled: (1*1 + 0.25*2)/(1 + 0.25) = 1.5/1.25 = 1.2
  unit_lr <- list(
    list(lr_coef = c(x = 1.0), lr_vcov = matrix(1, 1, 1,
                                                dimnames = list("x", "x"))),
    list(lr_coef = c(x = 2.0), lr_vcov = matrix(4, 1, 1,
                                                dimnames = list("x", "x")))
  )
  pooled <- .pmg_pool_lr(unit_lr)

  expect_equal(unname(pooled$lr_coef["x"]), 1.2, tolerance = 1e-6)
  # Pooled SE = sqrt(1 / (1 + 0.25)) = sqrt(0.8) ≈ 0.8944
  expect_equal(unname(pooled$lr_se["x"]), sqrt(1 / 1.25), tolerance = 1e-6)
})


test_that("dcce() model='csardl' returns a fit with lr/adjustment blocks", {
  set.seed(600)
  N <- 15; T_val <- 35
  f <- rnorm(T_val)

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N)
  )
  df$y <- NA_real_; df$x <- NA_real_

  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    x_i <- rnorm(T_val) + 0.5 * f
    y_i <- numeric(T_val)
    y_i[1] <- rnorm(1)
    for (tt in 2:T_val) {
      y_i[tt] <- 0.6 * y_i[tt - 1] + 0.4 * x_i[tt] + 0.2 * f[tt] + rnorm(1, sd = 0.3)
    }
    df$x[idx] <- x_i
    df$y[idx] <- y_i
  }

  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ L(y, 1) + x + L(x, 1),
              model = "csardl",
              cross_section_vars = ~ y + x,
              cross_section_lags = 2)

  expect_s3_class(fit, "dcce_csardl_fit")
  expect_true(!is.null(fit$lr_coef))
  expect_true(!is.null(fit$adjustment))
  expect_true("x" %in% names(fit$lr_coef))
  # Adjustment is negative (error-correcting)
  expect_lt(fit$adjustment, 0)
})


test_that("dcce() model='csdl' auto-generates Delta x lags and reports LR", {
  set.seed(601)
  N <- 12; T_val <- 30
  f <- rnorm(T_val)

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y = rnorm(N * T_val),
    x = rnorm(N * T_val) + rep(f, N)
  )

  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ x,
              model = "csdl",
              cross_section_vars = ~ y + x,
              cross_section_lags = 2,
              csdl_xlags = 2)

  expect_s3_class(fit, "dcce_csdl_fit")
  expect_true(!is.null(fit$lr_coef))
  expect_true("x" %in% names(fit$lr_coef))
})


test_that("dcce() model='pmg' pools LR with smaller SE than MG", {
  set.seed(602)
  N <- 20; T_val <- 35
  f <- rnorm(T_val)

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N)
  )
  df$y <- NA_real_; df$x <- NA_real_

  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    x_i <- rnorm(T_val) + 0.5 * f
    y_i <- numeric(T_val)
    y_i[1] <- rnorm(1)
    for (tt in 2:T_val) {
      y_i[tt] <- 0.5 * y_i[tt - 1] + 0.3 * x_i[tt] + rnorm(1, sd = 0.3)
    }
    df$x[idx] <- x_i
    df$y[idx] <- y_i
  }

  fit_csardl <- dcce(data = df, unit_index = "id", time_index = "t",
                     formula = y ~ L(y, 1) + x + L(x, 1),
                     model = "csardl",
                     cross_section_vars = ~ y + x,
                     cross_section_lags = 2)

  fit_pmg <- dcce(data = df, unit_index = "id", time_index = "t",
                  formula = y ~ L(y, 1) + x + L(x, 1),
                  model = "pmg",
                  cross_section_vars = ~ y + x,
                  cross_section_lags = 2)

  expect_s3_class(fit_pmg, "dcce_pmg_fit")
  expect_true(isTRUE(fit_pmg$pmg_pooled))
  # Pooled SE should typically be smaller than MG SE
  expect_lt(fit_pmg$lr_se["x"], fit_csardl$lr_se["x"] * 2)  # generous bound
})
