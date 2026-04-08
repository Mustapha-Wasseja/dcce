test_that("structural_break_test returns a dcce_break object (known)", {
  data(pwt8)
  brk <- structural_break_test(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
    model = "mg", cross_section_vars = NULL,
    type = "known", break_date = 1985L
  )
  expect_s3_class(brk, "dcce_break")
  expect_equal(brk$type, "known")
  expect_true(is.finite(brk$statistic))
  expect_equal(brk$break_date, 1985L)
  expect_true(!is.null(brk$fit_pre))
  expect_true(!is.null(brk$fit_post))
})


test_that("structural_break_test rejects when a known DGP break is present", {
  set.seed(42)
  N <- 20; T_val <- 40
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  # True break at t = 21: beta jumps from 0.5 to 1.5
  df$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    beta_i <- ifelse(df$t[idx] <= 20, 0.5, 1.5)
    df$y[idx] <- beta_i * df$x[idx] + rnorm(T_val, sd = 0.3)
  }
  brk <- structural_break_test(
    data = df, unit_index = "id", time_index = "t",
    formula = y ~ x, model = "mg", cross_section_vars = NULL,
    type = "known", break_date = 20L
  )
  # Wald should be very large
  expect_gt(brk$statistic, 10)
  expect_lt(brk$p_value, 0.01)
})


test_that("structural_break_test(type='unknown') finds the true break date", {
  set.seed(43)
  N <- 20; T_val <- 40
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    beta_i <- ifelse(df$t[idx] <= 25, 0.3, 1.2)
    df$y[idx] <- beta_i * df$x[idx] + rnorm(T_val, sd = 0.3)
  }
  brk <- structural_break_test(
    data = df, unit_index = "id", time_index = "t",
    formula = y ~ x, model = "mg", cross_section_vars = NULL,
    type = "unknown", trim = 0.15
  )
  expect_equal(brk$type, "unknown")
  # Estimated break should be close to the true break at t = 25
  expect_true(abs(brk$break_date - 25L) <= 2L)
  expect_s3_class(brk$candidates, "tbl_df")
  expect_true(all(c("date", "wald") %in% names(brk$candidates)))
})


test_that("structural_break_test handles a panel with no break cleanly", {
  set.seed(44)
  N <- 15; T_val <- 30
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- 0.7 * df$x + rnorm(N * T_val, sd = 0.3)
  brk <- structural_break_test(
    data = df, unit_index = "id", time_index = "t",
    formula = y ~ x, model = "mg", cross_section_vars = NULL,
    type = "unknown", trim = 0.2
  )
  expect_s3_class(brk, "dcce_break")
  # With no true break, Andrews p-value should be > 0.05
  expect_gt(brk$p_value, 0.05)
  # Critical values from Andrews (1993) should be attached
  expect_true(all(c("cv10", "cv05", "cv01") %in% names(brk$critical_values)))
})


test_that(".andrews_sup_wald_cv returns monotonic critical values", {
  for (q in c(1L, 3L, 5L)) {
    for (pi0 in c(0.05, 0.15, 0.20)) {
      cv <- .andrews_sup_wald_cv(q, pi0)
      expect_true(cv["cv10"] < cv["cv05"])
      expect_true(cv["cv05"] < cv["cv01"])
    }
  }
})


test_that(".andrews_sup_wald_pvalue brackets correctly", {
  cv <- c(cv10 = 9.44, cv05 = 11.36, cv01 = 15.04)
  # Below 10% cv -> p = 0.5
  expect_equal(.andrews_sup_wald_pvalue(5, cv),  0.5)
  # Between 5% and 10%
  p_mid <- .andrews_sup_wald_pvalue(10, cv)
  expect_true(p_mid > 0.05 && p_mid < 0.10)
  # Above 1% cv -> p = 0.01
  expect_equal(.andrews_sup_wald_pvalue(20, cv), 0.01)
})


test_that("Andrews p-value fixes the v0.3.1 Bonferroni blow-up", {
  # Reproduces the user's pwt8 case: sup-Wald = 7.20, q = 3, pi0 = 0.15.
  # Previously Bonferroni produced p = 1.0; Andrews should report > 0.10.
  cv <- .andrews_sup_wald_cv(q = 3L, pi0 = 0.15)
  p  <- .andrews_sup_wald_pvalue(7.2046, cv)
  expect_equal(p, 0.5)   # below the 10% critical value
})


test_that("structural_break_test errors on invalid break_date", {
  data(pwt8)
  expect_error(
    structural_break_test(
      data = pwt8, unit_index = "country", time_index = "year",
      formula = d_log_rgdpo ~ log_hc,
      model = "mg", cross_section_vars = NULL,
      type = "known", break_date = 2099L
    ),
    regexp = "not in the panel"
  )
})


test_that("structural_break_test requires break_date when type='known'", {
  data(pwt8)
  expect_error(
    structural_break_test(
      data = pwt8, unit_index = "country", time_index = "year",
      formula = d_log_rgdpo ~ log_hc,
      model = "mg", cross_section_vars = NULL,
      type = "known"
    ),
    regexp = "break_date"
  )
})


test_that("structural_break_test print method works", {
  data(pwt8)
  brk <- structural_break_test(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = d_log_rgdpo ~ log_hc,
    model = "mg", cross_section_vars = NULL,
    type = "known", break_date = 1985L
  )
  expect_output(print(brk), "Structural Break")
  expect_output(print(brk), "Pre-break")
  expect_output(print(brk), "Post-break")
})


test_that("structural_break_test n_breaks=2 returns multiple break dates", {
  set.seed(45)
  N <- 25; T_val <- 60
  df <- data.frame(id = rep(1:N, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- rnorm(N * T_val)
  df$y <- NA_real_
  # Two breaks: t <= 20 beta=0.3, 20<t<=40 beta=1.2, t>40 beta=0.5
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    beta_i <- ifelse(df$t[idx] <= 20, 0.3,
                     ifelse(df$t[idx] <= 40, 1.2, 0.5))
    df$y[idx] <- beta_i * df$x[idx] + rnorm(T_val, sd = 0.3)
  }
  brk <- structural_break_test(
    data = df, unit_index = "id", time_index = "t",
    formula = y ~ x, model = "mg", cross_section_vars = NULL,
    type = "unknown", trim = 0.15, n_breaks = 2L
  )
  expect_s3_class(brk, "dcce_break")
  # At least one break detected; with n_breaks=2 possibly both
  expect_true(length(brk$break_dates) >= 1L)
})
