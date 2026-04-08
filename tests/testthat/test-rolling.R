test_that("dcce_rolling returns a dcce_rolling object with a coefficient tibble", {
  data(pwt8)
  roll <- dcce_rolling(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
    model      = "cce",
    cross_section_vars = ~ .,
    window     = 20L,
    step       = 4L
  )
  expect_s3_class(roll, "dcce_rolling")
  expect_true(!is.null(roll$coefficients))
  expect_s3_class(roll$coefficients, "tbl_df")
  expect_true(all(c("window_end", "window_start", "term", "estimate",
                    "std.error", "conf.low", "conf.high") %in%
                    names(roll$coefficients)))
})


test_that("dcce_rolling generates multiple windows with expected span", {
  data(pwt8)
  roll <- dcce_rolling(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "mg",
    cross_section_vars = NULL,
    window     = 15L,
    step       = 5L
  )
  expect_gt(roll$n_windows, 1L)
  # Window length should equal window - 1 in years (inclusive endpoints)
  diffs <- roll$coefficients$window_end - roll$coefficients$window_start
  expect_true(all(diffs == 14L))
})


test_that("dcce_rolling errors on window larger than panel", {
  data(dcce_sim)
  expect_error(
    dcce_rolling(
      data       = dcce_sim,
      unit_index = "unit",
      time_index = "time",
      formula    = y ~ x,
      model      = "mg",
      cross_section_vars = NULL,
      window     = 100L
    ),
    regexp = "time periods"
  )
})


test_that("dcce_rolling skips windows that fail and raises when none succeed", {
  data(dcce_sim)
  # Over-parameterised spec that will fail on narrow windows
  expect_error(
    dcce_rolling(
      data       = dcce_sim,
      unit_index = "unit",
      time_index = "time",
      formula    = y ~ L(y, 1) + x,
      model      = "dcce",
      cross_section_vars = ~ .,
      cross_section_lags = 10L,
      window     = 12L,
      step       = 5L
    )
  )
})


test_that("dcce_rolling print method works and coefficient path has content", {
  data(pwt8)
  roll <- dcce_rolling(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "mg",
    cross_section_vars = NULL,
    window     = 15L,
    step       = 5L
  )
  expect_output(print(roll), "Rolling-Window")
  expect_output(print(roll), "log_hc")
})


test_that("dcce_rolling works with model='dcce' when the spec fits", {
  data(pwt8)
  roll <- dcce_rolling(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ L(log_rgdpo, 1) + log_hc + log_ck,
    model      = "dcce",
    cross_section_vars = ~ .,
    cross_section_lags = 1L,
    window     = 25L,
    step       = 5L
  )
  expect_s3_class(roll, "dcce_rolling")
  expect_gt(roll$n_windows, 1L)
  # L(log_rgdpo,1) should appear in the coefficient table
  expect_true("L(log_rgdpo,1)" %in% roll$coefficients$term)
})
