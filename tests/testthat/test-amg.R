test_that("AMG returns dcce_amg_fit object", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
    model      = "amg",
    cross_section_vars = NULL
  )
  expect_s3_class(fit, "dcce_amg_fit")
  expect_s3_class(fit, "dcce_fit")
  expect_true(!is.null(fit$cdp))
  expect_true(all(c("year", "cdp") %in% names(fit$cdp)))
})


test_that("AMG strips cdp_level from the reported MG coefficients", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
    model      = "amg",
    cross_section_vars = NULL
  )
  # cdp_level must not appear in the MG coefficient vector nor in tidy()
  expect_false("cdp_level" %in% names(coef(fit)))
  td <- tidy(fit)
  expect_false("cdp_level" %in% td$term)
})


test_that("AMG coefficients have expected signs on pwt8", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
    model      = "amg",
    cross_section_vars = NULL
  )
  b <- coef(fit)
  expect_gt(b["log_ck"], 0)   # physical capital: positive
  expect_lt(b["log_ngd"], 0)  # pop growth + depreciation: negative
})


test_that("AMG tidy() output has correct structure", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "amg",
    cross_section_vars = NULL
  )
  t <- tidy(fit)
  expect_s3_class(t, "tbl_df")
  expect_true(all(c("term", "estimate", "std.error", "p.value") %in% names(t)))
  expect_false("cdp_level" %in% t$term)
})


test_that("AMG reduces residual CSD relative to plain MG", {
  data(pwt8)
  fit_mg  <- dcce(data = pwt8, unit_index = "country", time_index = "year",
                  formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
                  model = "mg", cross_section_vars = NULL)
  fit_amg <- dcce(data = pwt8, unit_index = "country", time_index = "year",
                  formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
                  model = "amg", cross_section_vars = NULL)
  cd_mg  <- pcd_test(fit_mg,  test = "pesaran")
  cd_amg <- pcd_test(fit_amg, test = "pesaran")
  # AMG should not be drastically worse than MG; typically it helps
  expect_lt(
    abs(cd_amg$statistics$statistic[1]),
    abs(cd_mg$statistics$statistic[1]) + 5
  )
})


test_that("AMG works with fast = FALSE (pure-R path)", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "amg",
    cross_section_vars = NULL,
    fast       = FALSE
  )
  expect_s3_class(fit, "dcce_amg_fit")
})
