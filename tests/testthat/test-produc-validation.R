# Structural validation tests using Produc (US states panel)
# These test output structure and pipe compatibility on real data.
# Skip if plm is not available.

test_that("dcce works on Produc data (MG)", {
  skip_if_not_installed("plm")
  data("Produc", package = "plm")
  Produc$log_gsp  <- log(Produc$gsp)
  Produc$log_pcap <- log(Produc$pcap)
  Produc$log_pc   <- log(Produc$pc)
  Produc$log_emp  <- log(Produc$emp)
  Produc$state <- as.character(Produc$state)
  Produc$year  <- as.integer(as.character(Produc$year))

  fit <- dcce(
    data = Produc, unit_index = "state", time_index = "year",
    formula = log_gsp ~ log_pcap + log_pc + log_emp + unemp,
    model = "mg", cross_section_vars = NULL
  )

  expect_s3_class(fit, "dcce_mg_fit")
  expect_equal(fit$N, 48L)
  expect_equal(fit$n_obs, 816L)

  # log_emp should be positive and dominant
  expect_gt(coef(fit)["log_emp"], 0.5)
})

test_that("dcce CCE matches plm CCEMG on Produc", {
  skip_if_not_installed("plm")
  data("Produc", package = "plm")
  Produc$log_gsp  <- log(Produc$gsp)
  Produc$log_pcap <- log(Produc$pcap)
  Produc$log_pc   <- log(Produc$pc)
  Produc$log_emp  <- log(Produc$emp)
  Produc$state <- as.character(Produc$state)
  Produc$year  <- as.integer(as.character(Produc$year))

  fit_cce <- dcce(
    data = Produc, unit_index = "state", time_index = "year",
    formula = log_gsp ~ log_pcap + log_pc + log_emp + unemp,
    model = "cce", cross_section_vars = ~ ., cross_section_lags = 0
  )

  # pmg() internally calls plm() so we need the full namespace attached
  library(plm)
  plm_pdata <- pdata.frame(Produc, index = c("state", "year"))
  plm_fit <- pmg(
    log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
    data = plm_pdata, model = "cmg"
  )

  # dcce and plm CCEMG should match to 3 decimal places
  expect_equal(
    unname(coef(fit_cce)["log_emp"]),
    unname(coef(plm_fit)["log(emp)"]),
    tolerance = 0.001
  )
  expect_equal(
    unname(coef(fit_cce)["log_pcap"]),
    unname(coef(plm_fit)["log(pcap)"]),
    tolerance = 0.001
  )
})

test_that("tidy/glance output structure on Produc", {
  skip_if_not_installed("plm")
  data("Produc", package = "plm")
  Produc$log_gsp  <- log(Produc$gsp)
  Produc$log_pcap <- log(Produc$pcap)
  Produc$log_pc   <- log(Produc$pc)
  Produc$log_emp  <- log(Produc$emp)
  Produc$state <- as.character(Produc$state)
  Produc$year  <- as.integer(as.character(Produc$year))

  fit <- dcce(
    data = Produc, unit_index = "state", time_index = "year",
    formula = log_gsp ~ log_pcap + log_pc + log_emp + unemp,
    model = "cce", cross_section_vars = ~ ., cross_section_lags = 0
  )

  # tidy structure
  t_out <- tidy(fit)
  expect_s3_class(t_out, "tbl_df")
  expect_true(all(c("term", "estimate", "std.error", "statistic",
                     "p.value", "conf.low", "conf.high", "type") %in% names(t_out)))

  # glance structure
  g_out <- glance(fit)
  expect_s3_class(g_out, "tbl_df")
  expect_equal(g_out$n_units, 48L)
  expect_true(all(c("nobs", "n_units", "t_bar", "r.squared",
                     "estimator") %in% names(g_out)))

  # coef(unit) structure
  b_unit <- coef(fit, type = "unit")
  expect_s3_class(b_unit, "tbl_df")
  expect_true(all(c("unit", "term", "estimate") %in% names(b_unit)))
  # 48 states * 5 vars (intercept + 4 regressors)
  expect_equal(nrow(b_unit), 48L * 5L)
})

test_that("pipe chain works on Produc", {
  skip_if_not_installed("plm")
  data("Produc", package = "plm")
  Produc$log_gsp  <- log(Produc$gsp)
  Produc$log_pcap <- log(Produc$pcap)
  Produc$log_pc   <- log(Produc$pc)
  Produc$log_emp  <- log(Produc$emp)
  Produc$state <- as.character(Produc$state)
  Produc$year  <- as.integer(as.character(Produc$year))

  # fit |> tidy()
  result <- Produc |>
    dcce(unit_index = "state", time_index = "year",
         formula = log_gsp ~ log_pcap + log_pc + log_emp + unemp,
         model = "cce", cross_section_vars = ~ .) |>
    tidy()
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 5L)

  # fit |> pcd_test()
  cd_result <- Produc |>
    dcce(unit_index = "state", time_index = "year",
         formula = log_gsp ~ log_pcap + log_pc + log_emp + unemp,
         model = "cce", cross_section_vars = ~ .) |>
    pcd_test(test = "pesaran")
  expect_s3_class(cd_result, "dcce_cd")
})

test_that("CD test on CCE residuals shows low CSD on Produc", {
  skip_if_not_installed("plm")
  data("Produc", package = "plm")
  Produc$log_gsp  <- log(Produc$gsp)
  Produc$log_pcap <- log(Produc$pcap)
  Produc$log_pc   <- log(Produc$pc)
  Produc$log_emp  <- log(Produc$emp)
  Produc$state <- as.character(Produc$state)
  Produc$year  <- as.integer(as.character(Produc$year))

  fit_cce <- dcce(
    data = Produc, unit_index = "state", time_index = "year",
    formula = log_gsp ~ log_pcap + log_pc + log_emp + unemp,
    model = "cce", cross_section_vars = ~ ., cross_section_lags = 0
  )

  cd <- pcd_test(fit_cce, test = "pesaran")
  # After CCE, CD should not be extreme
  expect_true(abs(cd$statistics$statistic[1]) < 5)
})
