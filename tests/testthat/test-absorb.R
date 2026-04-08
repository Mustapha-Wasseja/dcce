test_that(".absorb_demean_one removes a single-factor mean exactly", {
  set.seed(1)
  g <- factor(rep(1:5, each = 10))
  x <- rnorm(50) + as.numeric(g)   # group-specific means
  x_dm <- .absorb_demean_one(x, g)
  # After demeaning, group means must be 0
  means <- tapply(x_dm, g, mean)
  expect_true(all(abs(means) < 1e-12))
})


test_that(".absorb_demean_multi converges and removes two-way means", {
  set.seed(2)
  g1 <- factor(rep(1:6, each = 12))
  g2 <- factor(rep(1:12, times = 6))
  x  <- rnorm(72) + 2 * as.numeric(g1) - 0.5 * as.numeric(g2)
  groups <- list(g1 = g1, g2 = g2)
  x_dm <- .absorb_demean_multi(x, groups, tol = 1e-12)
  m1 <- tapply(x_dm, g1, mean)
  m2 <- tapply(x_dm, g2, mean)
  expect_true(all(abs(m1) < 1e-8))
  expect_true(all(abs(m2) < 1e-8))
})


test_that(".absorb_resolve accepts formulas and character vectors", {
  df <- data.frame(a = 1:10, b = factor(rep(1:2, 5)), c = factor(rep(1:5, 2)))
  attr(df, "unit_var") <- "b"
  attr(df, "time_var") <- "c"

  g1 <- .absorb_resolve(~ b, df)
  expect_true(is.list(g1))
  expect_equal(names(g1), "b")

  g2 <- .absorb_resolve(c("b", "c"), df)
  expect_equal(names(g2), c("b", "c"))

  expect_null(.absorb_resolve(NULL, df))
  expect_error(.absorb_resolve(~ notacol, df), "not found")
})


test_that("dcce() accepts absorb = NULL and matches the no-absorb call", {
  data(pwt8)
  fit1 <- dcce(data = pwt8, unit_index = "country", time_index = "year",
               formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
               model = "mg", cross_section_vars = NULL)
  fit2 <- dcce(data = pwt8, unit_index = "country", time_index = "year",
               formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
               model = "mg", cross_section_vars = NULL, absorb = NULL)
  expect_equal(coef(fit1), coef(fit2))
})


test_that("dcce() absorb = ~ year demeans time dummies from the regression", {
  data(pwt8)
  # Absorbing 'year' is equivalent to adding year fixed effects.
  # Residuals should still be defined; coefficients will differ from
  # the no-absorb MG.
  fit_raw <- dcce(data = pwt8, unit_index = "country", time_index = "year",
                  formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
                  model = "mg", cross_section_vars = NULL)
  fit_abs <- dcce(data = pwt8, unit_index = "country", time_index = "year",
                  formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
                  model = "mg", cross_section_vars = NULL,
                  absorb = ~ year)
  expect_s3_class(fit_abs, "dcce_fit")
  expect_true(is.finite(coef(fit_abs)["log_ck"]))
  # At least one coefficient should change after absorbing year FE
  expect_false(isTRUE(all.equal(coef(fit_raw), coef(fit_abs),
                                tolerance = 1e-8)))
})


test_that("dcce() absorb rejects a variable that does not exist", {
  data(pwt8)
  expect_error(
    dcce(data = pwt8, unit_index = "country", time_index = "year",
         formula = d_log_rgdpo ~ log_hc + log_ck,
         model = "mg", cross_section_vars = NULL,
         absorb = ~ not_a_column),
    regexp = "not found"
  )
})


test_that("dcce() absorb works with CCE estimator", {
  data(pwt8)
  fit <- dcce(data = pwt8, unit_index = "country", time_index = "year",
              formula = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
              model = "cce", cross_section_vars = ~ .,
              absorb = ~ year)
  expect_s3_class(fit, "dcce_cce_fit")
  expect_true(all(is.finite(coef(fit))))
})
