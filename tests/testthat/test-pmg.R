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

test_that(".pmg_recover_lr computes long-run coefficients", {
  # phi = -0.3 (speed of adjustment)
  # omega = c(0.15, 0.39) (coefficients on levels of x1, x2)
  # w = -omega/phi = c(0.5, 1.3)
  b_phi <- -0.3
  b_omega <- c(0.15, 0.39)
  V <- diag(3) * 0.01

  result <- .pmg_recover_lr(b_phi, b_omega, V,
                            phi_idx = 1L, omega_idx = 2:3)

  expect_equal(result$lr_coef, c(-0.15 / -0.3, -0.39 / -0.3), tolerance = 0.01)
  expect_true(is.matrix(result$lr_vcov))
  expect_equal(dim(result$lr_vcov), c(2, 2))
})

test_that(".csardl_recover_lr computes long-run from ARDL", {
  # ARDL(1,1): y_t = b1*y_{t-1} + b2_0*x_t + b2_1*x_{t-1}
  # LR = (b2_0 + b2_1) / (1 - b1)
  b <- c(y_L1 = 0.5, x_L0 = 0.2, x_L1 = 0.1)
  V <- diag(3) * 0.01
  rownames(V) <- colnames(V) <- names(b)

  result <- .csardl_recover_lr(
    b, py = 1, px = 1,
    y_lag_names = "y_L1",
    x_lag_names = c("x_L0", "x_L1"),
    V = V
  )

  # LR = (0.2 + 0.1) / (1 - 0.5) = 0.6
  expect_equal(unname(result$lr_coef), 0.6, tolerance = 0.01)
})
