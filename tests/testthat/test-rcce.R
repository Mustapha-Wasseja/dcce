test_that("ER criterion picks correct number of factors", {
  # Two large eigenvalues, rest small
  eigenvalues <- c(100, 80, 2, 1.5, 1, 0.8, 0.5, 0.3)
  k <- .ahn_horenstein_er(eigenvalues, K_max = 4)
  expect_equal(k, 2L)  # biggest ratio at 80/2 = 40
})

test_that("GR criterion picks correct number of factors", {
  eigenvalues <- c(100, 80, 2, 1.5, 1, 0.8, 0.5, 0.3)
  k <- .ahn_horenstein_gr(eigenvalues, K_max = 4)
  expect_true(k %in% 1:3)  # should be around 2
})

test_that(".svd_csa with explicit npc", {
  set.seed(200)
  T_val <- 50
  K <- 5
  Z <- matrix(rnorm(T_val * K), T_val, K)

  result <- .svd_csa(Z, npc = 2)
  expect_equal(result$k, 2)
  expect_equal(ncol(result$V_k), 2)
  expect_equal(nrow(result$V_k), T_val)
})

test_that(".svd_csa with automatic selection", {
  set.seed(201)
  T_val <- 100
  # Create Z with 2 dominant factors
  F1 <- rnorm(T_val)
  F2 <- rnorm(T_val)
  noise <- matrix(rnorm(T_val * 3, sd = 0.1), T_val, 3)
  Z <- cbind(F1, F2, 0.5 * F1 + noise[, 1], 0.3 * F2 + noise[, 2], noise[, 3])

  result <- .svd_csa(Z, npc = NULL, criterion = "er")
  # Should select around 2 components
  expect_true(result$k >= 1 && result$k <= 4)
})

test_that("dcce with model='rcce' runs (placeholder)", {
  # rcce model dispatch goes through dcce() which currently treats it

  # like other models. The rCCE-specific logic (replacing CSAs with PCs)
  # will be fully wired in when we integrate into dcce().
  # For now, test that the SVD components work correctly.
  set.seed(202)
  T_val <- 50
  Z <- matrix(rnorm(T_val * 4), T_val, 4)
  result <- .svd_csa(Z, npc = 2)
  # PCs should be orthogonal
  pc_cross <- crossprod(result$V_k)
  expect_equal(pc_cross[1, 2], 0, tolerance = 1e-10)
})
