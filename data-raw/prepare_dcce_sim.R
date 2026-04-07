## Generate dcce_sim: synthetic panel with known DGP
## DGP: Chudik & Pesaran (2015), equation (1)

set.seed(20240101)

N <- 30L
T_val <- 50L

# True parameters (heterogeneous, drawn around means)
beta1 <- pmin(pmax(rnorm(N, mean = 0.50, sd = 0.10), -0.99), 0.99)
beta2 <- rnorm(N, mean = 1.00, sd = 0.20)
alpha <- rnorm(N, mean = 0.00, sd = 0.50)
gamma_i <- rnorm(N, mean = 1.00, sd = 0.30)
delta_i <- rnorm(N, mean = 1.00, sd = 0.30)
mu <- rnorm(N, mean = 0.00, sd = 0.30)

# Generate common factor (AR(1) with rho = 0.60)
rho_f <- 0.60
f <- numeric(T_val + 50L)
f[1] <- 0
for (t in 2:(T_val + 50L)) f[t] <- rho_f * f[t - 1] + rnorm(1, sd = 1)
f <- f[51:(T_val + 50L)]

# Generate x and y for each unit
records <- vector("list", N * T_val)
k <- 1L
for (i in seq_len(N)) {
  x <- numeric(T_val)
  y <- numeric(T_val)
  x[1] <- mu[i] + delta_i[i] * f[1] + rnorm(1, sd = 0.5)
  y[1] <- alpha[i] / (1 - beta1[i]) + beta2[i] * x[1] +
          gamma_i[i] * f[1] + rnorm(1, sd = 0.5)
  for (t in 2:T_val) {
    x[t] <- mu[i] + delta_i[i] * f[t] + rnorm(1, sd = 0.5)
    y[t] <- alpha[i] + beta1[i] * y[t - 1] + beta2[i] * x[t] +
            gamma_i[i] * f[t] + rnorm(1, sd = 0.5)
  }
  for (t in seq_len(T_val)) {
    records[[k]] <- data.frame(
      unit = i, time = t, y = y[t], x = x[t],
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

dcce_sim <- do.call(rbind, records)
rownames(dcce_sim) <- NULL

dcce_sim_truth <- list(
  beta1_mg = mean(beta1),
  beta2_mg = mean(beta2),
  N = N,
  T = T_val,
  seed = 20240101L
)

usethis::use_data(dcce_sim, overwrite = TRUE)
usethis::use_data(dcce_sim_truth, overwrite = TRUE)
