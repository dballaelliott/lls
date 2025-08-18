# Script to generate data/info.sim.rda
# Run interactively or via devtools::load_all(); then source this file.

library(data.table)
# Simulate data 
set.seed(617)
n <- 500

# True APE is 1 - this is what we want to recover
tau <- runif(n)
tau <- tau / mean(tau)  # Normalize so mean is 1

# Make learning rates negatively correlated with belief effects
alpha <- 1 / (tau^2 + 1)
sigma2 <- 1 / (0.5 * tau + 1)

# Introduce endogeneity that affects both priors and outcomes
V <- rnorm(n)
U <- V + rnorm(n)

# Generate the experimental data
dt <- data.table(
  tau = tau,
  alpha = alpha,
  Z = runif(n) > 0.5
)

dt[, signal := ifelse(Z, 1, -1)]
dt[, prior := V + sqrt(sigma2) * rnorm(n) / 5]
dt[, posterior := alpha * (signal - prior) + prior]
dt[, Y := tau * posterior + U]
dt[, Y0 := tau * prior + U]

# Changes
dt[, dX := posterior - prior]
dt[, dY := Y - Y0]

info.sim <- dt

# Save using usethis for package data
if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(info.sim, overwrite = TRUE)
} else {
  dir.create("data", showWarnings = FALSE)
  save(info.sim, file = "data/info.sim.rda")
}
