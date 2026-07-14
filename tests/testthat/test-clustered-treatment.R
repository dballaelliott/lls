# rank normalization with tied / group-scalar treatments:
# ranks must span (0,1] within sign groups even when the same treatment value
# repeats across rows, and downstream aggregation must not depend on row counts

# helper: simulate a group-scalar (firm-level) dose stamped onto worker rows
sim.clustered <- function(n.firms = 20, workers = 25, seed = 42, negative = FALSE) {
  set.seed(seed)
  doses <- runif(n.firms, 0.1, 1)
  if (negative) doses <- doses * rep(c(-1, 1), length.out = n.firms)
  dose <- doses[rep(1:n.firms, each = workers)]
  n <- n.firms * workers
  data.frame(
    alpha = dose,
    posterior = dose + rnorm(n, sd = 0.1),
    Y = 2 * dose + rnorm(n)
  )
}

test_that("unique-per-row treatment: ranks match the row-count normalization exactly", {
  set.seed(1)
  n <- 300
  dose <- runif(n, 0.1, 1) # continuous, unique per row
  dat <- data.frame(alpha = dose, posterior = dose + rnorm(n, sd = 0.1), Y = 2 * dose + rnorm(n))

  est <- iv.lls(dat, y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1)

  # with unique values, max dense rank == .N, so the old and new normalizations
  # coincide bit-for-bit: ranks are (1:n)/n
  expect_equal(sort(unique(est$micro.dt$r)), (1:n) / n, tolerance = 0)
})

test_that("tied treatment: ranks span [1/K, 1]", {
  dat <- sim.clustered(n.firms = 20, workers = 25)

  est <- iv.lls(dat, y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1)

  r <- est$micro.dt$r
  expect_equal(max(r, na.rm = TRUE), 1)
  expect_equal(min(r[r > 0], na.rm = TRUE), 1 / 20)
  expect_true(is.finite(est$coef))
})

test_that("ranks and estimate are invariant to row duplication", {
  dat <- sim.clustered(n.firms = 20, workers = 25)

  est1 <- iv.lls(dat, y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1)
  est2 <- iv.lls(rbind(dat, dat), y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1)

  expect_setequal(unique(est1$micro.dt$r), unique(est2$micro.dt$r))
  expect_equal(est1$coef, est2$coef)
})

test_that("worker-level rows match firm-level collapse with weights = rows-per-firm", {
  set.seed(7)
  n.firms <- 27 # 12 control (dose 0) + 15 treated, distinct doses
  doses <- c(rep(0, 12), seq(0.2, 3, length.out = 15))
  sizes <- sample(5:40, n.firms, replace = TRUE)
  firm <- rep(1:n.firms, times = sizes)
  worker <- data.frame(dX = doses[firm])
  worker$dY <- 1.5 * worker$dX + rnorm(nrow(worker))

  # bandwidth/2 = 0.05 < 1/15 rank spacing: each window isolates a single dose
  est.worker <- panel.lls(worker, dy = "dY", dx = "dX",
                          normalize.x = TRUE, pointmass.zero = TRUE, bandwidth = 0.1)

  firm.dt <- aggregate(dY ~ dX + firm, data = cbind(worker, firm = firm), FUN = mean)
  firm.dt$n <- as.vector(table(firm))
  est.firm <- panel.lls(firm.dt, dy = "dY", dx = "dX", weights = "n",
                        normalize.x = TRUE, pointmass.zero = TRUE, bandwidth = 0.1)

  expect_equal(est.worker$coef, est.firm$coef, tolerance = 1e-8)
})

test_that("explicit r.support.points above the distinct-value count takes the weighted branch", {
  dat <- sim.clustered(n.firms = 20, workers = 25)

  est.default <- iv.lls(dat, y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1)
  # K + 1 grid points: no subsampling occurs, so this must be identical to the
  # default (row-weighted merge), not silently fall into the unweighted grid mean
  est.explicit <- iv.lls(dat, y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1,
                         r.support.points = 21)

  expect_equal(est.default$coef, est.explicit$coef, tolerance = 0)
})

test_that("repeated quantiles from point masses stay in the subsampled grid", {
  # one dose holds ~45% of rows: its rank must appear repeatedly among the grid
  # quantiles, and the grid mean must count that multiplicity (that is how the
  # quantile grid encodes row weights)
  set.seed(3)
  doses <- seq(0.5, 3, length.out = 11)
  dose <- doses[rep(1:11, times = c(rep(30, 5), 450, rep(30, 5)))]
  n <- length(dose)
  dat <- data.frame(alpha = dose, posterior = dose + rnorm(n, sd = 0.1),
                    Y = 2 * dose + rnorm(n))

  est <- iv.lls(dat, y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1,
                r.support.points = 5)

  grid <- est$micro.dt
  expect_equal(nrow(grid), 5) # grid rows, not data rows
  expect_true(sum(grid$r == 6 / 11) > 1) # heavy dose repeated, not deduplicated
  expect_equal(est$coef, mean(grid$tau_r, na.rm = TRUE)) # mean counts multiplicity
})

test_that("negative sign group with ties spans to -1", {
  dat <- sim.clustered(n.firms = 20, workers = 25, negative = TRUE)

  est <- iv.lls(dat, y = "Y", x = "posterior", r = "alpha", bandwidth = 0.1)

  r <- est$micro.dt$r
  expect_equal(min(r, na.rm = TRUE), -1)
  expect_equal(max(r, na.rm = TRUE), 1)
})
