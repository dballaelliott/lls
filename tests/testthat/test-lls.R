# check that we can run with controls

test_that("lls runs on example data (IV mode)", {
  data(info.sim, package = "lls")

  est <- iv.lls(info.sim,
    y = "Y", x = "posterior", r = "alpha", control.fml = "prior",
    bandwidth = 0.05
  )

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))
})


test_that("lls runs on example data (panel mode)", {
  data(info.sim, package = "lls")

  est <- panel.lls(
    dat = info.sim,
    dx = "dX",           # Observed belief change
    dy = "dY"           # Observed outcome change
)

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))

  expect_output(print(est))
  expect_output(summary(est))

})


# check some of the more advanced options

test_that("bayesian bootstrap runs on example data", {
  data(info.sim, package = "lls")

  est <- lls(info.sim,
    y = "Y", x = "posterior", r = "alpha", mode = "iv", bandwidth = 0.05,
    bootstrap = TRUE, bootstrap.bayesian = TRUE, bootstrap.n = 2
  )

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))
})

test_that("classic bootstrap run on example data", {
  data(info.sim, package = "lls")

  est <- lls(info.sim,
    y = "Y", x = "posterior", r = "alpha", mode = "iv", bandwidth = 0.05,
    bootstrap = TRUE, bootstrap.bayesian = FALSE, bootstrap.n = 2, n.cores = 1
  )

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))
})
