
test_that("more advanced arguments work", {
  data(info.sim, package = "lls")

  est <- panel.lls(info.sim, dy = "dY", dx = "dX", pointmass.zero = TRUE)

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))

  est <- panel.lls(info.sim, dy = "dY", dx = "dX", pointmass.zero = TRUE, r.support.points = 10)

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))


  est <- panel.lls(info.sim, dy = "dY", dx = "dX",
                   r.support.points = 10, bootstrap = TRUE, bootstrap.n = 10, n.cores = 1)

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))

})

# check setting sd.trim = 0

test_that("sd.trim = 0 works", {
  data(info.sim, package = "lls")

  est <- panel.lls(info.sim, dy = "dY", dx = "dX", sd.trim = 0, bootstrap.n = 3, bootstrap = TRUE)

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))
})

# check with both !normalize.r & !normalize.x
test_that("lls works with !normalize.r & !normalize.x", {
  data(info.sim, package = "lls")


  est <- iv.lls(info.sim, y = "Y", x = "posterior", r = "alpha",
                control.fml = "prior", bandwidth = 0.05,
                normalize.r = FALSE)

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))
})

#test FE with and without controls

test_that("FE works with & without controls", {
  data(info.sim, package = "lls")
  setDT(info.sim)  # ensure data.table)
  # cut into 3 prior bins
  info.sim[, prior_bin := cut(prior, breaks = c(-Inf, -0.5, 0.5, Inf), labels = c("low", "medium", "high"))]


  est <- iv.lls(info.sim, y = "Y", x = "posterior", r = "alpha",
                FE.fml = "prior_bin",
                control.fml = "prior_bin:prior", bandwidth = 0.05)

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))

  # and without
  est <- iv.lls(info.sim, y = "Y", x = "posterior", r = "alpha",
                FE.fml = "prior_bin")

  expect_s3_class(est, "lls")
  expect_true(is.numeric(est$coef))

})

test_that("clustering works", {
  data(info.sim, package = "lls")
  setDT(info.sim)  # ensure data.table)
  # cut into 3 prior bins

  info.sim[, prior_bin := cut(prior,breaks = 10)]


    est <- iv.lls(info.sim, y = "Y", x = "posterior", r = "alpha",
                  bandwidth = 0.05, cluster = "prior_bin", bootstrap= TRUE, bootstrap.n = 3)

    expect_s3_class(est, "lls")
    expect_true(is.numeric(est$coef))
    expect_true(is.numeric(est$se))

    est <- iv.lls(info.sim, y = "Y", x = "posterior", r = "alpha",
                  bandwidth = 0.05, cluster = "prior_bin", bootstrap= TRUE,
                  bootstrap.bayesian = FALSE,
                  bootstrap.n = 3,
                  n.cores = 1)

    expect_s3_class(est, "lls")
    expect_true(is.numeric(est$coef))
    expect_true(is.numeric(est$se))



})
