# make sure we can handle weird arguments gracefully

test_that("lls handles edge cases in arguments", {
  data(info.sim, package = "lls")

  # Test with zero bandwidth
  est_zero_bandwidth <- lls(info.sim, y = "Y", x = "posterior", r = "alpha", mode = "iv", bandwidth = 0)
  expect_s3_class(est_zero_bandwidth, "lls")
  expect_true(is.numeric(est_zero_bandwidth$coef))

  # Test with negative bandwidth
  expect_error(
    lls(info.sim, y = "Y", x = "posterior", r = "alpha", mode = "iv", bandwidth = -0.05)
  )

})
