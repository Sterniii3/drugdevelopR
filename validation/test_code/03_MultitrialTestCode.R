#' @editor Lukas D Sauer
#' @editDate 2022-12-20
test_that("03.01", {
  # Multi-trial time-to-event endpoints
  res <- optimal_multitrial(alpha = 0.025,
                            beta = 0.1,
                            hr1 = 0.69, hr2 = 0.88,
                            xi2 = 0.7, xi3 = 0.7,
                            d2min = 10, bd2max = 400, stepd2 = 2,
                            hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01,
               )
  expect_equal(res$n2, 122)
  expect_equal(res$n3, 200)
  expect_equal(res$n, 322)
  expect_equal(res$u, 78, tolerance = 0.005)
  expect_equal(res$HRgo, 0.78)
  expect_equal(res$d2, 85)
  expect_equal(res$d3, 140)
  expect_equal(res$d, 225)
  
})
