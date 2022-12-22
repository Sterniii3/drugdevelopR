#' @editor Lukas D Sauer
#' @editDate 2022-12-20
test_that("03.01", {
  # Multi-trial time-to-event endpoints
  res <- optimal_multitrial(
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
