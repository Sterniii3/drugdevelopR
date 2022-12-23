#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("04.01", {
  # Multiarm time-to-event endpoints
  res <- optimal_multiarm(alpha = 0.025,
                          beta = 0.1
                          
               )
  expect_equal(res$u, , tolerance = 0.005)
  expect_equal(res$n2, )
  expect_equal(res$n3, )
  expect_equal(res$n, )
  expect_equal(res$d2, )
  expect_equal(res$d3, )
  expect_equal(res$d, )
  expect_equal(res$HRgo, )
})
