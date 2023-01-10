test_that("Utility increases with lower hazard ratio", {
  expect_lte(utility_R2(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                          hr1 = 0.7, hr2 = 0.81, id1 = 280, id2 = 420,
                          xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                          c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                          K = Inf, N = Inf, S = -Inf,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000, fixed = TRUE)[1],
               utility_R2(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                          hr1 = 0.5, hr2 = 0.81, id1 = 280, id2 = 420,
                          xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                          c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                          K = Inf, N = Inf, S = -Inf,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000, fixed = TRUE)[1])
})

test_that("Utility increases with lower hazard ratio", {
  expect_lte(utility_R(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                          hr1 = 0.7, hr2 = 0.81, id1 = 280, id2 = 420,
                          xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                          c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                          K = Inf, N = Inf, S = -Inf,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000, fixed = FALSE)[1],
               utility_R(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                          hr1 = 0.5, hr2 = 0.81, id1 = 280, id2 = 420,
                          xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                          c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                          K = Inf, N = Inf, S = -Inf,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000, fixed = FALSE)[1])
})

test_that("Utility increases with lower hazard ratio", {
  expect_lte(utility_L2(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                         hr1 = 0.7, hr2 = 0.81, id1 = 280, id2 = 420,
                         xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                         c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                         K = Inf, N = Inf, S = -Inf,
                         steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                         b1 = 1000, b2 = 2000, b3 = 3000, fixed = TRUE)[1],
               utility_L2(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                         hr1 = 0.5, hr2 = 0.81, id1 = 280, id2 = 420,
                         xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                         c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                         K = Inf, N = Inf, S = -Inf,
                         steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                         b1 = 1000, b2 = 2000, b3 = 3000, fixed = TRUE)[1])
})

test_that("Utility increases with lower hazard ratio", {
  expect_lte(utility_L(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                          hr1 = 0.7, hr2 = 0.81, id1 = 280, id2 = 420,
                          xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                          c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                          K = Inf, N = Inf, S = -Inf,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000, fixed = FALSE)[1],
               utility_L(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3,
                          hr1 = 0.5, hr2 = 0.81, id1 = 280, id2 = 420,
                          xi2 = 0.7, xi3 = 0.7, alpha = 0.025, beta = 0.1,
                          c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                          K = Inf, N = Inf, S = -Inf,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000, fixed = FALSE)[1])
})
