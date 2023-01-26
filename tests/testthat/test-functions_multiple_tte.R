test_that("Utility increases with lower hazard ratio", {
  expect_lte(utility_multiple_tte(n2 = 50, HRgo = 0.8, alpha = 0.025, beta = 0.1,
                                  hr1 = 0.75, hr2 = 0.80, id1 = 300, id2 = 600,
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                  K = Inf, N = Inf, S = -Inf, 
                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                  b11 = 1000, b21 = 2000, b31 = 3000,
                                  b12 = 1000, b22 = 1500, b32 = 2000,
                                  fixed = TRUE, rho = 0.3)[1], 
             utility_multiple_tte(n2 = 50, HRgo = 0.8, alpha = 0.025, beta = 0.1,
                                  hr1 = 0.5, hr2 = 0.6, id1 = 300, id2 = 600,
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                  K = Inf, N = Inf, S = -Inf, 
                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                  b11 = 1000, b21 = 2000, b31 = 3000,
                                  b12 = 1000, b22 = 1500, b32 = 2000,
                                  fixed = TRUE, rho = 0.3)[1])
})
