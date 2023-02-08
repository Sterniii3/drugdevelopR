test_that("optimal_multiple_tte works for fixed = TRUE", {
  expect_equal(optimal_multiple_tte(hr1 = 0.75, hr2 = 0.80,  id1 = 210, id2 = 420,
                                    n2min = 30, n2max = 90, stepn2 = 6,
                                    hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,
                                    alpha = 0.05, beta = 0.1,
                                    c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                    K = Inf, N = Inf, S = -Inf,
                                    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                    b11 = 1000, b21 = 2000, b31 = 3000,
                                    b12 = 1000, b22 = 1500, b32 = 2000,
                                    rho = 0.6, fixed = TRUE, num_cl = 1)$u, 353)
})
