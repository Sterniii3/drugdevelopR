test_that("multiplication works", {
  skip_on_cran()
  expect_equal(optimal_multiarm_binary( p0 = 0.6, p11 =  0.3, p12 = 0.5,
                                        n2min = 20, n2max = 100, stepn2 = 4, 
                                        rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,
                                        alpha = 0.05, beta = 0.1, 
                                        c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                        K = Inf, N = Inf, S = -Inf,
                                        steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                        b1 = 1000, b2 = 2000, b3 = 3000,
                                        strategy = 1, num_cl = 2)$u, 1671.27)
})
