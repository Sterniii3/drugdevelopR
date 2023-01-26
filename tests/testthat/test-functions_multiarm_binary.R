test_that("Utility increases for lower risk ratio", {
  expect_lte(utility_multiarm_binary(n2 = 50, RRgo = 0.8,
                                     alpha = 0.05, beta = 0.1, 
                                     p0 = 0.6, p11 =  0.3, p12 = 0.5, strategy = 1,
                                     c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                     K = Inf, N = Inf, S = -Inf, 
                                     steps1 = 1, stepm1 = 0.95,   stepl1 = 0.85,
                                     b1 = 1000, b2 = 2000, b3 = 3000)[1],
             utility_multiarm_binary(n2 = 50, RRgo = 0.8,
                                     alpha = 0.05, beta = 0.1, 
                                     p0 = 0.8, p11 =  0.2, p12 = 0.5, strategy = 1,
                                     c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                     K = Inf, N = Inf, S = -Inf, 
                                     steps1 = 1, stepm1 = 0.95,   stepl1 = 0.85,
                                     b1 = 1000, b2 = 2000, b3 = 3000)[1] )
})

test_that("Utility increases for lower risk ratio", {
  expect_lte(utility_multiarm_binary(n2 = 50, RRgo = 0.8,
                                     alpha = 0.05, beta = 0.1, 
                                     p0 = 0.6, p11 =  0.3, p12 = 0.5, strategy = 2,
                                     c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                     K = Inf, N = Inf, S = -Inf, 
                                     steps1 = 1, stepm1 = 0.95,   stepl1 = 0.85,
                                     b1 = 1000, b2 = 2000, b3 = 3000)[1],
             utility_multiarm_binary(n2 = 50, RRgo = 0.8,
                                     alpha = 0.05, beta = 0.1, 
                                     p0 = 0.8, p11 =  0.2, p12 = 0.5, strategy = 2,
                                     c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                     K = Inf, N = Inf, S = -Inf, 
                                     steps1 = 1, stepm1 = 0.95,   stepl1 = 0.85,
                                     b1 = 1000, b2 = 2000, b3 = 3000)[1] )
})
