test_that("optimal_multiarm_normal works for strategy 1", {
  skip_on_cran()
  expect_equal(optimal_multiarm_normal(Delta1 = 0.375, Delta2 = 0.625, 
                                       n2min = 20, n2max = 100, stepn2 = 4,    
                                       kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,    
                                       alpha = 0.05, beta = 0.1,    
                                       c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,   
                                       K = Inf, N = Inf, S = -Inf,   
                                       steps1 = 0, stepm1 = 0.5, stepl1 = 0.8, 
                                       b1 = 3000, b2 = 8000, b3 = 10000,   
                                       strategy = 1, num_cl = 2)$u, 2924.13)
})

