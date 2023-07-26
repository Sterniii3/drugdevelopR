test_that("Optimal_multiarm works", {
  skip_on_cran()
  expect_equal(optimal_multiarm(hr1 = 0.75, hr2 = 0.80, ec = 0.6,   
                                n2min = 30, n2max = 90, stepn2 = 6,   
                                hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,  
                                alpha = 0.05, beta = 0.1,  
                                c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,     
                                K = Inf, N = Inf, S = -Inf,   
                                steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,   
                                b1 = 1000, b2 = 2000, b3 = 3000,  
                                strategy = 1,  num_cl = 2)[2], 
               drugdevelopResult(data.frame(u=241.46)))
})
