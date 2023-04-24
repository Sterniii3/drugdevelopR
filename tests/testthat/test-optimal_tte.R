test_that("Optimal_tte works with prior distribution", {
  skip_on_cran()
  expect_equal(optimal_tte(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                           id1 = 210, id2 = 420,     
                           d2min = 20, d2max = 100, stepd2 = 5,    
                           hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,   
                           alpha = 0.05, beta = 0.1, xi2 = 0.7, xi3 = 0.7,  
                           c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                           K = Inf, N = Inf, S = -Inf,   
                           steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,   
                           b1 = 1000, b2 = 2000, b3 = 3000,   gamma = 0, 
                           fixed = FALSE, skipII = FALSE, num_cl = 1)[1], 
               data.frame(u=164.01))
})

test_that("Optimal_tte works without distribution", {
  skip_on_cran()
  expect_equal(optimal_tte(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                           id1 = 210, id2 = 420,     
                           d2min = 20, d2max = 100, stepd2 = 5,    
                           hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,   
                           alpha = 0.05, beta = 0.1, xi2 = 0.7, xi3 = 0.7,  
                           c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                           K = Inf, N = Inf, S = -Inf,   
                           steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,   
                           b1 = 1000, b2 = 2000, b3 = 3000,   gamma = 0, 
                           fixed = TRUE, skipII = FALSE, num_cl = 1)[1], 
               data.frame(u=1020.73))
})

test_that("Optimal_tte works when skipping phase II", {
  skip_on_cran()
  expect_equal(optimal_tte(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                           id1 = 210, id2 = 420,     
                           d2min = 20, d2max = 100, stepd2 = 5,    
                           hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.01,   
                           alpha = 0.05, beta = 0.1, xi2 = 0.7, xi3 = 0.7,  
                           c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                           K = Inf, N = Inf, S = -Inf,   
                           steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,   
                           b1 = 1000, b2 = 2000, b3 = 3000,   gamma = 0, 
                           fixed = TRUE, skipII = TRUE, num_cl = 1)[[2]]$u, 
               1703.7)
})