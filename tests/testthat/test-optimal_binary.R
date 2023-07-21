test_that("Optimal_binary works with prior distribution", {
  skip_on_cran()
  expect_equal(optimal_binary(w = 0.3, p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 30, in2 = 60,  
                              n2min = 20, n2max = 100, stepn2 = 4,   
                              rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05, 
                              alpha = 0.05, beta = 0.1,   
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                              K = Inf, N = Inf, S = -Inf,  
                              steps1 = 1, stepm1 = 0.95,  stepl1 = 0.85,  
                              b1 = 1000, b2 = 2000, b3 = 3000,  
                              gamma = 0,  fixed = FALSE,
                              skipII = FALSE,num_cl = 1)[2], 
               data.frame(u=678.04))
})

test_that("Optimal_binary works with fixed effects", {
  skip_on_cran()
  expect_equal(optimal_binary(w = 0.3, p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 30, in2 = 60,  
                              n2min = 20, n2max = 100, stepn2 = 4,   
                              rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05, 
                              alpha = 0.05, beta = 0.1,   
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                              K = Inf, N = Inf, S = -Inf,  
                              steps1 = 1, stepm1 = 0.95,  stepl1 = 0.85,  
                              b1 = 1000, b2 = 2000, b3 = 3000,  
                              gamma = 0,  fixed = TRUE,
                              skipII = FALSE,num_cl = 1)[2], 
               data.frame(u=1806.86))
})

test_that("Optimal_binary works when skipping phase II", {
  skip_on_cran()
  expect_equal(optimal_binary(w = 0.3, p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 30, in2 = 60,  
                              n2min = 20, n2max = 100, stepn2 = 4,   
                              rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05, 
                              alpha = 0.05, beta = 0.1,   
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                              K = Inf, N = Inf, S = -Inf,  
                              steps1 = 1, stepm1 = 0.95,  stepl1 = 0.85,  
                              b1 = 1000, b2 = 2000, b3 = 3000,  
                              gamma = 0,  fixed = TRUE,
                              skipII = TRUE,num_cl = 1)[2, ]$u, 2234.78)
})

test_that("Optimal_binary works when skipping phase II with a prior distribution", {
  skip_on_cran()
  expect_equal(optimal_binary(w = 0.3, p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 30, in2 = 60,  
                              n2min = 20, n2max = 100, stepn2 = 4,   
                              rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05, 
                              alpha = 0.05, beta = 0.1,   
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                              K = Inf, N = Inf, S = -Inf,  
                              steps1 = 1, stepm1 = 0.95,  stepl1 = 0.85,  
                              b1 = 1000, b2 = 2000, b3 = 3000,  
                              gamma = 0,  fixed = FALSE,
                              skipII = TRUE,num_cl = 1)[2, ]$u, 1160.48)
})
