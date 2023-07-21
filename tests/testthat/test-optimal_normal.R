test_that("Optimal normal works", {
  skip_on_cran()
  expect_equal(optimal_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                              in1=300, in2=600,  a = 0.25, b = 0.75,   
                              n2min = 20, n2max = 100, stepn2 = 4,  
                              kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,
                              alpha = 0.05, beta = 0.1,  
                              c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                              K = Inf, N = Inf, S = -Inf,   
                              steps1 = 0,  stepm1 = 0.5,stepl1 = 0.8,  
                              b1 = 3000, b2 = 8000, b3 = 10000,  
                              gamma = 0,  fixed = FALSE, skipII = FALSE,  
                              num_cl = 2)[2], 
               data.frame(u=2272.13))
})

test_that("Skipping Phase II works", {
  skip_on_cran()
 expect_equal(optimal_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                              in1=300, in2=600,  a = 0.25, b = 0.75,   
                              n2min = 20, n2max = 100, stepn2 = 4,  
                              kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,
                              alpha = 0.05, beta = 0.1,  
                              c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                              K = Inf, N = Inf, S = -Inf,   
                              steps1 = 0,  stepm1 = 0.5,stepl1 = 0.8,  
                              b1 = 3000, b2 = 8000, b3 = 10000,  
                              gamma = 0,  fixed = TRUE, skipII = TRUE,  
                              num_cl = 1)[2, ]$u, 2526.28)
})

test_that("Parameters for prior distribution are irrelevant if treatment effects are fixed", {
  skip_on_cran()
  expect_equal(optimal_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                              in1=300, in2=600,  a = 0.25, b = 0.75,   
                              n2min = 20, n2max = 100, stepn2 = 4,  
                              kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,
                              alpha = 0.05, beta = 0.1,  
                              c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                              K = Inf, N = Inf, S = -Inf,   
                              steps1 = 0,  stepm1 = 0.5,stepl1 = 0.8,  
                              b1 = 3000, b2 = 8000, b3 = 10000,  
                              gamma = 0,  fixed = TRUE, skipII = FALSE,  
                              num_cl = 2)[2], 
               optimal_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.695, 
                              in1=350, in2=700,  a = 0.25, b = 0.75,   
                              n2min = 20, n2max = 100, stepn2 = 4,  
                              kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,
                              alpha = 0.05, beta = 0.1,  
                              c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                              K = Inf, N = Inf, S = -Inf,   
                              steps1 = 0,  stepm1 = 0.5,stepl1 = 0.8,  
                              b1 = 3000, b2 = 8000, b3 = 10000,  
                              gamma = 0,  fixed = TRUE, skipII = FALSE,  
                              num_cl = 2)[2])
})

