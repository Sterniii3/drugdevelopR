test_that("Optimal_bias works for additive adjustment method", {
  skip_on_cran()
  expect_equal(optimal_bias(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                            id1 = 210, id2 = 420,     
                            d2min = 20, d2max = 100, stepd2 = 5, 
                            hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,     
                            adj = "additive",   
                            lambdamin = 0.6, lambdamax = 1, steplambda = 0.05,  
                            alphaCImin = 0.25, alphaCImax = 0.5, stepalphaCI = 0.025,  
                            alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                            K = Inf, N = Inf, S = -Inf,   
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,  
                            b1 = 1000, b2 = 2000, b3 = 3000,  
                            fixed = TRUE, num_cl = 2)$u, 859.71)
})

test_that("Optimal_bias works for multiplicative method", {
  skip_on_cran()
  expect_equal(optimal_bias(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                            id1 = 210, id2 = 420,     
                            d2min = 20, d2max = 100, stepd2 = 5, 
                            hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,     
                            adj = "multiplicative",   
                            lambdamin = 0.6, lambdamax = 1, steplambda = 0.05,  
                            alphaCImin = 0.025, alphaCImax = 0.5, stepalphaCI = 0.025,  
                            alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                            K = Inf, N = Inf, S = -Inf,   
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,  
                            b1 = 1000, b2 = 2000, b3 = 3000,  
                            fixed = FALSE, num_cl = 2)$u, 98.28)
})

test_that("Optimal_bias works for both methods", {
  skip_on_cran()
  expect_equal(optimal_bias(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                            id1 = 210, id2 = 420,     
                            d2min = 50, d2max = 150, stepd2 = 5, 
                            hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,     
                            adj = "both",   
                            lambdamin = 0.6, lambdamax = 1, steplambda = 0.05,  
                            alphaCImin = 0.025, alphaCImax = 0.5, stepalphaCI = 0.025,  
                            alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                            K = Inf, N = Inf, S = -Inf,   
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,  
                            b1 = 1000, b2 = 2000, b3 = 3000,  
                            fixed = TRUE, num_cl = 2)$u, c(972.87,922.63))
})

test_that("Optimal_bias works for method `all`", {
  skip_on_cran()
  expect_equal(optimal_bias(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                            id1 = 210, id2 = 420,     
                            d2min = 50, d2max = 150, stepd2 = 5, 
                            hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,     
                            adj = "all",   
                            lambdamin = 0.6, lambdamax = 1, steplambda = 0.05,  
                            alphaCImin = 0.25, alphaCImax = 0.5, stepalphaCI = 0.025,  
                            alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                            K = Inf, N = Inf, S = -Inf,   
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,  
                            b1 = 1000, b2 = 2000, b3 = 3000,  
                            fixed = TRUE, num_cl = 2)$u, c(972.87,922.63,966.59,921.11))
})

test_that("Optimal_bias works for method `all` with prior distributions", {
  skip_on_cran()
  expect_equal(optimal_bias(w = 0.3, hr1 = 0.69, hr2 = 0.88, 
                            id1 = 210, id2 = 420,     
                            d2min = 50, d2max = 150, stepd2 = 5, 
                            hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,     
                            adj = "all",   
                            lambdamin = 0.6, lambdamax = 1, steplambda = 0.05,  
                            alphaCImin = 0.25, alphaCImax = 0.5, stepalphaCI = 0.025,  
                            alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   
                            K = Inf, N = Inf, S = -Inf,   
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,  
                            b1 = 1000, b2 = 2000, b3 = 3000,  
                            fixed = FALSE, num_cl = 2)$u, c(98.28,75.95,96.69,77.40))
})
