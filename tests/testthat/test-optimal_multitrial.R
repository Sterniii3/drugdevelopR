test_that("Optimal_multitrial work for Scenario 1/2", {
  expect_equal(optimal_multitrial(w = 0.3,   hr1 = 0.69, hr2 = 0.88, 
                                  id1 = 210, id2 = 420,   
                                  d2min = 20, d2max = 100, stepd2 = 5,   
                                  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,    
                                  alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,  
                                  K = Inf, N = Inf, S = -Inf,   
                                  b1 = 1000, b2 = 2000, b3 = 3000,  
                                  case = 1, strategy = 2,  
                                  fixed = FALSE,   num_cl = 2)[3], 
               data.frame(u=20.25))
})



test_that("Optimal_multitrial work for Scenario 3/4", {
  expect_equal(optimal_multitrial(w = 0.3,   hr1 = 0.69, hr2 = 0.88, 
                                  id1 = 210, id2 = 420,   
                                  d2min = 20, d2max = 100, stepd2 = 5,   
                                  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,    
                                  alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,  
                                  K = Inf, N = Inf, S = -Inf,   
                                  b1 = 1000, b2 = 2000, b3 = 3000,  
                                  case = 3, strategy = 4,  
                                  fixed = TRUE,   num_cl = 2)[3], 
               data.frame(u=-400.59))
})

test_that("Optimal_multitrial work for Scenario 2/2", {
  expect_equal(optimal_multitrial(w = 0.3,   hr1 = 0.69, hr2 = 0.88, 
                                  id1 = 210, id2 = 420,   
                                  d2min = 20, d2max = 100, stepd2 = 5,   
                                  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,    
                                  alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,  
                                  K = Inf, N = Inf, S = -Inf,   
                                  b1 = 1000, b2 = 2000, b3 = 3000,  
                                  case = 2, strategy = 2,  
                                  fixed = TRUE,   num_cl = 2)[3], 
               data.frame(u=-10.8))
})

test_that("Optimal_multitrial work for Case 3: Strategy TRUE", {
  expect_equal(optimal_multitrial(w = 0.3,   hr1 = 0.69, hr2 = 0.88, 
                                  id1 = 210, id2 = 420,   
                                  d2min = 20, d2max = 100, stepd2 = 5,   
                                  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,    
                                  alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,  
                                  K = Inf, N = Inf, S = -Inf,   
                                  b1 = 1000, b2 = 2000, b3 = 3000,  
                                  case = 3, strategy = TRUE,  
                                  fixed = TRUE,   num_cl = 2)$u, c(-182.55,-380.40,-400.59) 
               )
})

test_that("Optimal_multitrial work for Case 2: Strategy 23", {
  expect_equal(optimal_multitrial(w = 0.3,   hr1 = 0.69, hr2 = 0.88, 
                                  id1 = 210, id2 = 420,   
                                  d2min = 20, d2max = 100, stepd2 = 20,   
                                  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.1,    
                                  alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,   
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,  
                                  K = Inf, N = Inf, S = -Inf,   
                                  b1 = 1000, b2 = 2000, b3 = 3000,  
                                  case = 2, strategy = 23,  
                                  fixed = TRUE,   num_cl = 2)$u, -164.31 
  )
})



