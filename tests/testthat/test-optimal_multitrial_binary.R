test_that("optimal_multitrial_binary works for case 1", {
  skip_on_cran()
  expect_equal(optimal_multitrial_binary(w = 0.3, p0 = 0.6, p11 =  0.3, p12 = 0.5,
                                         in1 = 30, in2 = 60,
                                         n2min = 20, n2max = 100, stepn2 = 4,
                                         rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,
                                         alpha = 0.05, beta = 0.1,
                                         c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                         K = Inf, N = Inf, S = -Inf,   
                                         b1 = 1000, b2 = 2000, b3 = 3000,
                                         case = 1, strategy = TRUE, 
                                         fixed = TRUE,num_cl = 2)$u, c(1806.86,1900.29))
})

test_that("optimal_multitrial_binary works for case 2", {
  skip_on_cran()
  expect_equal(optimal_multitrial_binary(w = 0.3, p0 = 0.6, p11 =  0.3, p12 = 0.5,
                                         in1 = 30, in2 = 60,
                                         n2min = 20, n2max = 100, stepn2 = 4,
                                         rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,
                                         alpha = 0.05, beta = 0.1,
                                         c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                         K = Inf, N = Inf, S = -Inf,   
                                         b1 = 1000, b2 = 2000, b3 = 3000,
                                         case = 2, strategy = TRUE, 
                                         fixed = FALSE,num_cl = 2)$u, c(325.53,801.12,1034.92,822.46))
})

test_that("optimal_multitrial_binary works for case 3", {
  skip_on_cran()
  expect_equal(optimal_multitrial_binary(w = 0.3, p0 = 0.6, p11 =  0.3, p12 = 0.5,
                                         in1 = 30, in2 = 60,
                                         n2min = 150, n2max = 250, stepn2 = 4,
                                         rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,
                                         alpha = 0.05, beta = 0.1,
                                         c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                         K = Inf, N = Inf, S = -Inf,   
                                         b1 = 1000, b2 = 2000, b3 = 3000,
                                         case = 3, strategy = TRUE, 
                                         fixed = TRUE,num_cl = 2)$u, c(1209.72,819.66,981.51))
})

