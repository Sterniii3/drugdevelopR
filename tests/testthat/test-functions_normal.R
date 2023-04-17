test_that("prior_normal returns the same results", {
  expect_equal(prior_normal(x=0.5,w=0.3,Delta1=1.5,Delta2=0.95,
                            in1=200,in2=400,a=0.4,b=0.6),
               0.81469965)
})

test_that("expanding boundaries reduces value", {
  expect_gt(prior_normal(x=0.5,w=0.3,Delta1=1.5,Delta2=0.95,
                            in1=200,in2=400,a=0.4,b=0.6),
               prior_normal(x=0.5,w=0.3,Delta1=1.5,Delta2=0.95,
                            in1=200,in2=400,a=0.2,b=0.8))
})

test_that("box normal returns a vector with 100000 realizations", {
  expect_length(box_normal(w=0.3,Delta1=1.5,Delta2=0.95,
                           in1=200,in2=400,a=0.4,b=0.6), 1000000)
})



 test_that("Epgo_normal returns the same results", {
     expect_equal(Epgo_normal(kappa=0.2,n2=50,w=0.3,
                              Delta1=0.625, Delta2=0.325, 
                              in1=300, in2=600, 
                              a=0.25,b=0.75, fixed=FALSE),0.8622991)
    })
 
 test_that("Same result for fixed effects if irrelevant input values differ",{
          expect_equal(Epgo_normal(kappa=0.2,n2=50,w=0.3,
                                  Delta1=0.625, Delta2=0.325, 
                                  in1=300, in2=600, 
                                  a=0.25,b=0.75, fixed=TRUE),
                       Epgo_normal(kappa=0.2,n2=50,w=0.7,
                                  Delta1=0.625, Delta2=0.5, 
                                  in1=350, in2=400, 
                                  a=0.1,b=0.6, fixed=TRUE))
       })
 
 test_that("Different Calculation of Epgo",{
           expect_equal(Epgo_normal(kappa=0.2,n2=50,w=0.3,
                                     Delta1=0.625, Delta2=0.325, 
                                     in1=300, in2=600, 
                                     a=0.25,b=0.75, fixed=TRUE),
                         pnorm((0.625-0.2)/sqrt(4/50),0,1)           )
 })
 
 
 test_that("Higher Treatment effect reduces n3",{
            expect_gt(En3_normal(kappa=0.2,n2=50,w=0.3,
                                 alpha = 0.05, beta = 0.1,
                                 Delta1=0.625, Delta2=0.325, 
                                 in1=300, in2=600, 
                                 a=0.25,b=0.75, fixed=TRUE),
                      En3_normal(kappa=0.2,n2=50,w=0.3,
                                 alpha = 0.05, beta = 0.1,
                                 Delta1=0.825, Delta2=0.325, 
                                 in1=300, in2=600, 
                                 a=0.25,b=0.75, fixed=TRUE))
 })
 
 
 test_that("Higher Treatment effect increases probability of a success",{
     expect_lt(EPsProg_normal(kappa=0.2,n2=50,w=0.3,
                          alpha = 0.05, beta = 0.1,
                          step1 = 0.5, step2 = 1,
                          Delta1=0.625, Delta2=0.325, 
                          in1=300, in2=600, 
                          a=0.25,b=0.75, fixed=TRUE, gamma=0),
               EPsProg_normal(kappa=0.2,n2=50,w=0.3,
                          alpha = 0.05, beta = 0.1,
                          step1 = 0.5, step2 = 1,
                          Delta1=0.825, Delta2=0.325, 
                          in1=300, in2=600, 
                          a=0.25,b=0.75, fixed=TRUE, gamma=0))
 })
 
 test_that("Bigger success region increases probability of a success",{
     expect_lt(EPsProg_normal(kappa=0.2,n2=50,w=0.3,
                              alpha = 0.05, beta = 0.1,
                              step1 = 0.5, step2 = 1,
                              Delta1=0.625, Delta2=0.325, 
                              in1=300, in2=600, 
                              a=0.25,b=0.75, fixed=TRUE, gamma=0),
               EPsProg_normal(kappa=0.2,n2=50,w=0.3,
                              alpha = 0.05, beta = 0.1,
                              step1 = 0.3, step2 = 1,
                              Delta1=0.825, Delta2=0.325, 
                              in1=300, in2=600, 
                              a=0.25,b=0.75, fixed=TRUE, gamma=0))
 })
 
 test_that("Higher gains for each benefit category increase utility",{
    expect_lt(getElement(utility_normal(n2=60, kappa=0.1, w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                             in1=300, in2=600, a = 0.25, b = 0.75,
                             alpha = 0.05, beta = 0.1, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                             K=Inf, N=Inf, S=-Inf,
                             steps1=0, stepm1=0.5, stepl1=0.8,
                             b1=1000, b2=2000, b3=3000,
                             gamma=0, fixed=TRUE),1),
              getElement(utility_normal(n2=60, kappa=0.1, w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                             in1=300, in2=600, a = 0.25, b = 0.75,
                             alpha = 0.05, beta = 0.1, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                             K=Inf, N=Inf, S=-Inf,
                             steps1=0, stepm1=0.5, stepl1=0.8,
                             b1=2000, b2=4000, b3=6000,
                             gamma=0, fixed=TRUE),1))
 })
 
 
 test_that("Sample size constraint works",{
   expect_equal(getElement(utility_normal(n2=100, kappa=0.1, w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                                       in1=300, in2=600, a = 0.25, b = 0.75,
                                       alpha = 0.05, beta = 0.1, 
                                       c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                                       K=Inf, N=90, S=-Inf,
                                       steps1=0, stepm1=0.5, stepl1=0.8,
                                       b1=1000, b2=2000, b3=3000,
                                       gamma=0, fixed=TRUE),1),-9999)
 })
 
 test_that("Cost constraint works",{
   expect_equal(getElement(utility_normal(n2=60, kappa=0.1, w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                                          in1=300, in2=600, a = 0.25, b = 0.75,
                                          alpha = 0.05, beta = 0.1, 
                                          c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                                          K=50, N=Inf, S=-Inf,
                                          steps1=0, stepm1=0.5, stepl1=0.8,
                                          b1=1000, b2=2000, b3=3000,
                                          gamma=0, fixed=TRUE),1),-9999)
 })
 
 test_that("Lower treatment effects decrease utility",{
   expect_gt(getElement(utility_normal(n2=60, kappa=0.1, w=0.3, Delta1 = 0.675, Delta2 = 0.825, 
                                       in1=300, in2=600, a = 0.25, b = 0.75,
                                       alpha = 0.05, beta = 0.1, 
                                       c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                                       K=Inf, N=Inf, S=-Inf,
                                       steps1=0, stepm1=0.5, stepl1=0.8,
                                       b1=1000, b2=2000, b3=3000,
                                       gamma=0, fixed=FALSE),1),
             getElement(utility_normal(n2=60, kappa=0.1, w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                                       in1=300, in2=600, a = 0.25, b = 0.75,
                                       alpha = 0.05, beta = 0.1, 
                                       c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                                       K=Inf, N=Inf, S=-Inf,
                                       steps1=0, stepm1=0.5, stepl1=0.8,
                                       b1=1000, b2=2000, b3=3000,
                                       gamma=0, fixed=FALSE),1))
 })
 
 test_that("Higher costs decrease utility",{
   expect_gt(getElement(utility_normal(n2=60, kappa=0.1, w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                                       in1=300, in2=600, a = 0.25, b = 0.75,
                                       alpha = 0.05, beta = 0.1, 
                                       c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20, 
                                       K=Inf, N=Inf, S=-Inf,
                                       steps1=0, stepm1=0.5, stepl1=0.8,
                                       b1=1000, b2=2000, b3=3000,
                                       gamma=0, fixed=FALSE),1),
             getElement(utility_normal(n2=60, kappa=0.1, w=0.3, Delta1 = 0.375, Delta2 = 0.625, 
                                       in1=300, in2=600, a = 0.25, b = 0.75,
                                       alpha = 0.05, beta = 0.1, 
                                       c2 = 0.8, c3 = 0.75, c02 = 20, c03 = 25, 
                                       K=Inf, N=Inf, S=-Inf,
                                       steps1=0, stepm1=0.5, stepl1=0.8,
                                       b1=1000, b2=2000, b3=3000,
                                       gamma=0, fixed=FALSE),1))
 })
 
 test_that("n3_skipII_normal returns same results",{
   expect_equal(n3_skipII_normal(alpha=0.05,beta=0.1,median_prior=0.675),
                75.183296)
 })
 
 test_that("Higher treatment effect (i.e. median prior larger) leads to smaller number of n3",{
   expect_gt(n3_skipII_normal(alpha=0.05,beta=0.1,median_prior=0.675),
             n3_skipII_normal(alpha=0.05,beta=0.1,median_prior=0.875))
 })
 
 test_that("Higher Treatment effect increases probability of a success",{
   expect_lt(EPsProg_skipII_normal(alpha = 0.05, beta = 0.1,
                            median_prior = 0.625, step1 = 0.5, step2 = 1,
                            Delta1=0.625, Delta2=0.325, 
                            in1=300, in2=600, 
                            a=0.25,b=0.75, fixed=TRUE, gamma=0),
             EPsProg_skipII_normal(alpha = 0.05, beta = 0.1,
                            median_prior = 0.625,step1 = 0.5, step2 = 1,
                            Delta1=0.825, Delta2=0.325, 
                            in1=300, in2=600, 
                            a=0.25,b=0.75, fixed=TRUE, gamma=0))
 }) 
 
 test_that("Higher Treatment effect increases probability of a success with prior distribution",{
   expect_lt(EPsProg_skipII_normal(alpha = 0.05, beta = 0.1, w=0.5,
                                   median_prior = 0.625, step1 = 0.5, step2 = 1,
                                   Delta1=0.625, Delta2=0.325, 
                                   in1=300, in2=600, 
                                   a=0.25,b=0.75, fixed=FALSE, gamma=0),
             EPsProg_skipII_normal(alpha = 0.05, beta = 0.1, w=0.5,
                                   median_prior = 0.625,step1 = 0.5, step2 = 1,
                                   Delta1=0.825, Delta2=0.325, 
                                   in1=300, in2=600, 
                                   a=0.25,b=0.75, fixed=FALSE, gamma=0))
 })
 
 test_that("Utility decreases with higher costs",{
   expect_gt(utility_skipII_normal(alpha = 0.05, beta = 0.1, median_prior = 0.625,
                                   c03 = 0.72, c3=20, b1=1000, b2=2000, b3=3000,  
                                   K=Inf, N=Inf, S=-Inf,
                                   steps1=0, stepm1=0.5, stepl1=0.8,
                                   w=0.3, Delta1=0.625, Delta2=0.325,
                                   in1=300, in2=600,a=0.25,b=0.75, 
                                   fixed=TRUE, gamma=0)[1],
             utility_skipII_normal(alpha = 0.05, beta = 0.1, median_prior = 0.625,
                                   c03 = 0.75, c3=25, b1=1000, b2=2000, b3=3000,  
                                   K=Inf, N=Inf, S=-Inf,
                                   steps1=0, stepm1=0.5, stepl1=0.8,
                                   w=0.3, Delta1=0.625, Delta2=0.325,
                                   in1=300, in2=600,a=0.25,b=0.75, 
                                   fixed=TRUE, gamma=0)[1])
 })
 
 test_that("Utility increases with higher gains",{
   expect_lt(utility_skipII_normal(alpha = 0.05, beta = 0.1, median_prior = 0.625,
                                   c03 = 0.72, c3=20, b1=1000, b2=2000, b3=3000,  
                                   K=Inf, N=Inf, S=-Inf,
                                   steps1=0, stepm1=0.5, stepl1=0.8,
                                   w=0.3, Delta1=0.625, Delta2=0.325,
                                   in1=300, in2=600,a=0.25,b=0.75, 
                                   fixed=TRUE, gamma=0)[1],
             utility_skipII_normal(alpha = 0.05, beta = 0.1, median_prior = 0.625,
                                   c03 = 0.72, c3=20, b1=2000, b2=4000, b3=6000,  
                                   K=Inf, N=Inf, S=-Inf,
                                   steps1=0, stepm1=0.5, stepl1=0.8,
                                   w=0.3, Delta1=0.625, Delta2=0.325,
                                   in1=300, in2=600,a=0.25,b=0.75, 
                                   fixed=TRUE, gamma=0)[1])
 })
 
 
 test_that("Violating constraints leads to negative utility of -9999",{
   expect_equal(utility_skipII_normal(alpha = 0.05, beta = 0.1, median_prior = 0.625,
                                   c03 = 0.72, c3=20, b1=1000, b2=2000, b3=3000,  
                                   K=50, N=50, S=0.8,
                                   steps1=0, stepm1=0.5, stepl1=0.8,
                                   w=0.3, Delta1=0.625, Delta2=0.325,
                                   in1=300, in2=600,a=0.25,b=0.75, 
                                   fixed=TRUE, gamma=0)[1],
             -9999)
 })
