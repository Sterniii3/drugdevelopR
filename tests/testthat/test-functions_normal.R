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