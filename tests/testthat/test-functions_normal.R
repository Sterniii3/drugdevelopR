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