test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

 test_that("Epgo_normal returns the same results", {
     expect_equal(Epgo_normal(kappa=0.2,n2=50,w=0.3,
                              Delta1=0.625, Delta2=0.325, 
                              in1=300, in2=600, 
                              a=0.25,b=0.75, fixed=FALSE),0.8622991)
    })
 
 test_that("Same results for Fixed Effects if irrelevant Input Values differ",{
          expect_equal(Epgo_normal(kappa=0.2,n2=50,w=0.3,
                                  Delta1=0.625, Delta2=0.325, 
                                  in1=300, in2=600, 
                                  a=0.25,b=0.75, fixed=TRUE),
                       Epgo_normal(kappa=0.2,n2=50,w=0.7,
                                  Delta1=0.625, Delta2=0.5, 
                                  in1=350, in2=400, 
                                  a=0.1,b=0.6, fixed=TRUE))
       })