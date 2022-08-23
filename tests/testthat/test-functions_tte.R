test_that("prior_tte returns the same results", {
  expect_equal(prior_tte(x=0.5, w=0.3, 
                         hr1=0.9, hr2=0.7, 
                         id1=200, id2=400), 1.01711002)
})

test_that("parameter value for w irrelevant if treatment effect and amount of information are the same", {
  expect_equal(prior_tte(x=0.5, w=0.3, hr1=0.7, hr2=0.7, id1=400, id2=400), 
               prior_tte(x=0.5, w=0.5, hr1=0.7, hr2=0.7, id1=400, id2=400))
})

test_that("box_tte returns a vector with 100000 realizations", {
  expect_length(box_tte(w=0.3,hr1=1.5,hr2=0.95,
                           id1=200,id2=400), 1000000)
})

test_that("Higher hazard ratio leads to lower probabiltiy to go to phase III", {
   expect_lt(Epgo_tte(HRgo = 0.5, d2 = 50, w=NULL, hr1=0.6, hr2=NULL,
                      id1=NULL, id2=NULL, fixed=TRUE),
             Epgo_tte(HRgo = 0.5, d2 = 50, w=NULL, hr1=0.4, hr2=NULL,
                      id1=NULL, id2=NULL, fixed=TRUE))       
  })

test_that("Hazard ratio equal to go decision leads to 50% chance of going to phase III", {
  expect_lt(Epgo_tte(HRgo = 0.4, d2 = 50, w=NULL, hr1=0.6, hr2=NULL,
                     id1=NULL, id2=NULL, fixed=TRUE),
            Epgo_tte(HRgo = 0.4, d2 = 50, w=NULL, hr1=0.4, hr2=NULL,
                     id1=NULL, id2=NULL, fixed=TRUE))       
})

test_that("Higer value for second hazard ratio increases number of expected events in phase III",{
  expect_gt(Ed3_tte(HRgo = 0.5, d2 = 50, alpha = 0.05, beta = 0.1, 
                    w = 0.5, hr1 = 0.3, hr2 = 0.5, id1 = 200, id2 = 400, fixed = FALSE),
            Ed3_tte(HRgo = 0.5, d2 = 50, alpha = 0.05, beta = 0.1, 
                    w = 0.5, hr1 = 0.3, hr2 = 0.7, id1 = 200, id2 = 400, fixed = FALSE))
})
