library(nbpMatching)
context("Run Matching")

cov <- data.frame(
  id = paste0('subject-', seq(10)),
  gender = c(0,0,0,0,1,1,1,1,1,1),
  age = c(36, 35, 40, 40, 37, 31, 36, 23, 36, 33),
  cars = c(1,3,2,1,7,4,1,3,2,0),
  accidents = c(1,0,2,2,4,2,0,0,2,2),
  nophantom = c(1,1,1,1,0,0,0,0,0,0),
  Group2.Row = c(NA,NA,NA,NA,NA,NA,NA,NA,10,9)
)

test_that("nonbimatch creates matches", {
  x <- gendistance(cov[,-6], idcol='id', force='Group2.Row')
  dm <- distancematrix(x)
  mt <- nonbimatch(dm, precision=9)
  expect_equal(mt$matches[,4], c(8,7,4,3,6,5,2,1,10,9))
  x <- gendistance(cov, idcol='id', force='Group2.Row', ndiscard=2, talisman='nophantom')
  dm <- distancematrix(x)
  mt <- nonbimatch(dm, precision=9)
  expect_equal(mt$matches[,4], c(2,1,4,3,12,8,11,6,10,9,7,5))
})
