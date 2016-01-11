library(nbpMatching)
context("Generate Distance Matrix")

cov <- data.frame(
  id = paste0('subject-', seq(10)),
  gender = c(0,0,0,0,1,1,1,1,1,1),
  age = c(36, 35, 40, 40, 37, 31, 36, 23, 36, 33),
  cars = c(1,3,2,1,7,4,1,3,2,0),
  accidents = c(1,0,2,2,4,2,0,0,2,2),
  nophantom = c(1,1,1,1,0,0,0,0,0,0),
  Group2.Row = c(NA,NA,NA,NA,NA,NA,NA,NA,10,9)
)

maxval <- Inf

test_that("gendistance fails on empty data.frame", {
  testthat::skip_on_cran()
  expect_error(gendistance(cov[FALSE,]))
  expect_error(gendistance(cov[11,]))
  expect_error(gendistance(cov[11:12,]))
})

test_that("gendistance creates NxN matrix", {
  testthat::skip_on_cran()
  x <- gendistance(cov[,1:5])
  expect_equal(dim(x$dist), rep(nrow(cov), 2))
})

test_that("gendistance creates NxN matrix with phantoms", {
  testthat::skip_on_cran()
  x <- gendistance(cov[,1:5], ndiscard=2)
  expect_equal(dim(x$dist), rep(nrow(cov)+2, 2))
  expect_equal(x$ndiscard, 2)
})

test_that("gendistance takes column index or name as input", {
  testthat::skip_on_cran()
  x <- gendistance(cov[,1:5], idcol=1, prevent=2, rankcols=3:4)
  expect_equal(names(x$ignored), c('id','gender'))
  expect_equal(x$prevent, 'gender')
  expect_equal(x$rankcols, 1:2)
  y <- gendistance(cov[,1:5], idcol='id', prevent='gender', rankcols=c('age','cars'))
  expect_equal(x, y)
})

test_that("gendistance respects idcol", {
  testthat::skip_on_cran()
  x <- gendistance(cov[,1:5], idcol='id')
  expect_equal(rownames(x$cov), rownames(x$dist), colnames(x$dist), cov$id)
})

test_that("gendistance respects prevent", {
  testthat::skip_on_cran()
  x <- gendistance(cov[,1:5], idcol=1, prevent=2)
  expect_equal(x$prevent, 'gender')
  expect_true(all(x$dist[1:4,1:4] == maxval))
  expect_true(all(x$dist[5:10,5:10] == maxval))
})

test_that("gendistance respects mates", {
  testthat::skip_on_cran()
  x <- gendistance(cov[,-6], idcol='id', force='Group2.Row')
  expect_equal(x$mates, cov[,'Group2.Row'])
  expect_true(all(x$dist[9:10,1:8] == maxval))
  expect_true(all(x$dist[1:8,9:10] == maxval))
  expect_equal(x$dist[9,10], x$dist[10,9], 0)
})

test_that("gendistance respects talisman", {
  testthat::skip_on_cran()
  x <- gendistance(cov[,1:6], idcol=1, prevent=2, ndiscard=2, talisman='nophantom')
  expect_true(all(x$dist[1:4,11:12] == maxval))
  expect_true(all(x$dist[5:10,11:12] == 0))
})

test_that("gendistance imputes missing values", {
  testthat::skip_on_cran()
  cv <- cov[rep(1:10,3),2:5]
  cv[11:20,'age'] <- cv[11:20,'age'] - cv[11:20,'cars'] + cv[11:20,'accidents']
  cv[11:20,'cars'] <- cv[11:20,'cars'] - cv[11:20,'gender'] + cv[11:20,'accidents']
  cv[21:30,'age'] <- cv[21:30,'age'] - cv[21:30,'cars'] - cv[21:30,'accidents']
  cv[21:30,'cars'] <- cv[21:30,'age'] + cv[21:30,'gender']
  cv[7,'age'] <- NA
  x <- gendistance(cv[,1:3], weights=c(4,2,1), missing.weight=0.5)
  expect_equal(x$weights, c(4,2,1,1))
  expect_equal(names(x$cov)[4], 'age.missing')
  expect_equal(x$rankcols, 4)
  expect_equal(x$missing.weight, 0.5)
})

test_that("distancematrix can use gendistance", {
  testthat::skip_on_cran()
  x <- gendistance(cov, idcol='id', prevent='gender', force='Group2.Row', ndiscard=2, talisman='nophantom')
  dm <- distancematrix(x)
  expect_true(class(dm) == 'distancematrix')
})
