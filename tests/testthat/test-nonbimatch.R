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
    testthat::skip_on_cran()
    x <- gendistance(cov[,-6], idcol='id', force='Group2.Row')
    dm <- distancematrix(x)
    mt <- nonbimatch(dm, precision=9)
    expect_equal(mt$matches[,4], c(8,7,4,3,6,5,2,1,10,9))
    x <- gendistance(cov, idcol='id', force='Group2.Row', ndiscard=2, talisman='nophantom')
    dm <- distancematrix(x)
    mt <- nonbimatch(dm, precision=9)
    expect_equal(mt$matches[,4], c(2,1,4,3,12,8,11,6,10,9,7,5))
})

test_that("nonbimatch handles zeroed distancematrix", {
    testthat::skip_on_cran()
    dm <- distancematrix(matrix(0, nrow=6, ncol=6))
    expect_silent(nonbimatch(dm))
})

test_that("nonbimatch receives precision warnings", {
    testthat::skip_on_cran()
    dm <- distancematrix(matrix(1e10, nrow=16, ncol=16))
    expect_warning(nonbimatch(dm, precision=12))
})

test_that("nonbimatch object is valid", {
    testthat::skip_on_cran()
    x <- gendistance(cov[,-6], idcol='id', force='Group2.Row')
    dm <- distancematrix(x)
    mt <- nonbimatch(dm, precision=8)
    expect_s4_class(mt, 'nonbimatch')
    expect_true(is.list(mt))
    expect_equal(names(mt), c('matches','halves','total','mean'))
})

test_that("nonbimatch can be pruned with subsetMatches", {
    testthat::skip_on_cran()
    dat <- cbind(cov[,1:5], block=c(rep(1,6),rep(0,4)))
    x <- gendistance(dat, idcol='id', prevent='block')
    dm <- distancematrix(x)
    mt <- nonbimatch(dm)
    expect_equal(nrow(subsetMatches(mt)), sum(!is.infinite(mt$halves[,'Distance'])))
    expect_equal(subsetMatches(mt, infinite=FALSE), mt$halves)
    x <- gendistance(dat, idcol='id', ndiscard=2, prevent='block')
    dm <- distancematrix(x)
    mt <- nonbimatch(dm)
    expect_equal(nrow(subsetMatches(mt)), sum(!grepl('phantom', mt$halves[,'Group2.ID'])))
    expect_equal(subsetMatches(mt, phantom=FALSE), mt$halves)
    expect_warning(dm <- distancematrix(x$dist[-1,-1]))
    mt <- nonbimatch(dm, threshold=2)
    nr <- nrow(mt$halves)
    nphant <- sum(grepl('phantom', mt$halves[,'Group2.ID']))
    nchame <- sum(grepl('chameleon', mt$halves[,'Group2.ID']))
    nghost <- sum(grepl('ghost', mt$halves[,'Group2.ID']))
    expect_equal(subsetMatches(mt, phantom=FALSE, chameleon=FALSE, ghost=FALSE), mt$halves)
    expect_equal(nrow(subsetMatches(mt)), nr-nphant-nchame-nghost)
    expect_equal(nrow(subsetMatches(mt, phantom=FALSE, chameleon=FALSE)), nr-nghost)
    expect_equal(nrow(subsetMatches(mt, ghost=FALSE, phantom=FALSE)), nr-nchame)
    expect_equal(nrow(subsetMatches(mt, ghost=FALSE, chameleon=FALSE)), nr-nphant)
    expect_equal(nrow(subsetMatches(mt, ghost=FALSE)), nr-nphant-nchame)
    expect_equal(nrow(subsetMatches(mt, phantom=FALSE)), nr-nchame-nghost)
    expect_equal(nrow(subsetMatches(mt, chameleon=FALSE)), nr-nphant-nghost)
})
