library(nbpMatching)

mdm <- distancematrix('sample_dist_mat.csv')
res <- nonbimatch(mdm)
expect_equal(length(get.sets(res)), 28)
expect_equal(subsetMatches(res, halvesOnly = TRUE), res$halves)
expect_equal(subsetMatches(res, halvesOnly = FALSE), res$matches)

xx <- scalar.dist(1:3)
expect_equal(xx, matrix(c(0,1,2,1,0,1,2,1,0), 3))

expect_equal(length(quantile(mdm)), 101)
