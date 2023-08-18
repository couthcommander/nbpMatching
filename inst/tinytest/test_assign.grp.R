library(nbpMatching)

mdm <- distancematrix('sample_dist_mat.csv')
res <- nonbimatch(mdm)
g1 <- assign.grp(res, seed = 10)
g2 <- assign.grp(res$matches, seed = 10)

expect_equal(dim(g1), c(28, 6))
expect_equal(g1, g2)
expect_equal(g1[1:10,'treatment.grp'], c('A','A','B','B','B','A','B','B','A','A'))
