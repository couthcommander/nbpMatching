library(nbpMatching)

mdm <- distancematrix('sample_dist_mat.csv')
res <- nonbimatch(mdm)

expect_true(inherits(res, 'nonbimatch'))
expect_equal(length(res), 4)
expect_equal(dim(res$matches), c(28, 5))
expect_equal(dim(res$halves), c(14, 5))
expect_equivalent(res$total, 4.16792951)
expect_equivalent(res$mean, 4.16792951 / 14)

res1 <- nonbimatch(mdm, threshold = 0.5)
expect_warning(res2 <- nonbimatch(mdm, precision = 0))
res3 <- nonbimatch(mdm, precision = 8)
expect_warning(res4 <- nonbimatch(mdm, precision = 10))
expect_equal(grep('chameleon', res1$matches[,1]), 29:34)
expect_equal(res2, res)
expect_equal(res3, res)
expect_equal(res4, res)

set.seed(1)
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
res5 <- runner(df, mate.random = TRUE)
expect_equal(length(res5), 6)
expect_equal(length(full.qom(df)), 6)

qm <- qom(res5$setup$cov, res5$matches$matches, use.se = TRUE, all.vals = TRUE)
expect_equal(length(qm), 3)
