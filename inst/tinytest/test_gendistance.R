library(nbpMatching)

# need confirmable static Mahal dist matrix

sdm <- read.csv('sample_dist_mat.csv')

set.seed(1)
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
df.dist <- gendistance(df, idcol=1, ndiscard=2)

expect_equal(length(df.dist), 9)
expect_equal(dim(df.dist$dist), c(28, 28))
expect_equal(dim(df.dist$cov), c(25, 3))
expect_equivalent(df.dist$weights, c(1, 1, 0.1))
expect_equal(df.dist$ndiscard, 3)
expect_equal(sum(abs(df.dist$dist-sdm) > 1e-8, na.rm = TRUE), 0)
expect_true(is.infinite(df.dist$dist[1,1]))

df.weighted <- gendistance(df, idcol=1, weights=c(1,2,1), ndiscard=2, missing.weight=0.25)
expect_equal(dim(df.weighted$dist), c(28, 28))
expect_equal(dim(df.weighted$cov), c(25, 3))
#expect_equivalent(df.weighted$weights, c(2, 1, 0.25))
expect_equivalent(df.weighted$weights, c(2, 1, 0.5))
expect_equal(df.weighted$ndiscard, 3)
vals <- c(2.857567, 5.416773, 1.799399, 2.57978)
expect_true(all(abs(unlist(df.weighted$dist[1,2:5])-vals) < 1e-6))
expect_true(all(abs(unlist(df.weighted$dist[2:5,1])-vals) < 1e-6))

df[,3] <- df[,2]*2
expect_warning(df.sing.solve <- gendistance(df, idcol=1, ndiscard=2))
vals <- c(0.6042711, 2.079626, 1.628953, 0.4282725)
expect_true(all(abs(unlist(df.sing.solve$dist[1,2:5])-vals) < 1e-6))
expect_true(all(abs(unlist(df.sing.solve$dist[2:5,1])-vals) < 1e-6))

expect_warning(df.sing.ginv <- gendistance(df, idcol=1, ndiscard=2, singular.method="ginv"))
vals <- c(0.4527432, 1.461176, 1.212691, 0.3519085)
expect_true(all(abs(unlist(df.sing.ginv$dist[1,2:5])-vals) < 1e-6))
expect_true(all(abs(unlist(df.sing.ginv$dist[2:5,1])-vals) < 1e-6))
