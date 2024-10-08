% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendistance.R
\name{gendistance}
\alias{gendistance}
\alias{gendistance,data.frame-method}
\title{Generate a Distance Matrix}
\usage{
gendistance(
  covariate,
  idcol = NULL,
  weights = NULL,
  prevent = NULL,
  force = NULL,
  rankcols = NULL,
  missing.weight = 0.1,
  ndiscard = 0,
  singular.method = "solve",
  talisman = NULL,
  prevent.res.match = NULL,
  outRawDist = FALSE,
  ...
)
}
\arguments{
\item{covariate}{A data.frame object, containing the covariates of the data
set.}

\item{idcol}{An integer or column name, providing the index of the column
containing row ID's.}

\item{weights}{A numeric vector, the length should match the number of
columns.  This value determines how much weight is given to each column when
generating the distance matrix.}

\item{prevent}{A vector of integers or column names, providing the index of
columns that should be used to prevent matches.  When generating the distance
matrix, elements that match on these columns are given a maximum distance.}

\item{force}{An integer or column name, providing the index of the column
containing information used to force pairs to match.}

\item{rankcols}{A vector of integers or column names, providing the index of
columns that should have the rank function applied to them before generating
the distance matrix.}

\item{missing.weight}{A numeric value, or vector, used to generate the weight
of missingness indicator columns.  Missingness indicator columns are created
if there is missing data within the data set.  Defaults to 0.1.  If a single
value is supplied, weights are generating by multiplying this by the original
columns' weight.  If a vector is supplied, it's length should match the
number of columns with missing data, and the weight is used as is.}

\item{ndiscard}{An integer, providing the number of elements that should be
allowed to match phantom values.  The default value is 0.}

\item{singular.method}{A character string, indicating the function to use
when encountering a singular matrix.  By default, \code{\link{solve}} is called.
The alternative is to call \code{\link[MASS]{ginv}} from the \pkg{MASS} package.}

\item{talisman}{An integer or column name, providing location of talisman column.
The talisman column should only contains values of 0 and 1.  Records with zero
will match phantoms perfectly, while other records will match phantoms at max distance.}

\item{prevent.res.match}{An integer or column name, providing location of the column
containing assigned treatment groups.  This is useful in some settings, such as
trickle-in randomized trials.  When set, non-NA values from this column are
replaced with the value 1.  This prevents records with previously assigned
treatments (the \sQuote{reservior}) from matching each other.}

\item{outRawDist}{a logical, indicating if the raw distance matrix should also
be returned.  The raw form is before distance modifiers such as \sQuote{prevent}
take effect.}

\item{\dots}{Additional arguments, not used at this time.}
}
\value{
a list object with several elements

 \item{dist}{generated distance matrix}

 \item{cov}{covariate matrix used to generate distances}

 \item{ignored}{ignored columns from original covariate matrix}

 \item{weights}{weights applied to each column in covariate matrix}

 \item{prevent}{columns used to prevent matches}

 \item{mates}{index of rows that should be forced to match}

 \item{rankcols}{index of columns that should use rank}

 \item{missing.weight}{weight to apply to missingness indicator columns}

 \item{ndiscard}{number of elements that will match phantoms}

 \item{rawDist}{raw distance matrix, only provided if \sQuote{outRawDist} is TRUE}
}
\description{
The gendistance function creates an \eqn{(N+K)}x\eqn{(N+K)} distance matrix
from an \eqn{N}x\eqn{P} covariates matrix, where \eqn{N} is the number
of subjects, \eqn{P} the number of covariates, and \eqn{K} the number of
phantom subjects requested (see \code{ndiscard} option). Provided the
covariates' covariance matrix is invertible, the distances computed are
Mahalanobis distances, or if covariate weights are provided, Reweighted
Mahalanobis distances (see \code{weights} option and Greevy, et al.,
Pharmacoepidemiology and Drug Safety 2012).
}
\details{
Given a data.frame of covariates, generate a distance matrix.  Missing values
are imputed with \code{\link{fill.missing}}.  For each column with missing
data, a missingness indicator column will be added.  Phantoms are fake
elements that perfectly match all elements.  They can be used to discard a
certain number of elements.
}
\examples{

set.seed(1)
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
# add some missing data
df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
df.dist <- gendistance(df, idcol=1, ndiscard=2)
# up-weight the second column
df.weighted <- gendistance(df, idcol=1, weights=c(1,2,1), ndiscard=2, missing.weight=0.25)
df[,3] <- df[,2]*2
df.sing.solve <- gendistance(df, idcol=1, ndiscard=2)
df.sing.ginv <- gendistance(df, idcol=1, ndiscard=2, singular.method="ginv")

}
\seealso{
\code{\link{distancematrix}}
}
\author{
Cole Beck
}
