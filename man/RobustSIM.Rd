\name{RobustSIM}
\alias{RobustSIM}
\title{Robust inference for single index model}
\description{
This package provides a robust inference procedure for high-dimensional single index model with distribution transformation.
}
\usage{
RobustSIM(x, y, interest, outlier.prop, outlier.multi, penalty, nfold)
}
\arguments{
  \item{x}{the data matrix.}
  \item{y}{the response.}
  \item{interest}{coefficient of interest.}
  \item{outlier.prop}{the proportion of outliers.}
  \item{outlier.multi}{a pre-set constant standing for outlier strength.}
  \item{nfolds}{the number of folds for cross validation.}
}
\value{
  \item{pvals}{the p-values for coefficients of interest.}
  \item{teststats}{the test satstistics for coefficients of interest.}
}
\examples{
library(MASS)
library(Matrix)
n <- 200; p <- 400
p.sub <- p/10
list.temp <- NULL
for (ii in 1:10) {
  list.temp[[ii]] <- toeplitz((0.1*ii-0.1)^seq(0, p.sub-1))
}
Sig <- as.matrix(bdiag(list.temp))
x <- mvrnorm(n, mu=rep(0, p), Sigma=Sig)
error <- rnorm(n)
beta.true <- c(1, 1, 1, 1, rep(0, p-4))
y <- x\%*\%beta.true + 1*error
interest <- c(1:4, p); err.type <- 1; model <- 1;
outlier.prop <- 0; outlier.multi <- 0;
penalty <- "lasso"; nfolds <- 10;
RobustSIM(x, y, interest, outlier.prop, outlier.multi, penalty, nfolds)
}
