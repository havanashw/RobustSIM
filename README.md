# RobustSIM
The goal of RobustSIM is to provide robust inference procedures for high-dimensional single index model with distribution transformation.


## Installation


## Example
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
y <- x%*%beta.true + 1*error
interest <- c(1:4, p); err.type <- 1; model <- 1; 
outlier.prop <- 0; outlier.multi <- 0;
penalty <- "lasso"; nfolds <- 10;

RSIM(x, y, interest, err.type,
     model, outlier.prop, outlier.multi,
     penalty, nfolds)


