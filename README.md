# RobustSIM
R code for "Robust group and simultaneous inferences for high-dimensional single index model".


## Installation
```{R}
install.packages("remotes")
remotes::install_github("havanashw/RobustSIM")
```
```{R}
install.packages("devtools")
devtools::install_github("havanashw/RobustSIM")
```
## Example of statistical inference based on RobustSIM
```{R}
library(RobustSIM)
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
y <- x%*%beta.true + 1*error
interest <- c(1:4, p); err.type <- 1; model <- 1;
outlier.prop <- 0; outlier.multi <- 0;
penalty <- "lasso"; nfolds <- 10;
RobustSIM(x, y, interest, outlier.prop, outlier.multi, penalty, nfolds)
```

## Example of fdr control based on RobustSIM
```{R}
library(MASS)
library(Matrix)
library(ncvreg)
library(RobustSIM)
library(snowfall)

# the snowfall package is used for parallel computing on a multi-CPU computer 
sfInit(parallel=TRUE, cpus=4, type="SOCK")
sfLibrary(MASS)
sfLibrary(Matrix)
sfLibrary(ncvreg)
sfLibrary(RobustSIM)

FDR_func <- function(n=200, p=400, q=4, err.type=1, 
                     outlier.prop=0, outlier.multi=0,
                     model=1, penalty="lasso", nfolds=10,
                     alpha.fdr=0.1, i=1) {
  # n: the data size
  # p: the data dimension
  # q: sparsity level
  # err.type: error distribution
  #           1: standard normal distribution
  #           2: cauchy distribution
  # outlier.prop: the proportion of outliers
  # outlier.multi: a pre-set constant standing for outlier strength
  # model: 1: linear model
  #        2: generalized linear model
  #        3: non-linear model
  # penalty: the choice of penalty
  # nfolds: number of folds for cross validation
  # alpha.fdr: the fdr level
  
  beta.true <- c(rep(1, q), rep(0, p-q))
  act.set <- (1:p)[beta.true != 0]
  nonact.set <- (1:p)[-act.set]
  
  ## generate data
  # x
  genedata <- function(i) {
    y <- Inf
    while (length(which(y == Inf)) != 0) {
      p.sub <- p/10
      list.temp <- NULL
      for (ii in 1:10) {
        list.temp[[ii]] <- toeplitz((0.1*ii-0.1)^seq(0, p.sub-1))
      }
      Sig <- as.matrix(bdiag(list.temp))
      x <- mvrnorm(n, mu=rep(0, p), Sigma=Sig)
      
      # error
      if(err.type == 1) { error <- rnorm(n) }
      if(err.type == 2) { error <- rt(n, 1) }
      
      # linear model
      if(model == 1) { y <- x%*%beta.true + error } 
      # generalized model
      if(model == 2) { y <- rpois(n, exp(x%*%beta.true)) }
      # non-linear model
      if(model == 3) { 
        y <- sin(0.5*x%*%beta.true) * exp(x%*%beta.true + 1) + error
      }
      
      ## add outliers to y
      outlier.pos <- sample(1:n, n*outlier.prop)
      y[outlier.pos] <- y[outlier.pos] + outlier.multi*max(y)
      Fy <- (rank(y))/length(y)
    }
    return(list(x=x, y=y, Fy=Fy))
  }
  
  # set.seed(i)
  outdata <- genedata(i)
  x <- outdata$x
  y <- outdata$y
  
  sfExport("x", "y", "outlier.prop", "outlier.multi", "penalty", "nfolds")
  sim.temp <- sfLapply(1:p, function(l) {
    out.res <- RobustSIM(x, y, l, outlier.prop, outlier.multi, penalty, nfolds)
    stat.record <- out.res$teststats
    return(stat.record)
  })
  stat.res <- unlist(sim.temp)
  
  ## SIM fdr control
  FDP.hat <- function(t, stat.res) {
    num <- 2*p*(1-pnorm(t))
    den <- max(length((1:p)[(abs(stat.res) >= t)]), 1)
    return(num/den)
  }
  
  FDP <- function(t, stat.res) {
    num <- length(nonact.set[(abs(stat.res[nonact.set]) >= t)])
    den <- max(length((1:p)[(abs(stat.res) >= t)]), 1)
    return(num/den)
  }
  
  bp <- sqrt(2*log(p))
  t <- seq(0, sqrt(2*log(p) - (log(log(p)))), 0.001)
  t.vec <- c(); k <- 1
  for (ti in t) {
    if (FDP.hat(ti, stat.res) <= alpha.fdr) {
      t.vec[k] <-  ti
      k <- k + 1
    }
  }
  t.hat <- ifelse(length(t.vec) != 0, min(t.vec), bp)
  fdp.SIM <- FDP(t.hat, stat.res)
  choose.act <- intersect((1:p)[abs(stat.res) >= t.hat], act.set)
  power.SIM <-  length(choose.act)/length(act.set)
  return(c(fdp.SIM, power.SIM))
}

## example
n=100; p=30; q=4; err.type=1;
outlier.prop=0; outlier.multi=0;
model=1; penalty="lasso"; nfolds=10;
alpha.fdr=0.2; i=1

FDR_func(n, p, q, err.type, outlier.prop, outlier.multi,
         model, penalty, nfolds, alpha.fdr, i)
```

