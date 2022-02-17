RobustSIM <- function(x, y, interest=c(1:2, p-1, p),
                      outlier.prop=0, outlier.multi=10,
                      penalty="lasso", nfolds=10) {
  need.pakgs <- c("ncvreg", "glmnet")
  has <- need.pakgs %in% rownames(installed.packages())
  if (any(!has)) {
    cat("Please install and library the dependent packages:", need.pakgs, "\n")
    install.packages(need.pakgs[!has])
    library(need.pakgs[!has])
  }


  ## the empirical distribution transformation of response
  # add outliers to y
  n <- nrow(x); p <- ncol(x)
  outlier.pos <- sample(1:n, n*outlier.prop)
  y[outlier.pos] <- y[outlier.pos] + outlier.multi*max(y)
  Fy <- (rank(y))/length(y)

  ## penalized estimator
  fit_func <- function(penalty, est="theta") {
    if(est == "theta") {
      model.x <- cv.ncvreg(z, xj, penalty=penalty, nfolds=nfolds)
      lambda.x <- model.x$lambda.min
      model.x <- ncvfit(z, xj, penalty=penalty, lambda=lambda.x)
      theta.hat <- as.numeric(model.x$beta)
      return(theta.hat)
    }
    if(est == "beta") {
      model.Fy <- cv.ncvreg(x, Fy-1/2, penalty=penalty, nfolds=nfolds)
      lambda.Fy <- model.Fy$lambda.min
      model.Fy <- ncvfit(x, Fy-1/2, penalty=penalty, lambda=lambda.Fy)
      beta.hat <- as.numeric(model.Fy$beta)
      return(beta.hat)
    }
    if(est == "beta.y") {
      model.y <- cv.ncvreg(x, y, penalty=penalty, nfolds=nfolds)
      lambda.y <<- model.y$lambda.min
      model.y <- ncvfit(x, y, penalty=penalty, lambda=lambda.y)
      beta.y.hat <- as.numeric(model.y$beta)
      return(beta.y.hat)
    }
  }
  beta.hat <- fit_func(penalty, est="beta")
  pval.record <- stat.record <- c()
  for (l in 1:length(interest)) {
    cat(l, "\r")
    j <- interest[l]
    xj <<- x[,j]
    z <<- x[,-j]
    theta.hat <- fit_func(penalty, est="theta")
    gamma.hat <- beta.hat[-j]
    e.hat <- Fy-1/2-z%*%gamma.hat
    my.hat.func <- function(yi) {
      1/n * sum((xj-z%*%theta.hat)*(ifelse(y>=yi, 1, 0)-Fy))
    }
    stat.T <- 1/sqrt(n) * sum(e.hat*(xj-z%*%theta.hat))
    my.hat <- c()
    for (k in 1:n) { my.hat[k] <- my.hat.func(y[k]) }
    Sig.hat <- 1/n * sum(((xj-z%*%theta.hat)*e.hat + my.hat)^2)

    stat.S <- stat.T/sqrt(Sig.hat)
    stat.record[l] <- stat.S
    pval.record[l] <- 2*(1-pnorm(abs(stat.S)))
  }
  return(list(pvals=pval.record, teststats=stat.record))
}
