#' Iteratively Reweighted Least Squares
#'
#' This function estimates parameters using the iteratively reweighted least squares method.
#'
#' @param data a data frame
#' @param yi either a character string or the column number of the outcome variable in the data frame.
#' @param xi a vector of either character strings or the column numbers of the covariates.
#' @param niter numeric for number of iterations (optional argument)
#'
#' @return a matrix of the parameter estimates (includes intercept), standard errors, 95% confidence intervals, t statistics, and p-values
#'
#' @examples
#'
#' est <- irls(data=data, yi="y", xi=c("x1","x2","x3"))
#' 
#' est <- irls(data=data, yi=1, xi=c(1:3))
#'
#' @export

irls = function(data, yi, xi, niter=200) {
  data <- data[complete.cases(data),]               # remove missing values
  n <- length(data[,1])                             # total length of data set
  x <- cbind(replicate(n,1),as.matrix(data[,xi]))   # combine with column of 1s (intercept)
  m <- length(x[1,])                                # number of parameters
  y <- data[,yi]                                    # outcome variable
  beta <- rnorm(m)                                  # generate initial beta
  
  logit = function(x) {
    return(exp(x) / (1 + exp(x)))
  }
  
  iter = 0
  maxit = niter
  eps = Inf
  tol = 10^-4
  
  while (eps > tol & iter < maxit) {
    
    beta0 = beta
    p <- numeric(0)
    mu <- sapply(x %*% beta, logit)
    for (i in 1:n)
      p[i] <- max(mu[i] * (1 - mu[i]),0.005)
    W = diag(p)
    z = x %*% beta + solve(W) %*% (y - mu)
    
    beta = solve(t(x) %*% W %*% x) %*% t(x) %*% W %*% z
    
    eps = sqrt(sum((beta - beta0)^2))
    
    iter = iter + 1
    if (iter == maxit)
      warning("Iteration limit reached without convergence")
    # print out info to keep track
    cat(sprintf(
        "Iter: %d beta0: %.3f beta1: %.3f beta2: %.3f eps:%f\n",
        iter,
        beta[1],
        beta[2],
        beta[3],
        eps)
    )
    
  }
  var = diag(solve(t(x) %*% W %*% x))
  se = sqrt(var)
  lower95ci = beta - qt(1-0.025, df=n)*sqrt(var)
  upper95ci = beta + qt(1-0.025, df=n)*sqrt(var)
  t = beta/sqrt(var)
  Pvalue = 2*pt(-abs(t),df=n)
  tab <- cbind(beta,se,lower95ci,upper95ci,t,Pvalue)
  colnames(tab) <- c("Parameters","Standard Error","Lower 95% CI","Upper 95% CI","t Statistic","P-value")
  return(tab)
}
