#' Calculate linear mixed effects coefficients and covariance matrix (identity link)
#' Code can be found in: http://sia.webpopix.org/EMlme.html
#'
#' This code calculates the coefficients of the fixed effects and
#' covariance matrix of the random effects specified in the model.
#'
#' @param y a vector of the outcome variable
#' @param id a vector of the individual level id variable
#' @param fixed a matrix of fixed effects (intercept term optional)
#' @param mixed a matrix of random effects (intercept term optional)
#' @param niter numeric for number of EM iterations
#'
#' @return a list of fixed effects coefficients, random effects cov. matrix,
#' and model selection metrics (logLik, AIC, BIC).
#'
#' @examples
#'
#' X <- cbind(seq(1,100,1), seq(1,200,2), seq(1,500,5))
#' glmm(X[,3], X[,1], X[,2], X[,2])
#'
#' @export
glmm <- function(y,id,fixed,mixed=fixed,niter=50) {
  uid <- unique(id)
  if(sum(sapply(fixed, is.factor)) > 0 |sum(sapply(fixed, is.character)) > 0){
    warning("fixed or random effects matrix contains non-numeric variable(s)")
  }
  X <- as.matrix(fixed)
  A <- as.matrix(mixed)
  y <- as.matrix(y)
  N <- length(uid)
  n <- length(y)
  nb.eta <- ncol(A)

  beta <- as.vector(solve(t(X)%*%X)%*%t(X)%*%y)
  Omega <- diag(rep(1,nb.eta))
  sigma2 <- 1
  z <- as.vector(y - X%*%beta)
  for (k in 1:niter) {
    iO <- solve(Omega)
    T <- R <- C <- 0
    mu <- u <- NULL
    for (i in uid ) {
      row.i <- which(id==i)
      Xi <- X[row.i,]
      Ai <- A[row.i,]
      AAi <- t(Ai)%*%Ai
      zi <- z[row.i]
      Gammai <- solve(AAi/sigma2 + iO)
      mui <- (Gammai%*%t(Ai)%*%zi)/sigma2
      mu <- c(mu, mui)
      u <- c(u, Ai%*%mui)
      Si <- Gammai + mui%*%t(mui)
      R <- R + Si
      T <- T + sum(diag(Si%*%AAi))
      C <- C + t(mui)%*%t(Ai)%*%zi
    }
    beta <- as.vector(solve(t(X)%*%X)%*%t(X)%*%(y-u))
    z <- as.vector(y - X%*%beta)
    sigma2 <- (sum(z^2) -2*C[1] + T)/n
    Omega <- as.matrix(R/N)
  }
  z <- as.vector(y - X%*%beta)
  LL <- -0.5*n*log(2*pi)
  for (i in uid ) {
    row.i <- which(id==i)
    Ai <- A[row.i,]
    zi <- z[row.i]
    Gi <- Ai%*%Omega%*%t(Ai) + diag(sigma2, nrow=length(row.i))
    LL <- LL -0.5*log(det(Gi)) -0.5*t(zi)%*%solve(Gi)%*%zi
  }
  nb.param <- length(beta) + nb.eta*(nb.eta+1)/2 + 1
  AIC <- -2*LL + 2*nb.param
  BIC <- -2*LL + log(n)*nb.param
  names(beta) <- colnames(X)
  return(list(beta=beta, Omega=Omega, sigma2=sigma2, LL=c(logLik=LL, AIC=AIC, BIC=BIC)))
}
