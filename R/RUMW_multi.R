
#' Estimation of the Regression Quantile Unit-Modified Weibull (RQ-UMW) Model
#'
#' Estimates the parameters of the Regression Quantile UMW model via
#' numerical maximization of the log-likelihood function. The score vector
#' and the Hessian matrix are computed analytically.
#'
#' @param y Numeric vector of observed responses. Values must lie in
#'   the open interval \eqn{(0,1)}.
#' @param X Numeric matrix of covariates associated with the
#'   \eqn{\mu} parameter. If an intercept is required, it must be explicitly
#'   included as a column in \code{X}.
#' @param Z Numeric matrix of covariates associated with the \eqn{\lambda}
#'   parameter. If an intercept is required, it must be explicitly
#'   included as a column in \code{Z}.
#' @param method Optimization method used by \code{\link[stats]{optim}}.
#'   Possible values are \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"} and \code{"SANN"}.
#' @param g_mu Link function for the \eqn{\mu} parameter.
#' @param g_lambda Link function for the \eqn{\lambda} parameter.
#' @param tau Quantile level, must be in \eqn{(0,1)} (default is 0.5).
#' @param applic Logical. If \code{TRUE}, the full output of the estimation
#'   procedure is returned. If \code{FALSE}, only the vector of parameter
#'   estimates is returned.
#' @param ginv_mu Inverse link function for the \eqn{\mu} parameter.
#' @param ginv_lambda Inverse link function for the \eqn{\lambda} parameter.
#' @param start.theta Optional numeric vector of initial values for the
#'   regression parameters. If \code{NULL} (default), suitable starting
#'   values are internally determined.
#'
#' @return
#' If \code{applic = TRUE}, returns a list containing:
#' \itemize{
#'   \item \code{par}: numeric vector of parameter estimates;
#'   \item \code{value}: maximized log-likelihood value;
#'   \item \code{counts}: number of function and gradient evaluations;
#'   \item \code{convergence}: convergence code returned by
#'     \code{\link[stats]{optim}} (0 indicates successful convergence);
#'   \item \code{hessian}: observed information matrix evaluated at the
#'     parameter estimates.
#' }
#'
#' If \code{applic = FALSE}, returns a numeric vector with the estimated
#' regression parameters.
#'
#' @examples
#' library(UMW)
#'
#' n <- 100
#' y<-runif(n = n,min = 0.001,max = 0.999)
#' X<-cbind(intercept = 1,X1 = runif(n),X2 = runif(n))
#' Z<-matrix(1,n,1)
#' # Logit:
#' g_mu <- function(mu) {log(mu / (1 - mu))}
#' ginv_mu <- function(eta) {1 / (1 + exp(-eta))}
#' # Identity:
#' g_lambda <- ginv_lambda <- function(x) x
#'
#' EST_RQUMW(y=y,X=X,Z=Z,g_mu = g_mu,g_lambda = g_lambda,
#' ginv_mu = ginv_mu,ginv_lambda = ginv_lambda,
#' method = "BFGS",tau = 0.5)
#'
#' @export
EST_RQUMW<-function(y, X, Z, method="BFGS", g_mu, g_lambda,tau=0.5,
                         applic=T,ginv_mu, ginv_lambda,start.theta=NULL)
{
  if (!is.numeric(y) || !all(y > 0 & y < 1))
    stop("`y` must be numeric values in the interval (0,1).")
  if (!is.numeric(tau) || !all(tau > 0 & tau < 1))
    stop("`tau` must be numeric values in the interval (0,1).")
  if (!method %in% c("Nelder-Mead", "BFGS", "CG", "SANN"))
    stop("`method` must be one of: 'Nelder-Mead', 'BFGS', 'CG', 'SANN'.")
  if (!is.matrix(X) || !is.matrix(Z)) stop("X and Z must be matrices.")
  if (any(c(ncol(X), ncol(Z)) < 1))
    stop("X and Z must each have at least one column.")
  if (!is.null(start.theta) && (length(start.theta) != ncol(X) + 1 + ncol(Z))) {
    stop("start.theta must be NULL or a numeric vector values with length equal to the number of parameters.")
  }
  if (!inherits(g_mu, "function") || !inherits(g_lambda, "function")) {
    stop("g_mu and g_lambda must be functions.")
  }
  if (!inherits(ginv_mu, "function") || !inherits(ginv_lambda, "function")) {
    stop("ginv_mu and ginv_lambda must be functions.")
  }
  if(is.null(start.theta)){
    # lm: Ordinary Least Squares
    df1 <- data.frame(ynew=g_mu(y), X)
    if(any((X[,1])== 1) & var(X[,1]) == 0){
      mod.ols1<-try(quantreg::rq(ynew~.,data=df1[,-2],tau = tau),T)
    }else{mod.ols1<-try(quantreg::rq(ynew~.-1,data=df1,tau = tau),T)}
    startbeta1<-mod.ols1$coefficients
    start.theta<-c(startbeta1,rep(1,1+ncol(Z)))
  }
  mod<-suppressWarnings(try(optim(par=start.theta,fn=llike_RQUMW,y=y,X=X,Z=Z,
                 gr = vscore_RQUMW,ginv_mu = ginv_mu, ginv_lambda = ginv_lambda,
                 g_mu=g_mu, g_lambda=g_lambda,method=method,tau=tau,hessian=F,
                 m.optim=1.0,control=list(fnscale=-1,reltol=1e-12,
                                          maxit = 2000)),T))
  if(class(mod)=="list"){
    if(ncol(Z)== 1 & var(Z[,1]) == 0){mod$par[ncol(X)+2] <- abs(mod$par[ncol(X)+2])}
    mod$par[ncol(X)+1] <- abs(mod$par[ncol(X)+1])
    mod$hessian<-hessian_RQUMW(theta = mod$par,y = y,X = X,Z = Z,tau = tau,g_mu=g_mu,
                               g_lambda=g_lambda,ginv_mu = ginv_mu,ginv_lambda = ginv_lambda)
  }
  tmp2<-test.fun(mod)
  if(class(tmp2)=="numeric"){
    if(applic==F){tmp2}else{
      mod$message<-NULL
      return(mod)}
  }else{if(applic==F){rep(NA,length(start.theta))}else{return(cat("ALGORITHM DID NOT CONVERGE!"))}}
}

## Log-vero RQUMW --------------------------------------------------------------

llike_RQUMW <- function(theta, y, X, Z, tau, ginv_mu, ginv_lambda,
                             g_mu, g_lambda, m.optim=1) {
  n_beta_mu <- ncol(X)
  n_beta_lambda <- ncol(Z)
  #
  beta_mu <- theta[1:n_beta_mu]
  beta_gamma <- theta[(n_beta_mu+1)]
  beta_lambda <- theta[(n_beta_mu+2):(n_beta_mu+1+n_beta_lambda)]
  if(ncol(Z)== 1 & var(Z[,1]) == 0){beta_lambda <- abs(beta_lambda)}
  gamma_i <- abs(beta_gamma)
  #
  mu_i    <- ginv_mu(as.vector(X %*% beta_mu))
  lambda_i<- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i)*log(tau))/((-log(mu_i))^(gamma_i))
  #
  log_like <- sum(log(alpha_i) - (lambda_i+1)*log(y) + (gamma_i-1)*log(-log(y)) +
                    log(gamma_i - lambda_i*log(y)) - alpha_i*((-log(y))^gamma_i)*y^(-lambda_i))
  if(m.optim==-1){return(-log_like)}
  else{return(log_like)}
}

## Score Function MLE ----------------------------------------------------------

vscore_RQUMW <- function(theta, y, X, Z, tau, g_mu, g_lambda,
                              ginv_mu, ginv_lambda, m.optim=1,vsmatrix=F) {
  D1_mu <- Deriv::Deriv(g_mu)
  D1_lambda <- Deriv::Deriv(g_lambda)
  if(D1_lambda(10)==1){D1_lambda <-function(x) rep(1, length(x))}
  #
  n_beta_mu    <- ncol(X)
  n_beta_lambda<- ncol(Z)
  #
  beta_mu    <- theta[1:n_beta_mu]
  beta_gamma <- theta[(n_beta_mu+1)]
  beta_lambda<- theta[(n_beta_mu+2):(n_beta_mu+1+n_beta_lambda)]
  if(ncol(Z)== 1 & var(Z[,1]) == 0){beta_lambda <- abs(beta_lambda)}
  gamma_i <- abs(beta_gamma)
  #
  mu_i    <- ginv_mu(as.vector(X %*% beta_mu))
  lambda_i<- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i)*log(tau))/((-log(mu_i))^(gamma_i))
  # mu
  wt <- ((y^lambda_i * (-log(mu_i))^gamma_i + log(tau)*(-log(y))^gamma_i * mu_i^lambda_i)*(lambda_i*log(mu_i)-gamma_i)) /
    (y^lambda_i*mu_i*(-log(mu_i))^gamma_i*log(mu_i))
  # gamma
  rt <- 1/(gamma_i - log(y)*lambda_i) + log(-log(y)) - log(-log(mu_i)) +
    (mu_i^lambda_i * log(tau) * (-log(y))^gamma_i * (log(-log(y)) - log(-log(mu_i)))) /
    ((-log(mu_i))^gamma_i * y^lambda_i)
  # lambda
  st <- -log(y)/(gamma_i - log(y)*lambda_i) - log(y) + log(mu_i) +
    (mu_i^lambda_i * log(tau) * (log(mu_i)-log(y)) * (-log(y))^gamma_i) /
    ((-log(mu_i))^gamma_i * y^lambda_i)
  if(vsmatrix == F){
    Ubeta <- as.vector(t(X) %*% (diag(1/D1_mu(mu_i))) %*% wt)
    Ugamma<-sum(rt)
    Ulambda<-as.vector(t(Z) %*% (diag(1/D1_lambda(lambda_i))) %*% st)
    vetor_score<-c(Ubeta,Ugamma,Ulambda)
    if(m.optim==-1){return(-vetor_score)}
    if(m.optim==1){return(vetor_score)}
  }else{
    Ubeta_m <- X * ((1/D1_mu(mu_i)) * wt)
    Ugamma_m <- rt
    Ulambda_m <- Z * ((1/D1_lambda(lambda_i)) * st)
    m_score<-as.matrix(cbind(Ubeta_m,Ugamma_m,Ulambda_m))
    if(m.optim==-1){return(-m_score)}
    if(m.optim==1){return(m_score)}
  }
  return(vetor_score)
}

## Hessian RQUMW ---------------------------------------------------------------

hessian_RQUMW <- function(theta, y, X, Z, tau, g_mu, g_lambda,
         ginv_mu, ginv_lambda)
{
  n <- length(y)
  n_beta_mu    <- ncol(X)
  n_beta_lambda<- ncol(Z)
  #
  beta_mu    <- theta[(1:n_beta_mu)]
  beta_gamma <- theta[(n_beta_mu+1)]
  beta_lambda<- theta[(n_beta_mu+2):(n_beta_mu+1+n_beta_lambda)]
  if(ncol(Z)== 1 & var(Z[,1]) == 0){beta_lambda <- abs(beta_lambda)}
  gamma_i <- abs(beta_gamma)
  #
  mu_i    <- ginv_mu(as.vector(X %*% beta_mu))
  lambda_i<- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i)*log(tau))/((-log(mu_i))^(gamma_i))
  #
  D1_mu <- Deriv::Deriv(g_mu)
  D2_mu <- Deriv::Deriv(D1_mu)
  D1_lambda <- Deriv::Deriv(g_lambda)
  if(D1_lambda(10)==1){
    D1_lambda <-function(x) rep(1, length(x))
    D2_lambda <-function(x) rep(0, length(x))
  }else{D2_lambda <- Deriv::Deriv(D1_lambda)}
  #
  T1_mu <- diag(1/D1_mu(mu_i))
  T2_mu <- diag(-(D2_mu(mu_i)/(D1_mu(mu_i)^2)))
  T1_lambda <- diag(1/D1_lambda(lambda_i))
  T2_lambda <- diag(-(D2_lambda(lambda_i)/(D1_lambda(lambda_i)^2)))
  A_1 <- mu_i^lambda_i*log(tau)*(-log(y))^gamma_i
  A_2 <- log(-log(y))-log(-log(mu_i))
  #
  d_GG<-((A_1*A_2^2)/((-log(mu_i))^gamma_i*y^lambda_i)-1/(gamma_i-log(y)*lambda_i)^2)
  d_LL<-((A_1*(log(mu_i)-log(y))^2/((-log(mu_i))^gamma_i*y^lambda_i))-log(y)^2/(gamma_i-log(y)*lambda_i)^2)
  d_GL<-(log(y)/(gamma_i-log(y)*lambda_i)^2+(A_1*A_2*(log(mu_i)-log(y)))/((-log(mu_i))^gamma_i*y^lambda_i))
  d_Gb<-((-((A_1*(A_2*(gamma_i-log(mu_i)*lambda_i)+1))/((-log(mu_i))^gamma_i*y^lambda_i))-1)/(mu_i*log(mu_i)))
  d_Lb<-(((-log(mu_i))^gamma_i*y^lambda_i+A_1)/(mu_i*(-log(mu_i))^gamma_i*y^lambda_i)-(A_1*(log(y)-log(mu_i))*(log(mu_i)*lambda_i-gamma_i))/(mu_i*(-log(mu_i))^gamma_i*log(mu_i)*y^lambda_i))
  d_bb<-(gamma_i/(mu_i^2*log(mu_i)^2)-(log(mu_i)*gamma_i*(A_1*(2*lambda_i-1)-(-log(mu_i))^gamma_i*y^lambda_i)-A_1*gamma_i*(gamma_i+1))/(mu_i^2*(-log(mu_i))^gamma_i*log(mu_i)^2*y^lambda_i)) -
    ((A_1*(1-lambda_i)+(-log(mu_i))^gamma_i*y^lambda_i)*lambda_i)/(mu_i^2*(-log(mu_i))^gamma_i*y^lambda_i)
  # mu
  wt <- ((y^lambda_i * (-log(mu_i))^gamma_i + log(tau)*(-log(y))^gamma_i * mu_i^lambda_i)*(lambda_i*log(mu_i)-gamma_i)) /
    (y^lambda_i*mu_i*(-log(mu_i))^gamma_i*log(mu_i))
  # gamma
  rt <- 1/(gamma_i - log(y)*lambda_i) + log(-log(y)) - log(-log(mu_i)) +
    (mu_i^lambda_i * log(tau) * (-log(y))^gamma_i * (log(-log(y)) - log(-log(mu_i)))) /
    ((-log(mu_i))^gamma_i * y^lambda_i)
  # lambda
  st <- -log(y)/(gamma_i - log(y)*lambda_i) - log(y) + log(mu_i) +
    (mu_i^lambda_i * log(tau) * (log(mu_i)-log(y)) * (-log(y))^gamma_i) /
    ((-log(mu_i))^gamma_i * y^lambda_i)
  #
  J_bb<-t(X)%*%diag(c(d_bb%*%T1_mu+(wt)%*%T2_mu))%*%T1_mu%*%X
  J_GG<-(d_GG)%*%rep(1,n)
  J_LL<-t(Z)%*%diag(c(d_LL%*%T1_lambda+(st)%*%T2_lambda))%*%T1_lambda%*%Z
  J_Gb<-t(X)%*%T1_mu%*%d_Gb
  J_Lb<-t(X)%*%diag(c(d_Lb%*%T1_mu%*%T1_lambda))%*%Z
  J_GL<-t(Z)%*%T1_lambda%*%d_GL
  #
  hessian <- rbind(
    cbind(J_bb, J_Gb,J_Lb),
    cbind(t(J_Gb), J_GG, t(J_GL)),
    cbind(t(J_Lb), J_GL, J_LL)
  )
  return(hessian)
}

