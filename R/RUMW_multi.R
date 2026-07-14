
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
#' @param method Optimization method used for parameter estimation.
#'   Possible values are \code{"nlminb"} and \code{"BFGS"}.
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
EST_RQUMW<-function(y, X, Z, method=c("nlminb","BFGS"), g_mu, g_lambda,tau=0.5,
                         applic=T,ginv_mu, ginv_lambda,start.theta=NULL)
{
  method <- match.arg(method)
  if (!is.numeric(y) || !all(y > 0 & y < 1))
    stop("`y` must be numeric values in the interval (0,1).")
  if (!is.numeric(tau) || !all(tau > 0 & tau < 1))
    stop("`tau` must be numeric values in the interval (0,1).")
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
  n_z<-ncol(Z)
  if(is.null(start.theta)){
    df1 <- data.frame(ynew=g_mu(y), X)
    mod.ols1 <- try(lm(ynew~.-1,data=df1), TRUE)
    startbeta1 <- mod.ols1$coefficients
    if(ncol(Z)== 1 & var(Z[,1]) == 0){repZ<-1}else{repZ<-rep(0,n_z)}
    start.theta <- c(startbeta1, 1, repZ)
  }
  if(method=="BFGS"){
    mod<-suppressWarnings(try(optim(par=start.theta,fn=llike_RQUMW,y=y,X=X,Z=Z,
                gr = vscore_RQUMW,ginv_mu = ginv_mu, ginv_lambda = ginv_lambda,
                g_mu=g_mu, g_lambda=g_lambda,method="BFGS",tau=tau,hessian=F,
                m.optim=1L,control=list(fnscale=-1,reltol=1e-8)),T))
  }
  if(method=="nlminb"){
    mod<-suppressWarnings(try(nlminb(start=start.theta, objective=llike_RQUMW,
                   gradient=vscore_RQUMW,y=y, X=X, Z=Z,
                   ginv_mu=ginv_mu, ginv_lambda=ginv_lambda, g_mu=g_mu,
                   g_lambda=g_lambda, tau=tau, m.optim=-1L,
                   control=list(rel.tol = 1e-8,trace=FALSE)),T))
  }
  if (!inherits(mod, "try-error")){
    if(method=="nlminb"){mod$value<--mod$objective}
    mod$hessian<-hessian_RQUMW(theta = mod$par,y = y,X = X,Z = Z,tau = tau,g_mu=g_mu,
                       g_lambda=g_lambda,ginv_mu = ginv_mu,ginv_lambda = ginv_lambda)
    tmp2<-test.fun(mod)
    if(is.numeric(tmp2)){
      if(applic==F){
        return(tmp2)
        }else{
        mod$message<-NULL
        return(mod)}
    }else{if(applic==F){return(rep(NA,length(start.theta)*2))
      }else{
      cat(tmp2)
      return(tmp2)}}
  }else{if(applic==F){return(rep(NA,length(start.theta)*2))
    }else{return(mod)}
  }
}

## Log-vero RQUMW --------------------------------------------------------------

llike_RQUMW <- function(theta, y, X, Z, tau, ginv_mu, ginv_lambda,
                             g_mu, g_lambda, m.optim=1) {
  n_beta_mu <- ncol(X)
  n_beta_lambda <- ncol(Z)
  beta_mu <- theta[1:n_beta_mu]
  gamma_i <- theta[(n_beta_mu + 1)]
  beta_lambda <- theta[(n_beta_mu + 2):(n_beta_mu + 1 + n_beta_lambda)]
  if (gamma_i <= 0 || (ncol(Z) == 1 && var(Z[, 1]) == 0 && beta_lambda <= 0)) {
    return(m.optim * (-1e+10))
  }
  mu_i <- ginv_mu(as.vector(X %*% beta_mu))
  log_y <- log(y)
  lambda_i <- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i) * log(tau))/((-log(mu_i))^(gamma_i))
  log_like <- sum(log(alpha_i) - (lambda_i + 1) * log_y +
                    (gamma_i - 1) * log(-log_y) + log(gamma_i - lambda_i *
                    log_y) - alpha_i * ((-log_y)^gamma_i) * y^(-lambda_i))
  if (is.na(log_like) || is.infinite(log_like)) {return(m.optim * (-1e+10))}
  return(m.optim * log_like)
}

## Score Function MLE ----------------------------------------------------------

vscore_RQUMW <- function(theta, y, X, Z, tau, g_mu, g_lambda,
                              ginv_mu, ginv_lambda, m.optim=1,vsmatrix=F) {
  n <- length(y)
  D1_mu <- Deriv::Deriv(g_mu)
  if (isTRUE(all.equal(g_lambda(c(1,5,10)), c(1,5,10)))) {
    D1_lambda <- function(x) rep(1, length(x))
  }else{D1_lambda <- Deriv::Deriv(g_lambda)}
  n_beta_mu <- ncol(X)
  n_beta_lambda <- ncol(Z)
  n_t <- n_beta_mu + 1 + n_beta_lambda
  beta_mu <- theta[1:n_beta_mu]
  gamma_i <- theta[(n_beta_mu + 1)]
  beta_lambda <- theta[(n_beta_mu + 2):(n_t)]
  if (gamma_i <= 0 || (ncol(Z) == 1 && var(Z[, 1]) == 0 && beta_lambda <= 0)) {
    if (vsmatrix == FALSE) {
      return(m.optim * rep(-1e+05, n_t))
    }
    else {
      m_score_penalizado <- matrix(-1e+05, nrow = n, ncol = n_t)
      return(m.optim * m_score_penalizado)
    }
  }
  mu_i <- ginv_mu(as.vector(X %*% beta_mu))
  log_y <- log(y)
  log_mu <- log(mu_i)
  lambda_i <- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i) * log(tau))/((-log_mu)^(gamma_i))
  wt <- ((y^lambda_i * (-log_mu)^gamma_i + log(tau) * (-log_y)^gamma_i *
            mu_i^lambda_i) * (lambda_i * log_mu - gamma_i))/(y^lambda_i *
                                                                  mu_i * (-log_mu)^gamma_i * log_mu)
  rt <- 1/(gamma_i - log_y * lambda_i) + log(-log_y) - log(-log_mu) +
    (mu_i^lambda_i * log(tau) * (-log_y)^gamma_i * (log(-log_y) -
                                                       log(-log_mu)))/((-log_mu)^gamma_i * y^lambda_i)
  st <- -log_y/(gamma_i - log_y * lambda_i) - log_y + log_mu +
    (mu_i^lambda_i * log(tau) * (log_mu - log_y) * (-log_y)^gamma_i)/((-log_mu)^gamma_i *
                                                                             y^lambda_i)
  if (vsmatrix == F) {
    Ubeta <- as.vector(t(X) %*% (wt / D1_mu(mu_i)))
    Ugamma <-  rep(1,n) %*% rt
    Ulambda <- as.vector(t(Z) %*% (st / D1_lambda(lambda_i)))
    vetor_score <- c(Ubeta, Ugamma, Ulambda)
    vetor_score[is.na(vetor_score) | is.infinite(vetor_score)] <- -1e+05
    return(m.optim * vetor_score)
  }
  else {
    Ubeta_m <- X * ((1/D1_mu(mu_i)) * wt)
    Ugamma_m <- rt
    Ulambda_m <- Z * ((1/D1_lambda(lambda_i)) * st)
    m_score <- as.matrix(cbind(Ubeta_m, Ugamma_m, Ulambda_m))
    m_score[is.na(m_score) | is.infinite(m_score)] <- -1e+05
    return(m.optim * m_score)
  }
}

## Hessian RQUMW ---------------------------------------------------------------

hessian_RQUMW <- function(theta, y, X, Z, tau, g_mu, g_lambda,
         ginv_mu, ginv_lambda, m.optim=1)
{
  n <- length(y)
  n_beta_mu <- ncol(X)
  n_beta_lambda <- ncol(Z)
  #
  beta_mu <- theta[(1:n_beta_mu)]
  gamma_i <- theta[(n_beta_mu+1)]
  beta_lambda <- theta[(n_beta_mu+2):(n_beta_mu+1+n_beta_lambda)]
  #
  mu_i <- ginv_mu(as.vector(X %*% beta_mu))
  log_mu <- log(mu_i)
  log_y <- log(y)
  lambda_i<- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i)*log(tau))/((-log_mu)^(gamma_i))
  #
  D1_mu <- Deriv::Deriv(g_mu)
  D2_mu <- Deriv::Deriv(D1_mu)
  if(isTRUE(all.equal(g_lambda(c(1,5,10)), c(1,5,10)))){
    D1_lambda <-function(x) rep(1, length(x))
    D2_lambda <-function(x) rep(0, length(x))
  }else{
    D1_lambda <- Deriv::Deriv(g_lambda)
    D2_lambda <- Deriv::Deriv(D1_lambda)
  }
  #
  t1_mu <- D1_mu(mu_i)
  t2_mu <- -(D2_mu(mu_i) / D1_mu(mu_i)^2)
  t1_lambda <- D1_lambda(lambda_i)
  t2_lambda <- -(D2_lambda(lambda_i) / D1_lambda(lambda_i)^2)
  A_1 <- mu_i^lambda_i*log(tau)*(-log_y)^gamma_i
  A_2 <- log(-log_y)-log(-log_mu)
  #
  d_GG<-((A_1*A_2^2)/((-log_mu)^gamma_i*y^lambda_i)-1/(gamma_i-log_y*lambda_i)^2)
  d_LL<-((A_1*(log_mu-log_y)^2/((-log_mu)^gamma_i*y^lambda_i))-log_y^2/(gamma_i-log_y*lambda_i)^2)
  d_GL<-(log_y/(gamma_i-log_y*lambda_i)^2+(A_1*A_2*(log_mu-log_y))/((-log_mu)^gamma_i*y^lambda_i))
  d_Gb<-((-((A_1*(A_2*(gamma_i-log_mu*lambda_i)+1))/((-log_mu)^gamma_i*y^lambda_i))-1)/(mu_i*log_mu))
  d_Lb<-(((-log_mu)^gamma_i*y^lambda_i+A_1)/(mu_i*(-log_mu)^gamma_i*y^lambda_i)-(A_1*(log_y-log_mu)*(log_mu*lambda_i-gamma_i))/(mu_i*(-log_mu)^gamma_i*log_mu*y^lambda_i))
  d_bb<-(gamma_i/(mu_i^2*log_mu^2)-(log_mu*gamma_i*(A_1*(2*lambda_i-1)-(-log_mu)^gamma_i*y^lambda_i)-A_1*gamma_i*(gamma_i+1))/(mu_i^2*(-log_mu)^gamma_i*log_mu^2*y^lambda_i)) -
    ((A_1*(1-lambda_i)+(-log_mu)^gamma_i*y^lambda_i)*lambda_i)/(mu_i^2*(-log_mu)^gamma_i*y^lambda_i)
  # mu
  wt <- ((y^lambda_i * (-log_mu)^gamma_i + log(tau)*(-log_y)^gamma_i * mu_i^lambda_i)*(lambda_i*log_mu-gamma_i)) /
    (y^lambda_i*mu_i*(-log_mu)^gamma_i*log_mu)
  # gamma
  rt <- 1/(gamma_i - log_y*lambda_i) + log(-log_y) - log(-log_mu) +
    (mu_i^lambda_i * log(tau) * (-log_y)^gamma_i * (log(-log_y) - log(-log_mu))) /
    ((-log_mu)^gamma_i * y^lambda_i)
  # lambda
  st <- -log_y/(gamma_i - log_y*lambda_i) - log_y + log_mu +
    (mu_i^lambda_i * log(tau) * (log_mu-log_y) * (-log_y)^gamma_i) /
    ((-log_mu)^gamma_i * y^lambda_i)
  #
  J_bb <- t(X) %*% (((d_bb / t1_mu + wt * t2_mu) / t1_mu) * X)
  J_GG <- rep(1, n) %*%  d_GG
  J_LL <- t(Z) %*% (((d_LL / t1_lambda + st * t2_lambda) / t1_lambda) * Z)
  J_Gb <- t(X) %*% (d_Gb / t1_mu)
  J_Lb <- t(X) %*% ((d_Lb / (t1_mu * t1_lambda)) * Z)
  J_GL <- t(Z) %*% (d_GL / t1_lambda)
  #
  hessian <- rbind(
    cbind(J_bb, J_Gb,J_Lb),
    cbind(t(J_Gb), J_GG, t(J_GL)),
    cbind(t(J_Lb), J_GL, J_LL)
  )
  return(m.optim * hessian)
}

