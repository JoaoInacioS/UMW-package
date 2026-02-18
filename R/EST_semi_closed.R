
# Est UMW ======================================================================
#' Parameter Estimation for the Unit-Modified Weibull(UMW) Distribution
#'
#' Estimates the parameters of the UMW distribution by maximum likelihood.
#'
#' The estimation is performed using a semi-closed form approach, where the
#' parameters \eqn{\gamma} and \eqn{\lambda} are obtained numerically via
#' optimization, while the parameter \eqn{\alpha} is estimated analytically
#' from its likelihood equation. The score vector and the Hessian matrix are
#' computed analytically.
#'
#' @param x Numeric vector with values in the unit interval (0, 1).
#' @param method Optimization method used by \code{\link[stats]{optim}}.
#'   Possible values are \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"} and \code{"SANN"}.
#' @param applic Logical; if TRUE returns the full optimization object.
#' @param start.theta Numeric vector of initial values for
#'   \eqn{(\gamma, \lambda)}.
#'
#' @return
#' If \code{applic = TRUE}, returns a list containing:
#' \itemize{
#'   \item \code{par}: numeric vector of parameter estimates
#'     \eqn{(\hat{\alpha}, \hat{\gamma}, \hat{\lambda})};
#'   \item \code{value}: maximized log-likelihood value;
#'   \item \code{counts}: number of function and gradient evaluations;
#'   \item \code{convergence}: convergence code returned by
#'     \code{\link[stats]{optim}} (0 indicates successful convergence);
#'   \item \code{hessian}: observed information matrix evaluated at the
#'     parameter estimates.
#' }
#'
#' If \code{applic = FALSE}, returns the parameter estimates
#' \eqn{(\hat{\alpha},\hat{\gamma},\hat{\lambda})} if the algorithm converges,
#' or a vector of \code{NA}s otherwise.
#'
#' @examples
#' library(UMW)
#'
#' x <- runif(n = 50,min = 0.001,max = 0.999)
#' Est_UMW(x,method="BFGS")
#'
#' @export
Est_UMW<-function(x,method="BFGS",applic = T,start.theta=c(1,1))
{
  if (!is.numeric(x) || !all(x > 0 & x < 1))
    stop("`x` must be numeric values in the interval (0,1).")
  if (!is.numeric(start.theta) || any(start.theta <= 0))
    stop("`start.theta` must be a numeric vector with positive values.")
  if (!method %in% c("Nelder-Mead", "BFGS", "CG", "SANN"))
    stop("`method` must be one of: 'Nelder-Mead', 'BFGS', 'CG', 'SANN'.")
  mod1<-suppressWarnings(try(optim(par=start.theta,fn=llike_UMW,x=x,
                  method=method,gr=vscore_UMW,hessian=F,m.optim=1.0,
                  control=list(fnscale=-1,reltol=1e-12,
                               maxit = 2000)),T))
  if(class(mod1)=="list"){
    mod1$par[2]<-abs(mod1$par[2])
    alpha<-length(x)/sum((-log(x))^mod1$par[1]/x^mod1$par[2])
    mod1$par<-c(alpha,mod1$par)
    mod1$hessian<-hessian_UMW(theta=mod1$par,x=x)
    tmp2<-test.fun(mod1)
  }
  if(class(tmp2)=="numeric"){
    if(applic==F){tmp2}else{
      mod1$message<-NULL
      return(mod1)}
  }else{if(applic==F){rep(NA,3)}else{return(cat("ALGORITHM DID NOT CONVERGE!"))}}
}

## Log-vero UMW ----------------------------------------------------------------

llike_UMW<-function(theta,x,m.optim = 1)
{
  gamma<-theta[1]
  lambda<-abs(theta[2])
  alpha<-length(x)/sum((-log(x))^gamma/x^lambda)
  log_like<-sum(log(alpha)-(lambda+1)*log(x)+(gamma-1)*log(-log(x))+log(gamma-lambda*log(x))-
                    alpha*((-log(x))^gamma)*x^(-lambda))
  if(m.optim==-1){return(-log_like)}
  if(m.optim==1){return(log_like)}
}

## Score Function MLE ----------------------------------------------------------

vscore_UMW<-function(theta,x,m.optim = 1)
{
  gamma<-theta[1]
  lambda<-abs(theta[2])
  alpha<-length(x)/sum((-log(x))^gamma/x^lambda)
  # alpha
  wt<- (1/alpha)-((-log(x))^gamma/x^lambda)
  # gamma
  rt<- 1/(gamma-log(x)*lambda)-(alpha*(-log(x))^gamma*log(-log(x)))/x^lambda+log(-log(x))
  # lambda
  st <- -log(x)/(gamma-log(x)*lambda)+(alpha*(-log(x))^gamma*log(x))/x^lambda-log(x)
  Ugamma<-sum(rt); Ulambda<-sum(st)
  vetor_score<-c(Ugamma,Ulambda)
  if(m.optim==-1){return(-vetor_score)}
  if(m.optim==1){return(vetor_score)}
}

## Hessian UMW -----------------------------------------------------------------

hessian_UMW<-function(theta,x)
{
  alpha<-theta[1]
  gamma<-theta[2]
  lambda<-theta[3]
  #
  n=length(x)
  I_n<-(rep(1,n))
  d_aa<-(-n/alpha^2)
  d_gg<-(-1/(gamma-log(x)*lambda)^2-(alpha*(-log(x))^gamma*log(-log(x))^2)/x^lambda)
  d_ll<-(-log(x)^2/(gamma-log(x)*lambda)^2-(alpha*(-log(x))^gamma*log(x)^2)/x^lambda)
  d_ag<-(-((-log(x))^gamma*log(-log(x)))/x^lambda)
  d_al<-(((-log(x))^gamma*log(x))/x^lambda)
  d_gl<-(log(x)/(gamma-log(x)*lambda)^2+(alpha*(-log(x))^gamma*log(x)*log(-log(x)))/x^lambda)
  #
  J_aa<-d_aa
  J_gg<-t(d_gg)%*%I_n
  J_ll<-t(d_ll)%*%I_n
  J_ag<-t(d_ag)%*%I_n
  J_al<-t(d_al)%*%I_n
  J_gl<-t(d_gl)%*%I_n
  #
  (matrix(c(J_aa,J_ag,J_al,J_ag,J_gg,J_gl,J_al,J_gl,J_ll),ncol = 3,byrow = T))
}






