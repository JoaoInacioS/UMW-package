
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
#' @param method Optimization method used for parameter estimation.
#'   Possible values are \code{"nlminb"} and \code{"BFGS"}.
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
Est_UMW<-function(x,method=c("nlminb","BFGS"),applic = T,start.theta=c(1,1))
{
  method <- match.arg(method)
  if (!is.numeric(x) || !all(x > 0 & x < 1))
    stop("`x` must be numeric values in the interval (0,1).")
  if (!is.numeric(start.theta) || any(start.theta <= 0))
    stop("`start.theta` must be a numeric vector with positive values.")
  if(method=="BFGS"){
    mod1<-suppressWarnings(try(optim(par=start.theta,fn=llike_UMW,x=x,
                  method="BFGS",gr=vscore_UMW,hessian=F,m.optim=1L,
                  control=list(fnscale=-1,reltol=1e-8)),T))
  }
  if(method=="nlminb"){
    mod1<-suppressWarnings(try(nlminb(start=start.theta, objective=llike_UMW,
                       gradient=vscore_UMW,x=x,m.optim=-1L,
                       control=list(rel.tol = 1e-8,trace=FALSE)),T))
  }
  if(!inherits(mod1, "try-error")){
    if(method=="nlminb"){mod1$value<--mod1$objective}
    alpha<-length(x)/sum((-log(x))^mod1$par[1]/x^mod1$par[2])
    mod1$par<-c(alpha,mod1$par)
    mod1$hessian<-hessian_UMW(theta=mod1$par,x=x)
    tmp2<-test.fun(mod1)
    if(is.numeric(tmp2)){
      if(applic==F){tmp2}else{
        mod1$message<-NULL
        return(mod1)}
    }else{if(applic==F){rep(NA,6)}else{
      cat(tmp2)
      return(tmp2)}}
  }else{if(applic==F){rep(NA,6)}else{return(mod1)}}
}

## Log-vero UMW ----------------------------------------------------------------

llike_UMW<-function(theta,x,m.optim = 1)
{
  gamma<-theta[1]
  lambda<-theta[2]
  if (gamma <= 0 || lambda <= 0) {return(m.optim * (-1e+10))}
  log_x<-log(x)
  alpha<-length(x)/sum((-log_x)^gamma/x^lambda)
  log_like<-sum(log(alpha)-(lambda+1)*log_x+(gamma-1)*log(-log_x)+log(gamma-lambda*log_x)-
                    alpha*((-log_x)^gamma)*x^(-lambda))
  if (is.na(log_like) || is.infinite(log_like)) {return(m.optim * (-1e+10))}
  return(m.optim * log_like)
}

## Score Function MLE ----------------------------------------------------------

vscore_UMW<-function(theta,x,m.optim = 1)
{
  gamma<-theta[1]
  lambda<-theta[2]
  if (gamma <= 0 || lambda <= 0) {return(m.optim * rep(-1e+10,2))}
  log_x<-log(x)
  alpha<-length(x)/sum((-log_x)^gamma/x^lambda)
  # alpha
  wt<- (1/alpha)-((-log_x)^gamma/x^lambda)
  # gamma
  rt<- 1/(gamma-log_x*lambda)-(alpha*(-log_x)^gamma*log(-log_x))/x^lambda+log(-log_x)
  # lambda
  st <- -log_x/(gamma-log_x*lambda)+(alpha*(-log_x)^gamma*log_x)/x^lambda-log_x
  Ugamma<-sum(rt); Ulambda<-sum(st)
  vetor_score<-c(Ugamma,Ulambda)
  vetor_score[is.na(vetor_score) | is.infinite(vetor_score)] <- -1e+05
  return(m.optim * vetor_score)
}

## Hessian UMW -----------------------------------------------------------------

hessian_UMW<-function(theta,x,m.optim = 1)
{
  alpha<-theta[1]
  gamma<-theta[2]
  lambda<-theta[3]
  #
  n=length(x)
  I_n<-(rep(1,n))
  d_aa<-(-n/alpha^2)
  log_x<-log(x)
  d_gg<-(-1/(gamma-log_x*lambda)^2-(alpha*(-log_x)^gamma*log(-log_x)^2)/x^lambda)
  d_ll<-(-log_x^2/(gamma-log_x*lambda)^2-(alpha*(-log_x)^gamma*log_x^2)/x^lambda)
  d_ag<-(-((-log_x)^gamma*log(-log_x))/x^lambda)
  d_al<-(((-log_x)^gamma*log_x)/x^lambda)
  d_gl<-(log_x/(gamma-log_x*lambda)^2+(alpha*(-log_x)^gamma*log_x*log(-log_x))/x^lambda)
  #
  J_aa<-d_aa
  J_gg<-t(d_gg)%*%I_n
  J_ll<-t(d_ll)%*%I_n
  J_ag<-t(d_ag)%*%I_n
  J_al<-t(d_al)%*%I_n
  J_gl<-t(d_gl)%*%I_n
  #
  hessian<-matrix(c(J_aa,J_ag,J_al,J_ag,J_gg,J_gl,J_al,J_gl,J_ll),ncol = 3,byrow = T)
  return(m.optim * hessian)
}







