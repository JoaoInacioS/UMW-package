# Density and Cumulative =======================================================

#' Probability Density Function of the Unit-Modified Weibull(UMW) Distribution
#'
#' Computes the probability density function (PDF) of the UMW distribution
#' for values in the unit interval.
#'
#' The PDF is given by:
#'
#' \deqn{
#' f(x;\alpha,\gamma,\lambda) =
#' \frac{\alpha}{\log(x)\,x^{\lambda+1}}
#' (-\log x)^{\gamma}
#' (\lambda \log x - \gamma)
#' \exp\left\{-\alpha (-\log x)^{\gamma} x^{-\lambda}\right\},
#' \quad 0 < x < 1,
#' }
#' where \eqn{\log(\cdot)} denotes the natural logarithm and
#' \eqn{\alpha}, \eqn{\gamma}, \eqn{\lambda > 0} are real-valued parameters.
#'
#' @param x Numeric vector with values in the unit interval (0, 1).
#' @param theta Numeric vector of length 3 containing the parameter values
#'   \eqn{(\alpha,\gamma,\lambda)}, all assumed to be positive.
#'
#' @return A numeric vector of the same length as `x` containing the values
#'   of the probability density function evaluated at `x`.
#'
#' @examples
#' library(UMW)
#'
#' x <- runif(n = 50,min = 0.001,max = 0.999)
#' theta <- c(alpha = 0.7, gamma = 1.3, lambda = 0.5)
#' dUMW(x, theta)
#'
#' @export
dUMW<-function(x,theta)
{
  if (!is.numeric(x) || !all(x > 0 & x < 1))
    stop("`x` must be numeric values in the interval (0,1).")
  if (!is.numeric(theta) || length(theta) != 3 || !all(theta > 0))
    stop("`theta` must be a numeric vector of length 3 with positive values.")
  alpha<-theta[1]
  gamma<-theta[2]
  lambda<-theta[3]
  f<- (alpha/(log(x)*x^(lambda+1)))*((-log(x))^gamma)*(lambda*log(x)-gamma)*
    exp(-alpha*((-log(x))^gamma)*x^(-lambda))
  return(f)
}

#' Cumulative Distribution Function of the Unit-Modified Weibull(UMW) Distribution
#'
#' Computes the cumulative distribution function (CDF) of the UMW distribution
#' for values in the unit interval.
#'
#' The CDF is given by:
#'
#' \deqn{
#' F(x;\alpha,\gamma,\lambda) =
#' \exp\left\{-\alpha (-\log x)^{\gamma} x^{-\lambda}\right\},
#' \quad 0 < x < 1,
#' }
#'
#' where \eqn{\log(\cdot)} denotes the natural logarithm and
#' \eqn{\alpha}, \eqn{\gamma}, \eqn{\lambda > 0} are real-valued parameters.
#'
#' @param x Numeric vector with values in the unit interval \eqn{(0,1)}.
#' @param theta Numeric vector of length 3 containing the parameters
#'   \eqn{(\alpha,\gamma,\lambda)}.
#'
#' @return A numeric vector containing the value of the cumulative distribution
#'   function evaluated at each element of \code{x}.
#'
#' @examples
#' library(UMW)
#'
#' x <- runif(n = 50,min = 0.001,max = 0.999)
#' theta <- c(alpha = 0.7, gamma = 1.3, lambda = 0.5)
#' pUMW(x, theta)
#'
#' @export
pUMW<-function(x,theta)
{
  if (!is.numeric(x) || !all(x > 0 & x < 1))
    stop("`x` must be numeric values in the interval (0,1).")
  if (!is.numeric(theta) || length(theta) != 3 || !all(theta > 0))
    stop("`theta` must be a numeric vector of length 3 with positive values.")
  alpha<-theta[1]
  gamma<-theta[2]
  lambda<-theta[3]
  f<- exp(-alpha*((-log(x))^gamma)*(x^(-lambda)))
  return(f)
}

# generate UMW observations ====================================================

#' Random Generation for the Unit-Modified Weibull(UMW) Distribution
#'
#' Generates random observations from the UMW distribution using
#' numerical inversion based on the cumulative distribution function.
#'
#' @param n Integer indicating the number of observations to be generated.
#' @param theta Numeric vector of length 3 containing the parameter values
#'   \eqn{(\alpha,\gamma,\lambda)}, all assumed to be positive.
#' @param set_seed Integer. Optional random seed for reproducibility.
#'
#' @return A numeric vector of length \code{n} containing random observations
#'   from the UMW distribution.
#'
#' @details
#' The random variates are generated via the inversion method, where
#' each observation is obtained by numerically solving the defining
#' equation of the UMW cumulative distribution function using
#' \code{\link[stats]{uniroot}}.
#'
#' @examples
#' theta <- c(alpha = 0.7, gamma = 1.3, lambda = 0.5)
#' rUMW(n = 50,theta = theta,set_seed = 123)
#'
#' @export
rUMW<-function(n,theta,set_seed=NULL)
{
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n))
    stop("`n` must be a single positive integer.")
  if (!is.numeric(theta) || length(theta) != 3 || !all(theta > 0))
    stop("`theta` must be a numeric vector of length 3 with positive values.")
  if (!is.null(set_seed)) {set.seed(set_seed)}
  alpha<-theta[1]
  gamma<-theta[2]
  lambda<-theta[3]
  eq<-function(x,parms)
    {((-log(x))^parms[3])*(x^(-parms[4]))+(1/parms[2])*log(parms[1])}
  y<-c()
  j<-1
  while(j<=n){
    u<-runif(1)
    parms<-c(u,alpha,gamma,lambda)
    tmp<-try(suppressWarnings(stats::uniroot(f=eq,lower=0.0001,upper=0.9999,
                                             parms=parms)),T)
    if(class(tmp)=="list"){
      y[j]<-tmp$root
      j<-j+1
    }
  }
  return(y)
}

# quantile function ============================================================

#' Quantile Function of the Unit-Modified Weibull(UMW) Distribution
#'
#' Computes the quantile function (inverse CDF) of the UMW distribution
#' for given probability levels.
#'
#' The quantile function is defined implicitly as the solution
#' \eqn{x \in (0,1)} of the equation:
#'
#' \deqn{
#' \exp\left\{-\alpha (-\log x)^{\gamma} x^{-\lambda}\right\} = p,
#' \quad 0 < p < 1.
#' }
#'
#' Equivalently, the quantile \eqn{Q(p)} satisfies:
#'
#' \deqn{
#' (-\log x)^{\gamma} x^{-\lambda} + \frac{1}{\alpha}\log(p)
#' = 0,
#' }
#'
#' where \eqn{\log(\cdot)} denotes the natural logarithm and
#' \eqn{\alpha}, \eqn{\gamma}, \eqn{\lambda > 0} are real-valued parameters.
#'
#' The quantiles are obtained numerically using
#' \code{\link[stats]{uniroot}}.
#'
#' @param p Numeric vector of probabilities with values in the interval
#'   \eqn{(0,1)}.
#' @param theta Numeric vector of length 3 containing the parameters
#'   \eqn{(\alpha,\gamma,\lambda)}, all positive.
#'
#' @return A numeric vector of the same length as \code{p} containing
#'   the quantiles of the UMW distribution.
#'
#' @examples
#' library(UMW)
#'
#' theta <- c(alpha = 0.7, gamma = 1.3, lambda = 0.5)
#' p <- c(0.1,0.5,0.9)
#' qUMW(p,theta)
#'
#' @export
qUMW <- function(p, theta) {
  if (!is.numeric(p) || !all(p > 0 & p < 1))
    stop("`p` must be numeric values in the interval (0,1).")
  if (!is.numeric(theta) || length(theta) != 3 || !all(theta > 0))
    stop("`theta` must be a numeric vector of length 3 with positive values.")
  alpha <- theta[1]
  gamma <- theta[2]
  lambda <- theta[3]
  sapply(p, function(pp) {
    eq <- function(x) {
      (-log(x))^gamma * x^(-lambda) +
        (1 / alpha) * log(pp)
    }
    stats::uniroot(eq, lower = 1e-4, upper = 1 - 1e-4)$root
  })
}

# Order statistics =============================================================

#' Density of the r-th Order Statistic of the Unit-Modified Weibull(UMW) Distribution
#'
#' Computes the probability density function of the r-th order statistic
#' from a random sample of size \eqn{n} drawn from the UMW distribution.
#'
#' The probability density function of the \eqn{r}-th order statistic is given by
#'
#' \deqn{
#' f_{Y_{(r)}}(x) =
#' \frac{n!}{(r-1)!(n-r)!}
#' \frac{\alpha (-\log x)^{\gamma}}{\log(x)\, x^{\lambda+1}}
#' (\lambda \log x - \gamma)
#' \left[1 - \exp\{-\alpha (-\log x)^{\gamma} x^{-\lambda}\}\right]^{n-r}
#' \exp\left\{- r \alpha (-\log x)^{\gamma} x^{-\lambda}\right\},
#' \quad 0 < x < 1,
#' }
#'
#' where \eqn{\log(\cdot)} denotes the natural logarithm and
#' \eqn{\alpha}, \eqn{\gamma}, \eqn{\lambda > 0}
#' are real-valued parameters.
#'
#' @param x Numeric vector with values in the unit interval \eqn{(0,1)}.
#' @param theta Numeric vector of length 3 containing the parameters
#'   \eqn{(\alpha,\gamma,\lambda)}.
#' @param r Integer indicating the order of the statistic
#'   \eqn{(1 \le r \le n)}.
#' @param n Integer indicating the sample size.
#'
#' @return A numeric vector containing the values of the probability density
#'   function of the \eqn{r}-th order statistic evaluated at \code{x}.
#'
#' @examples
#' library(UMW)
#'
#' theta <- c(alpha = 0.7, gamma = 1.3, lambda = 0.5)
#' n <- 10
#' r <- 3
#' x <- runif(n = 50,min = 0.001,max = 0.999)
#' dUMW_order(x, theta, r, n)
#'
#' @export
dUMW_order <- function(x, theta, r, n) {
  if (!is.numeric(x) || !all(x > 0 & x < 1))
    stop("`x` must be numeric values in the interval (0,1).")
  if (!is.numeric(theta) || length(theta) != 3 || !all(theta > 0))
    stop("`theta` must be a numeric vector of length 3 with positive values.")
  if (!is.numeric(r) || length(r) != 1 || r < 1 || r > n || r != floor(r))
    stop("`r` must be a single integer between 1 and n.")
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n))
    stop("`n` must be a single positive integer.")
  alpha  <- theta[1]
  gamma  <- theta[2]
  lambda <- theta[3]
  order_r<-(factorial(n) / (factorial(r - 1) * factorial(n - r))) *
    (alpha * (-log(x))^gamma / (log(x) * x^(lambda + 1)) *
       (lambda * log(x) - gamma) * (1 - exp(-alpha * (-log(x))^gamma *
       x^(-lambda)))^(n - r) * exp(-r * alpha * (-log(x))^gamma * x^(-lambda)))
  return(order_r)
}


# Hazard function ==============================================================

#' Hazard Function of the Unit-Modified Weibull(UMW) Distribution
#'
#' Computes the hazard rate function of the UMW distribution.
#'
#' The hazard function of the UMW distribution is defined as
#'
#' \deqn{
#' h(x) =
#' \frac{
#' \alpha\, x^{-(\lambda+1)} (-\log x)^{\gamma}
#' (\lambda \log x - \gamma)
#' }{
#' \log(x)\left[\exp\left(\alpha (-\log x)^{\gamma} x^{-\lambda}\right)-1\right]
#' },
#' \quad 0 < x < 1,
#' }
#'
#' where \eqn{\log(\cdot)} denotes the natural logarithm and
#' \eqn{\alpha}, \eqn{\gamma}, \eqn{\lambda > 0}
#' are real-valued parameters.
#'
#' @param x Numeric vector with values in the unit interval \eqn{(0, 1)}.
#' @param theta Numeric vector of length 3 containing the parameters
#'   \eqn{(\alpha,\gamma,\lambda)}.
#'
#' @return A numeric vector containing the values of the hazard function
#'   evaluated at \code{x}.
#'
#' @examples
#' library(UMW)
#'
#' x <- runif(n = 50,min = 0.001,max = 0.999)
#' theta <- c(alpha = 0.7, gamma = 1.3, lambda = 0.5)
#' hUMW(x, theta)
#'
#' @export
hUMW <- function(x, theta) {
  if (!is.numeric(x) || !all(x > 0 & x < 1))
    stop("`x` must be numeric values in the interval (0,1).")
  if (!is.numeric(theta) || length(theta) != 3 || !all(theta > 0))
    stop("`theta` must be a numeric vector of length 3 with positive values.")
  alpha  <- theta[1]
  gamma  <- theta[2]
  lambda <- theta[3]
  h <- (alpha * x^(-lambda - 1) * (-log(x))^gamma * (lambda * log(x) - gamma)) /
    (log(x) * (exp(alpha * (-log(x))^gamma / x^lambda) - 1))
  return(h)
}




