
# Fit UMW ======================================================================

#' Fit the Unit-Modified Weibull(UMW) Distribution
#'
#' Fits the UMW distribution to data via maximum likelihood estimation
#' and provides summary statistics, coefficient estimates, and
#' goodness-of-fit measures. The score vector
#' and the Hessian matrix are computed analytically.
#'
#' @param y Numeric vector of observed data. Values must lie in the open
#'   interval \eqn{(0,1)}. Missing values and observations outside this
#'   interval are removed prior to estimation.
#' @param method Character string specifying the optimization method to
#'   be passed to \code{\link[stats]{optim}} (default is \code{"BFGS"}).
#' @param start.theta Numeric vector of initial values for the parameters
#'   of the UMW distribution. Default is \code{c(1, 1)}.
#' @param print Logical; if \code{TRUE}, prints the results of the analysis
#'   to the console. If \code{FALSE}, no output is printed.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \strong{summary}: Named numeric vector with descriptive statistics
#'   and normality diagnostics of the sample, including minimum (Min.),
#'   median (Median), mean (Mean), maximum (Max.), standard deviation (sd),
#'   skewness (AC), kurtosis (K), and Anderson–Darling normality test p-value (p.AD).
#'   \item \strong{coef}: Matrix containing parameter estimates, standard errors
#'   and significance tests.
#'   \item \strong{metrics}: Numeric vector with information criteria and
#'   goodness-of-fit statistics: Akaike Information Criterion (AIC), Bayesian
#'   Information Criterion (BIC), corrected AIC (AICc), Kolmogorov–Smirnov statistic (KS),
#'   Anderson–Darling statistic (AD), and Cramér–von Mises statistic (CvM).
#'   \item \strong{convergence}: Logical indicating whether the optimization algorithm
#'   converged successfully.
#'   \item \strong{num_na}: Number of missing observations in the original data.
#'   \item \strong{num_omit}: Number of observations omitted for being outside the
#'   interval \eqn{(0,1)}.
#'   \item{n}{Number of effective observations in the sample (after removing missing
#'   values and values outside the interval (0,1)).}
#'   \item{method}{Optimization method used during parameter estimation.}
#' }
#'
#' @examples
#' library(UMW)
#'
#' x<-runif(50,0,1)
#' fit_UMW(x,method="BFGS")
#'
#' @export
fit_UMW<-function(x,method="BFGS",start.theta=c(1,1),print=T)
{
  out<-c()
  n_x <- length(x)
  if(any(is.na(x))) {
    num_na <- sum(is.na(x))
    out$num_na <- num_na
  }else {
    out$num_na <- 0
  }
  if(any(x <= 0) | any(x >= 1) | any(is.na(x))){
    x <- clean_vector(vector = x)
    n_x2 <- length(x)
    num_omit <- n_x - n_x2 - out$num_na
    out$num_omit <- num_omit
  }else {
    out$num_omit <- 0
  }
  out$n <- length(x)
  if (out$n <= 3) {
    stop("The sample must be larger and/or all values must lie within the interval (0,1).")
  }

  out$summary <- c(round(c(summary(x)[c(1, 3, 4, 6)],
                   sd = sd(x),AC = moments::skewness(x), K = moments::kurtosis(x),
                                         p.AD = nortest::ad.test(x)$p.value), 3))
  out$method<-method
  mod1<-suppressWarnings(try(Est_UMW(x=x,method = method,applic = TRUE),T))
  if(is.null(mod1)){stop("ALGORITHM DID NOT CONVERGE!")}
  out$metrics<-metrics(y=x,loglik = mod1$value,par = mod1$par,F_dist = pUMW,k=3)
  out$coef<-Coef_estim(par = mod1$par,hessian = mod1$hessian,
                       name_par = c("α", "γ", "λ"))
  out$convergence<-mod1$convergence==0
  #
  if(print==T){
  cat("
n:",out$n,"   NA's:",out$num_na,"   Omit(0 \u2265 x \u2265 1):",out$num_omit,"
")
  cat("
Summary:
")
  print(out$summary)
  cat("
Coefficients:
")
  print(out$coef)
  cat("
Metrics:
")
  tab<-round(as.numeric(out$metrics),3)
  names(tab)<-c("AIC","BIC","AICc","KS","AD","CvM")
  print(tab)
  }
  invisible(out)
}


## Cleaning function (0.1) -----------------------------------------------------

clean_vector <- function(vector)
{
  vector[vector <= 0] <- NA
  vector[vector >= 1] <- NA
  vector <- na.omit(vector)
  return(as.numeric(vector))
}

## Metrics function ------------------------------------------------------------

metrics<-function(y,loglik,par,F_dist,k)
{
  KS <- function(theta,y,n,fda){
    y <- sort(y)
    aux1 <- c()
    for(i in 1:n){
      aux1[i] <- abs(fda(theta=theta,x=y[i])-i/n)
    }
    return(max(aux1))
  }
  AD <- function(theta,y,n,fda){
    x1 = sort(y)
    u = fda(theta=theta,x=x1)
    aux=rep(0,n)
    for(i in 1:n){
      aux[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
    }
    A2 = -n -(1/n)*sum(aux)
    return(A2)
  }
  cramer <- function(theta,y,n,fda){
    x1 = sort(y)
    u = fda(theta=theta,x=x1)
    aux=rep(0,n)
    for(i in 1:n){
      aux[i]=(((2*i-1)/(2*n))-u[i])^2
    }
    W2=sum(aux)+1/(12*n)
    return(W2)
  }
  n<-length(y)
  AIC1 <- -2*loglik+2*k
  BIC1 <- -2*loglik+log(n)*k
  AICc1 <- AIC1 + ((2*k^2)+(2*k))/(n-k-1)
  KS1 <- KS(theta = par,y = y,n = n,fda = F_dist)
  AD1 <- AD(theta = par,y = y,n = n,fda = F_dist)
  CvM1 <- cramer(theta = par,y = y,n = n,fda = F_dist)
  return(list(AIC = AIC1, BIC = BIC1, AICc = AICc1, KS = KS1, AD = AD1, CvM = CvM1))
}

## coef ------------------------------------------------------------------------

#' Compute Coefficient Table with Standard Errors, Z-Values, and P-Values
#'
#' This function computes the standard errors, z-values, and p-values for a
#' set of parameter estimates given a Hessian matrix or a log-likelihood function.
#'
#' @param par Numeric vector of parameter estimates.
#' @param hessian Optional numeric matrix. The Hessian evaluated at \code{par}.
#'   If \code{NULL}, the Hessian is computed numerically from \code{f_dist} and \code{y}.
#' @param f_dist Optional function. Probability density function used to compute
#'   the log-likelihood. **The function must have arguments `theta` (parameters) and `x` (data).**
#' @param y Optional numeric vector. Observed data for computing the Hessian.
#' @param name_par Optional character vector. Names of the parameters.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item Estimate: estimated parameter
#'   \item Std. Error: standard error
#'   \item z value: z-statistic
#'   \item Pr(>|z|): p-value
#'   \item Significance: significance codes
#' }
#'
#' @examples
#' library(UMW)
#'
#' par <- c(0.7, 1.3, 0.5)
#' name_par <- c("alpha", "gamma", "lambda")
#'
#' # Example using the Hessian matrix:
#'
#' hessian_ex <- -matrix(c(10, 2, 1, 2, 8, 0.5, 1, 0.5, 5),
#'                       nrow = 3, byrow = TRUE) # fictitious Hessian
#'
#' Coef_estim(par = par, hessian = hessian_ex, name_par = name_par)
#'
#' # Example using the probability density function:
#'
#' set.seed(1)
#' y <- runif(50, min = 0.001, max = 0.999)
#'
#' \dontrun{
#' dUMW1<-function(x,theta) {
#'   alpha<-theta[1]
#'   gamma<-theta[2]
#'   lambda<-theta[3]
#'   (alpha/(log(x)*x^(lambda+1)))*((-log(x))^gamma)*(lambda*log(x)-gamma)*
#'     exp(-alpha*((-log(x))^gamma)*x^(-lambda))
#' }
#' }
#' Coef_estim(par = par, f_dist = dUMW1, y = y, name_par = name_par)
#'
#' @export
Coef_estim <- function(par, hessian = NULL, f_dist = NULL, y = NULL, name_par = NULL) {
  if (is.null(hessian)) {
    if (is.null(f_dist) | is.null(y)) {
      stop("Either 'hessian' must be provided, or both 'f_dist' and 'y' must be non-NULL.")
    }
    args_fd <- names(formals(f_dist))
    if (!all(c("theta", "x") %in% args_fd)) {
      stop("f_dist must have arguments named 'theta' and 'x'.")
    }
    hessian_dist <- function(par, y, f_dist) {
      loglik_func <- function(par, y) {
        sum(log(f_dist(theta = par, x = y)))
      }
      numDeriv::hessian(func = loglik_func, x = par, y = y, method = "Richardson")
    }
    hessian <- hessian_dist(par = par, y = y, f_dist = f_dist)
  }
  stderror<-sqrt(diag(solve(-hessian)))
  #-----------------------z-Values and P-values---------------------------------
  thetaH0<-as.vector(c(rep(0,length(par))))
  zstat<-c(abs((par-thetaH0)/stderror))
  pvalues<-2*(1-pnorm(zstat))
  #-----------------------Result Table------------------------------------------
  options(scipen = 1)
  model<-as.data.frame(cbind(round(par,3),round(stderror,3),round(zstat,3),round(pvalues,3)))
  model<-model |> dplyr::mutate(
    nivsig=dplyr::case_when(
      V4 <= 0.001 ~ "***",
      V4 > 0.001 & V4 <= 0.01 ~ "**",
      V4 > 0.01 & V4 <= 0.05 ~ "*",
      V4 > 0.05 & V4 <= 0.1 ~ "·",
      TRUE ~ ""),
    V4=dplyr::case_when(
      V4 < 0.001 ~ paste0("<",0.001),
      TRUE ~ paste0(V4)),
    V2=dplyr::case_when(
      V2 < 0.001 ~ paste0("<",0.001),
      TRUE ~ paste0(V2)))
  colnames(model)<-c("Estimate","Std. Error","z value","Pr(>|z|)","")
  if(!is.null(name_par)){row.names(model)<-name_par}
  return(model)
}

# Fit RQUMW ====================================================================

#' Fit the Regression Quantile Unit-Modified Weibull (RQ-UMW) Model
#'
#' Fits the RQ-UMW model by numerically maximizing the
#' log-likelihood function. The score vector and Hessian matrix are computed
#' analytically.
#'
#' @param f A symbolic formula describing the model structure. The response
#'   variable must lie in the open interval \eqn{(0,1)}. Additional model
#'   components may be specified using the \code{|} operator.
#' @param data Optional data frame containing the variables in the model.
#' @param tau Quantile level in \eqn{(0,1)} to be modeled. Default is \code{0.5}.
#' @param link_mu Character string or function specifying the link for the
#'   \eqn{\mu} parameter. Can be "logit", "probit", "cauchit", "cloglog", "loglog",
#'   or a user-supplied function implementing the desired link. Default is \code{"logit"}.
#' @param method Optimization method used by \code{\link[stats]{optim}}.
#'   Possible values are \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"} and \code{"SANN"}.
#' @param start.theta Optional numeric vector of initial values for the
#'   optimization algorithm.
#' @param printmodel Logical indicating whether a model summary should be
#'   printed to the console. Default is \code{TRUE}.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{formula}: Model formula supplied by the user.
#'   \item \strong{data}: Data frame used for fitting the model.
#'   \item \strong{y}: Numeric vector of observed responses.
#'   \item \strong{X, W, Z}: Design matrices for the \eqn{\mu}, \eqn{\gamma},
#'   and \eqn{\lambda} parameters, respectively. Each column represents a covariate.
#'   \item \strong{pars}: Estimated parameter vector.
#'   \item \strong{coef}: Matrix with parameter estimates, standard errors,
#'   z-statistics, p-values, and significance codes.
#'   \item \strong{fitted}: Vector of fitted values.
#'   \item \strong{residuals}: Quantile residuals.
#'   \item \strong{loglik}: Log-likelihood of the fitted model.
#'   \item \strong{metrics}: List containing information criteria (AIC, BIC, AICc),
#'   the generalized \eqn{R^2}, and the Anderson–Darling test for normality of
#'   the quantile residuals.
#'   \item \strong{convergence}: Logical indicating whether the optimization
#'   algorithm converged successfully.
#'   \item \strong{n}: Number of observations.
#'   \item \strong{method}: Optimization method used.
#'   \item \strong{alpha}: Estimated \eqn{\alpha} parameter of the UMW distribution.
#'   \item \strong{quantile}: Quantile level \eqn{\tau} used in the model.
#'   \item \strong{link}: List of all link functions used in the model:
#'         \code{g_mu, g_gamma, g_lambda} and their inverses
#'         \code{ginv_mu, ginv_gamma, ginv_lambda}.
#'   \item \strong{loglik0}: Log-likelihood of the null model (without regression structure).
#'   \item \strong{outer.iter}: Number of iterations performed by the optimizer.
#'   \item \strong{df.residual}: Residual degrees of freedom.
#'   \item \strong{st.res}: Standardized residual variance.
#'   \item \strong{zeta_hat}: Linear predictor for the \eqn{\mu} parameter (\eqn{X \beta}).
#'   \item \strong{hessian}: Hessian matrix of the log-likelihood at the estimated parameters.
#'   \item \strong{measures}: Numeric vector of error measures computed on the fitted values:
#'         Mean Squared Error (MSE), Root Mean Squared Error (RMSE),
#'         Mean Absolute Error (MAE), and Mean Absolute Percentage Error (MAPE).
#' }
#'
#' @examples
#' library(UMW)
#'
#' set.seed(34)
#' y<-runif(100, 0.01, 0.99)
#' X1<-rnorm(100)
#' X2<-rnorm(100)
#' data <- data.frame(y = y, X1 = X1, X2 = X2)
#'
#' fit1 <- fit_RQUMW(f = y ~ . | 1 | 1, data = data, tau = 0.5,
#'                   method = "BFGS",link_mu = "logit")
#'
#' fit2 <- fit_RQUMW(y ~ ., data = data, tau = 0.25,
#'                   method = "CG",link_mu = "probit")
#'
#' fit3 <- fit_RQUMW(y ~ X1 + X2 -1, data = NULL, tau = 0.75,
#'                   method = "Nelder-Mead",link_mu = "cloglog")
#'
#' fit4 <- fit_RQUMW(f = y ~ X1-1| 1 | X1, data = NULL, tau = 0.5,
#'                   method = "BFGS",link_mu = "loglog")
#'
#' fit5 <- fit_RQUMW(f = y ~ X1 |X2-1|X1, data = NULL, tau = 0.35,
#'                   method = "BFGS",link_mu = "cauchit")
#'
#' @export

fit_RQUMW<-function(f,data = NULL,tau=0.5,link_mu="logit",method="BFGS",
                          start.theta=NULL,printmodel=T)
{
  #----------------------------------------------------------------------------#
  if (!is.numeric(tau) || !all(tau > 0 & tau < 1))
    stop("`tau` must be numeric values in the interval (0,1).")
  if (!(is.null(data) || is.data.frame(data))) {
    stop("'data' must be either NULL or a data.frame.")
  }
  f <- as.formula(f, env = parent.frame())
  parts <- strsplit(deparse(f), "\\|")[[1]]
  parts <- trimws(parts)
  if (length(parts) == 1) parts <- c(parts, "1", "1")
  if (length(parts) == 2) parts <- c(parts, "1")
  fix_part <- function(p) {
    p <- trimws(p)
    if (grepl("y~", p, fixed = TRUE)) {
      return(p)
    }
    return(paste("y~", p))
  }
  f_mu <- as.formula(parts[1])
  y <- as.numeric(model.response(model.frame(f_mu, data = data)))
  X <- model.matrix(f_mu, data = data)
  W <- model.matrix(as.formula(fix_part(parts[2])), data = data)
  Z <- model.matrix(as.formula(fix_part(parts[3])), data = data)
  if (!is.numeric(y) || !all(y > 0 & y < 1))
    stop("`y` must be numeric values in the interval (0,1).")
  #-----------------------Output-----------------------------------------------#
  out<-c()
  if(is.null(data)){
    cbdata<-as.data.frame(cbind(y,X,W,Z))
    out$data<-cbdata[, sapply(cbdata, function(x) length(unique(x)) > 1) &
                       !duplicated(as.list(cbdata))]
  }else{out$data<-data}
  out$formula<-f;out$X<-as.matrix(X); out$W<-as.matrix(W); out$Z<-as.matrix(Z);
  out$y<-y; out$quantile<-tau; out$n<-length(y); out$method<-method
  if(class(link_mu)=="function"){out$link$link_mu<-"unknown"}else{out$link$link_mu<-link_mu}
  class(out)<-c("RQ-UMW")
  #-----------------------Choice of Link Function------------------------------#
  out$link<-modifyList(out$link, get_link_functions(link_mu = link_mu))
  out$link<-modifyList(out$link, func_linkWZ(W=W,Z=Z))
  out$link[c("n_W","n_Z","Z","W")] <- NULL
  #-----------------------colnames(X)------------------------------------------#
  beta0<-c(expression(beta[0]),expression(beta[1]),expression(beta[2]),expression(beta[3]),
           expression(beta[4]),expression(beta[5]),expression(beta[6]),expression(beta[7]))
  if(!is.null(colnames(X))){tmpnames1<-colnames(X)}else{tmpnames1<-c(beta0[1:ncol(as.matrix(X))])}
  if(!is.null(colnames(W))){tmpnames2<-colnames(W)}else{tmpnames2<-c(beta0[1:ncol(as.matrix(W))])}
  if(!is.null(colnames(Z))){tmpnames3<-colnames(Z)}else{tmpnames3<-c(beta0[1:ncol(as.matrix(Z))])}
  #-------------------------Estimate-------------------------------------------#
  tmp1<-suppressWarnings(EST_RQUMW(y = y, X = out$X, W = out$W, Z = out$Z,
        method = method,tau = tau,g_mu = out$link$g_mu,ginv_mu = out$link$ginv_mu,
        g_gamma = out$link$g_gamma,ginv_gamma = out$link$ginv_gamma,
        g_lambda = out$link$g_lambda,ginv_lambda = out$link$ginv_lambda,
        start.theta = start.theta,applic = T))
  if (class(tmp1) != "list") {stop("ALGORITHM DID NOT CONVERGE!")}
  out$loglik<-tmp1$value
  out$convergence<-tmp1$convergence==0
  out$pars<-tmp1$par
  out$outer.iter<-as.numeric(tmp1$counts[1])
  out$hessian<- tmp1$hessian
  #-----------------------Initial Part of Function Output----------------------#
  beta<-out$pars[1:ncol(X)]
  out$zetahat<-as.vector(X%*%as.matrix(beta))
  gamma<-W%*%as.matrix(out$pars[(ncol(X)+1):(ncol(X)+ncol(W))])
  lambda<-Z%*%as.matrix(out$pars[(ncol(X)+ncol(W)+1):(ncol(X)+ncol(W)+ncol(Z))])
  out$fitted<-as.vector(out$link$ginv_mu(out$zetahat))
  names(out$pars)<-c(paste0("μ:",tmpnames1),paste0("γ:",tmpnames2),paste0("λ:",tmpnames3))
  out$qrnames<-names(out$pars)
  out$alpha<- as.vector(-((out$fitted^lambda)*log(tau))/((-log(out$fitted)))^(gamma))
  #-------------------------Estimate0------------------------------------------#
  tmp0<-suppressWarnings(Est_UMW(x=y,method=method,applic = T))
  out$loglik0<-tmp0$value
  #-----------------Cumulative Distribution Function---------------------------#
  peged<-function(y,alpha,gamma,lambda){(exp(-(alpha*(-log(y))^gamma)/y^lambda))}
  #------Quantile Residue------------------------------------------------------#
  k_par <- length(out$pars)
  out$residuals<-as.numeric(qnorm(peged(y=y,alpha=out$alpha,gamma=gamma,lambda=lambda)))
  if (any(is.infinite(out$residuals), na.rm = TRUE)) {
    stop("Residuals contain infinite values.")
  }
  out$df.residual<-out$n-k_par
  out$st.res<-sum((as.vector(out$residuals)-mean(as.vector(out$residuals)))^2)/out$df.residual
  out$coef<-Coef_estim(par = out$pars,hessian = out$hessian,name_par = out$qrnames)
  #--------------------Model Selection Criteria & R^2--------------------------#
  out$metrics$aic <- -2*out$loglik+2*(k_par)
  out$metrics$bic <- -2*out$loglik+log(out$n)*(k_par)
  out$metrics$aicc <- out$metrics$aic+(((2*(k_par)^2)+(2*(k_par)))/(out$n-(k_par)-1))
  out$metrics$tAD_resid <- nortest::ad.test(out$residuals)
  out$metrics$R2G <- pmax(0,R2G(loglik = out$loglik,loglik0 = out$loglik0,n = out$n))
  #----------------------------Error measures----------------------------------#
  out$measures <- error_measures(actual = out$y,predicted = out$fitted)
  #----------------------------Result Table------------------------------------#
  if(printmodel==TRUE){
    summary_RQUMW(out)
  }
  invisible(out)
}


## Summary RQ-UMW --------------------------------------------------------------

#' Summary method for RQ-UMW models
#'
#' Produces a formatted summary of an object of class \code{"RQ-UMW"},
#' including parameter estimates, standard errors, goodness-of-fit
#' statistics, and residual diagnostics.
#'
#' @param fit An object of class \code{"RQ-UMW"} returned by
#'   \code{\link{fit_RQUMW}}.
#'
#' @return
#' Prints a summary of the fitted model to the console.
#'
#' @export
summary_RQUMW<-function(fit)
{
  if (!inherits(fit, "RQ-UMW")) {
    stop("Object must be of class 'RQ-UMW'.")
  }
  cat(c("============================================================="),fill = TRUE)
  cat(c("Link:",fit$link$link_mu," ","Quantile:",fit$quantile," ","Optimization:",fit$method),fill = TRUE)
  cat(c("============================================================="),fill = TRUE)
  print(fit$coef)
  cat(c("---"),fill = TRUE)
  cat(c("Significance code: 0 *** 0.001 ** 0.01 * 0.05 · 0.1 1"),fill = TRUE)
  cat(c("-------------------------------------------------------------"),fill = TRUE)
  cat(c("Loglike:",round(fit$loglik,3)," ","AD(p) Resid.:",round(fit$metrics$tAD_resid$p.value,3)," ","R2G:",round(fit$metrics$R2G,3)),fill=TRUE)
  cat(c("AIC:",round(fit$metrics$aic,3)," ","BIC:",round(fit$metrics$bic,3)," ","AICc:",round(fit$metrics$aicc,3)),fill=TRUE)
  cat(c("MSE:",round(fit$measures$MSE,3)," ","RMSE:",round(fit$measures$RMSE,3)," ","MAE:",round(fit$measures$MAE,3))," ","MAPE:",round(fit$measures$MAPE,3),fill=TRUE)
  cat(c("---"),fill = TRUE)
  cat(c("Number of Iterations in the Optimization Algorithm:",fit$outer.iter),fill=TRUE)
  cat((c("Residual standard error: ",round(fit$st.res,3)," on ",fit$df.residual," degrees of freedom")),fill=TRUE)
  cat("Residuals:",fill=TRUE)
  print(round(summary(as.vector(fit$residuals)),3))
  cat(c("============================================================="),fill = TRUE)
}

## Error measures --------------------------------------------------------------

error_measures <- function(actual, predicted)
{
  mse <- suppressWarnings(mean((actual - predicted)^2))
  rmse <- suppressWarnings(sqrt(mean((actual - predicted)^2)))
  mae <- suppressWarnings(mean(abs(actual - predicted)))
  mape <- suppressWarnings(mean(abs((actual - predicted) / actual)) * 100)
  results <- list(MSE = mse,RMSE = rmse,MAE = mae,MAPE = mape)
  return(results)
}

## R2 generalized coefficient of determination (nagelkerke) --------------------

R2G<-function(loglik,loglik0,n)
{
  R2<- 1-exp((-2/n)*((as.numeric((loglik)))-(as.numeric(loglik0))))
  return(c(R2G=R2))
}
