# Generate reparameterized sample RQ-UMW =======================================

rRQUMW<-function(n,theta,X,Z,tau=0.5,ginv_mu,ginv_lambda)
{
  n_X <- ncol(X); n_Z <- ncol(Z)
  #
  beta_mu <- theta[1:n_X]
  gamma_i <- theta[(n_X+1)]
  beta_lambda <- theta[(n_X+2):(n_X+1+n_Z)]
  #
  mu_i    <- ginv_mu(as.vector(X %*% beta_mu))
  lambda_i<- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i)*log(tau))/((-log(mu_i))^(gamma_i))
  #
  eq<-function(x,parms)
  {((-log(x))^parms[3])*(x^(-parms[4]))+(1/parms[2])*log(parms[1])}
  y<-c()
  j<-1
  while(j<=n){
    u<-runif(1)
    parms<-c(u,alpha_i[j],gamma_i,lambda_i[j])
    tmp<-try(suppressWarnings(stats::uniroot(f=eq,lower=0.0001,upper=0.9999,parms=parms)),T)
    if(class(tmp)=="list"){
      y[j]<-tmp$root
      j<-j+1
    }
  }
  return(y)
}

# Other Functions ==============================================================

## Function to save simulation -------------------------------------------------

save_cen_RQ2<-function(outputRQUMW,RF,tau,n,write="sim")
{
  l1<-nrow(outputRQUMW$sim$output)
  if(l1>=RF){
    save(outputRQUMW,file=paste0(write,"_R",RF,"_n",n,"_t",tau,".RData"))
  }
  if(l1<RF){
    save(outputRQUMW,file=paste0(write,"_R",l1,"_n",n,"_t",tau,".RData"))
  }
}

## function for connection -----------------------------------------------------

func_linkWZ<-function(Z)
{
  out<-c()
  if(is.null(Z)){out$Z<-matrix(1,n,1)}else{out$Z=as.matrix(Z)}
  out$n_Z <- ncol(Z)
  if(all(out$Z[,1] == 1, na.rm = TRUE) & out$n_Z==1){
    out$g_lambda <- out$ginv_lambda <- function(x) x
    out$link_lambda <- "identity"
  }else{
    out$g_lambda    <- function(l)   log(l)
    out$ginv_lambda <- function(eta) exp(eta)
    out$link_lambda <- "log"
  }
  return(out)
}

## generate covariates ---------------------------------------------------------

gen_covariates <- function(formula, n)
{
  p <- strsplit(deparse(formula), "\\|")[[1]]
  p <- c(trimws(p), "1")[1:2]
  get_info <- function(f) {
    tt <- terms(as.formula(f))
    list(intercept = attr(tt, "intercept"),
      k = length(attr(tt, "term.labels")))
  }
  info <- list(mu     = get_info(p[1]),
    lambda = get_info(paste("~", p[2])))
  gen <- function(info){
    X <- if (info$k > 0)
      matrix(runif(n * info$k), n) else NULL
    if (info$intercept)
      X <- cbind(rep(1,n), X)
    X
  }
  list(X = gen(info$mu),Z = gen(info$lambda))
}


## mu linking function ---------------------------------------------------------

get_link_functions <- function(link_mu)
{
  finv_num <- function(link_function, eta_value, start_value = rep(0.5,length(eta_value))) {
    inverse_function <- function(mu) {
      link_function(mu) - eta_value
    }
    result <- suppressWarnings(nleqslv::nleqslv(start_value, inverse_function)$x)
    return(result)
  }
  if(class(link_mu)=="character"){
    if(link_mu == "logit"){
      g_mu <- function(mu) {log(mu / (1 - mu))}
      ginv_mu <- function(eta) {1 / (1 + exp(-eta))}
    }
    if(link_mu=="probit"){
      g_mu <- function(mu) {qnorm(mu)}
      ginv_mu <- function(eta) {pnorm(eta)}
    }
    if(link_mu=="cauchit"){
      g_mu <- function(mu) {tan(pi * (mu - 0.5))}
      ginv_mu <- function(eta) {0.5 + (atan(eta) / pi)}
    }
    if(link_mu=="cloglog"){
      g_mu <- function(mu) {log(-log(1 - mu))}
      ginv_mu <- function(eta) { 1 - exp(-exp(eta))}
    }
    if(link_mu=="loglog"){
      g_mu <- function(mu) {log(-log(mu))}
      ginv_mu <- function(eta) {exp(-exp(eta))}
    }
  }
  if(class(link_mu)=="function"){
    g_mu <- link_mu
    ginv_mu <- function(eta) {finv_num(g_mu,eta_value = eta)}
  }
  if(class(link_mu)!="function" & class(link_mu)!="character"){
    stop("Return the link_mu function as a 'function' class or using the definitions
         'logit', 'probit', 'cauchit', 'loglog', 'cloglog'.")
  }
  return(list(g_mu=g_mu,ginv_mu=ginv_mu))
}


## RQ-UMW function sample replicates -------------------------------------------

samples_RQUMW <- function(theta, X, Z, n = 50, re = 100, tau = 0.5,
                          ginv_mu, ginv_lambda,set_seed = NULL,
                          n_cores = (parallel::detectCores() - 1))
{
  require(foreach)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  if (!is.null(set_seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(set_seed)
    parallel::clusterSetRNGStream(cl, iseed = set_seed)
    doRNG::registerDoRNG(set_seed)
  }
  opts <- progresso(iterations = re, sec = "[1/2]")
  samples_list <- try(foreach(j = 1:re,.packages = c("foreach"),.options.snow = opts,
                          .export = c("rRQUMW")) %dopar% {
                            as.numeric(rRQUMW(n = n, theta = theta, X = X, Z = Z,
                               tau = tau, ginv_mu = ginv_mu,ginv_lambda = ginv_lambda))
  },T)
  foreach::registerDoSEQ()
  parallel::stopCluster(cl)
  amostra <- as.data.frame(do.call(cbind, samples_list))
  return(amostra)
}


## RQ-UMW monte carlo simulation function --------------------------------------

sim_est_RQUMW <- function(sample, theta, X, Z, n, tau,
                          re = 100, RF = 100, method = "BFGS", g_mu,
                          ginv_mu, g_lambda, ginv_lambda, start.theta = NULL,
                          n_cores = (parallel::detectCores() - 1))
{
  require(foreach)
  opts <- progresso(iterations = re, sec = "[2/2]")
  len_k <- length(theta)
  estimate<-bigstatsr::FBM(nrow=re,ncol=len_k)
  n_NA<-bigstatsr::FBM(nrow=re,ncol=1)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  time <- system.time({try(foreach(k = 1:re,.packages = c("foreach"),
      .options.snow = opts,.export = c("EST_RQUMW", "llike_RQUMW", "is.positive",
                  "vscore_RQUMW", "hessian_RQUMW","test.fun", "which.NA")) %dopar% {
    output<-as.data.frame(estimate[])
    output[output==0]<-NA
    l_out<-nrow(which.NA(output))
    if(l_out<RF){
      estimate[k,]<-EST<- suppressWarnings(EST_RQUMW(y = sample[, k], X = X, Z = Z,
      tau = tau,g_mu = g_mu, method = method,ginv_mu = ginv_mu,
      g_lambda = g_lambda, ginv_lambda = ginv_lambda,
      start.theta = start.theta, applic = FALSE))
      if(any(is.na(EST))){n_NA[k]<-1}
    }
    if(l_out>=RF){
    }
  },T)})
  foreach::registerDoSEQ()
  parallel::stopCluster(cl)
  #
  n_NA_total<-sum(n_NA[])
  output1<-as.data.frame(estimate[])
  output1[output1==0]<-NA
  output1<-which.NA(output1)
  if (nrow(output1) > RF) output1 <- output1[1:RF, ]
  #
  tabela <- tab_est(output1, theta)
  row.names(tabela) <- c(param_names_RQUMW(X,Z))
  print(tabela)
  #
  sec <- round(time[["elapsed"]])
  time_str <- sprintf("%02dh %02dm %02ds",sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
  return(list(output = output1,table = tabela,no_converg = n_NA_total,time = time_str))
}

## parameter names -------------------------------------------------------------

param_names_RQUMW <- function(X, Z)
{
  has_intercept <- function(M) {ncol(M) > 0 && var(M[, 1]) == 0}
  idx <- function(M) {
    if (is.null(M) || ncol(M) == 0) return(character(0))
    p <- ncol(M)
    if (has_intercept(M))
      c("(icpt)", seq_len(p - 1))
    else
      seq_len(p)
  }
  c(paste0("β", idx(X)),"γ",paste0("δ", idx(Z)))
}



# Simulate RQ-UMW ==============================================================

#' Monte Carlo Simulation Study for the Regression Quantile Unit-Modified Weibull (RQ-UMW) Model
#'
#' Conducts a Monte Carlo simulation study for the (RQ-UMW) model under different
#' sample sizes and quantile levels.
#'
#' @param f A model formula that specifies the covariate structure of the RQUMW
#' model.
#' @param theta Numeric vector containing the true parameter values used
#'   in the data-generating process.
#' @param n Integer vector indicating the sample sizes to be considered.
#' @param re Integer indicating the number of Monte Carlo replications.
#' @param RF Integer number indicating the desired number of replicates.
#' @param tau Numeric vector with quantile levels in the interval \eqn{(0,1)}.
#' @param link_mu Character string or function specifying the link function
#'   for the location parameter. If a character string is provided, it must be
#'   one of \code{"logit"}, \code{"probit"}, \code{"cauchit"},
#'   \code{"loglog"}, or \code{"cloglog"}. Alternatively, a custom link
#'   function may be supplied as an object of class \code{"function"}.
#' @param save Logical; if \code{TRUE}, the simulation results are saved to disk.
#' @param n_cores Integer indicating the number of CPU cores to be used
#'   for parallel computation.
#' @param cen_name Character string used as a prefix for saved files.
#' @param method Optimization method used by \code{\link[stats]{optim}}.
#'   Possible values are \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"} and \code{"SANN"}.
#' @param start.theta Optional numeric vector of initial values for the
#' parameters. If \code{NULL}, default values are internally
#'   computed.
#' @param set_seed Integer. Optional random seed for reproducibility.
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{sample}: generated samples from the RQUMW model for each
#'     combination of \code{n} and \code{tau};
#'   \item \code{sim}: a list containing the results of the Monte Carlo study,
#'   with the following elements:
#'   \itemize{
#'     \item \code{output}: a numeric matrix of parameter estimates, where
#'     each row corresponds to one Monte Carlo replication and each column
#'     corresponds to one model parameter;
#'     \item \code{table}: a numeric matrix (or data frame) summarizing the
#'     Monte Carlo results, including the true parameter values, empirical
#'     mean, bias, mean squared error (MSE), asymmetry (AS) and kurtosis (K).
#'     This table is printed to the console as a side effect of the function;
#'     \item \code{no_converg}: an integer giving the number of Monte Carlo
#'     replications for which the estimation algorithm did not converge;
#'     \item \code{time}: a character string reporting the total computation
#'     time of the Monte Carlo study.
#'   \item \code{cen_cov}: generated covariates used in the simulations;
#'   \item \code{link_mu}: the link function specification used for the
#'     \eqn{\mu} parameter.
#' }}
#'
#' @examples
#' library(UMW)
#'
#' # Example without saving
#' theta1<-c(beta=c(0.2,-0.4),gamma=c(1.5),delta=c(0.8,1.4))
#' f1<-y~X|Z
#'
#' simulate_RQUMW(f = f1,theta = theta1,n = c(100),tau = c(0.5),
#'                re = 1100,RF = 1000,method = "BFGS",set_seed = 25)
#'
#' # Example saving
#' theta2<-c(beta=c(0.5,-0.6,0.2),gamma=c(1.5),delta=c(2.3))
#' f2<-y~X1+X2|1
#'
#' simulate_RQUMW(f = f2,theta = theta2,n = c(50,100),tau = c(0.25,0.5),re = 1100,
#'                RF = 1000,save = T,cen_name = paste0("sim_cen2"),
#'                method = "BFGS",set_seed = 25)
#'
#' @export
simulate_RQUMW <- function(f = y ~ X1 + X2|1,theta=c(0.5,-0.6,0.2,1.5,2.3),method="BFGS",
                           n=50,re=100,RF=100,tau=0.5,link_mu="probit",save=F,start.theta = NULL,
                           n_cores = (parallel::detectCores()-1),cen_name="sim",set_seed=NULL)
{
  if (!is.null(set_seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(set_seed)
  }
  glf<-get_link_functions(link_mu = link_mu)
  g_mu<-glf$g_mu
  ginv_mu<-glf$ginv_mu
  #
  nsamplerq <- sim_rumw <- cen_cov <- list()
  for (tauc in tau) {
    a1<-paste0("t", tauc)
    nsamplerq[[a1]] <- sim_rumw[[a1]] <- cen_cov[[a1]] <- list()
    for (nc in n){
      cat("
tau =",tauc,"& n =",nc,"----
     ")
      a2<-paste0("n", nc)
      cen_cov[[a1]][[a2]]=gen_covariates(formula = f,n = nc)
      X<-cen_cov[[a1]][[a2]]$X;Z<-cen_cov[[a1]][[a2]]$Z
      funcWZ<-func_linkWZ(Z=Z)
      #
      nsamplerq[[a1]][[a2]] <- samplerq <- suppressWarnings(
        samples_RQUMW(n=nc,theta=theta,X=X,Z=Z,ginv_mu=ginv_mu,tau=tauc,re=re,
                      ginv_lambda=funcWZ$ginv_lambda,n_cores = n_cores,set_seed=set_seed))
      sim_rumw[[a1]][[a2]] <- suppressWarnings(
        sim_est_RQUMW(sample = samplerq, X = X, Z = Z, n = nc,method = method,
                      theta = theta, tau = tauc, re = re, RF = RF,start.theta = start.theta,
                      n_cores =  n_cores,g_mu = g_mu,g_lambda = funcWZ$g_lambda,
                      ginv_mu = ginv_mu,ginv_lambda = funcWZ$ginv_lambda))
      # Save
      if(save==T){
        save_n<-list(sample = nsamplerq[[a1]][[a2]],sim = sim_rumw[[a1]][[a2]],
                     cen_cov = cen_cov[[a1]][[a2]])
        save_cen_RQ2(outputRQUMW = save_n, RF = RF,  n = nc, tau = tauc, write = paste0(cen_name,"_RQUMW"))
      }
    }
    cat("--------------------------------------------------------
")
  }
  invisible(list(sample = nsamplerq,sim = sim_rumw,cen_cov = cen_cov,link_mu=link_mu))
}


# theta1<-c(beta=c(0.2,-0.4),gamma=c(1.5),delta=c(0.8,1.4))
# f1<-y~X|Z
#
# simulate_RQUMW(f = f1,theta = theta1,n = c(100),tau = c(0.5),
#                re = 1100,RF = 1000,method = "BFGS",set_seed = 25)

