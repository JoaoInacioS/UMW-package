# Generate reparameterized sample RQ-UMW =======================================

rRQUMW<-function(n,theta,X,W,Z,tau=0.5,ginv_mu,ginv_gamma,ginv_lambda)
{
  n_X <- ncol(X); n_W <- ncol(W); n_Z <- ncol(Z)
  #
  beta_mu <- theta[1:n_X]
  beta_gamma <- theta[(n_X+1):(n_X+n_W)]
  beta_lambda <- theta[(n_X+n_W+1):(n_X+n_W+n_Z)]
  #
  mu_i    <- ginv_mu(as.vector(X %*% beta_mu))
  gamma_i <- ginv_gamma(as.vector(W %*% beta_gamma))
  lambda_i<- ginv_lambda(as.vector(Z %*% beta_lambda))
  alpha_i <- -((mu_i^lambda_i)*log(tau))/((-log(mu_i))^(gamma_i))
  #
  eq<-function(x,parms)
  {((-log(x))^parms[3])*(x^(-parms[4]))+(1/parms[2])*log(parms[1])}
  y<-c()
  j<-1
  while(j<=n){
    u<-runif(1)
    parms<-c(u,alpha_i[j],gamma_i[j],lambda_i[j])
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

func_linkWZ<-function(W,Z)
{
  out<-c()
  if(is.null(W)){out$W<-matrix(1,n,1)}else{out$W=as.matrix(W)}
  if(is.null(Z)){out$Z<-matrix(1,n,1)}else{out$Z=as.matrix(Z)}
  out$n_W <- ncol(W)
  out$n_Z <- ncol(Z)
  #
  if(all(out$W[,1] == 1, na.rm = TRUE) & out$n_W==1){
    out$g_gamma <- out$ginv_gamma <- function(x) x
    out$link_gamma <- "identity"
  }else{
    out$g_gamma   <- function(l)   log(l)
    out$ginv_gamma <- function(eta) exp(eta)
    out$link_gamma <- "log"
  }
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
  p <- c(trimws(p), "1", "1")[1:3]
  get_info <- function(f){
    tt <- terms(as.formula(f))
    list(intercept = attr(tt, "intercept"),k = length(attr(tt, "term.labels")))
  }
  info <- list(mu     = get_info(p[1]),gamma  = get_info(paste("~", p[2])),
               lambda = get_info(paste("~", p[3])))
  gen <- function(info){
    X <- if (info$k > 0)
      matrix(runif(n * info$k), n) else NULL
    if (info$intercept)
      X <- cbind(rep(1,n), X)
    X
  }
  return(list(X = gen(info$mu),W = gen(info$gamma),Z = gen(info$lambda)))
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

samples_RQUMW <- function(theta, X, W, Z, n = 50, re = 100, tau = 0.5,
                          ginv_mu, ginv_gamma, ginv_lambda,set_seed = NULL,
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
                            as.numeric(rRQUMW(n = n, theta = theta, X = X, W = W, Z = Z,
                                              tau = tau, ginv_mu = ginv_mu,
                                              ginv_gamma = ginv_gamma,
                                              ginv_lambda = ginv_lambda))
  },T)
  foreach::registerDoSEQ()
  parallel::stopCluster(cl)
  amostra <- as.data.frame(do.call(cbind, samples_list))
  return(amostra)
}


## RQ-UMW monte carlo simulation function --------------------------------------

sim_est_RQUMW <- function(sample, theta, X, W, Z, n, tau,
                          re = 100, RF = 100, method = "BFGS",
                          g_mu, ginv_mu, g_gamma, ginv_gamma,
                          g_lambda, ginv_lambda, start.theta = NULL,
                          n_cores = (parallel::detectCores() - 1))
{
  require(foreach)
  opts <- progresso(iterations = re, sec = "[2/2]")
  len_k <- length(theta)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  time <- system.time({try(estimate_list <- foreach(k = 1:re,.packages = c("foreach"),
      .options.snow = opts,.export = c("EST_RQUMW", "llike_RQUMW", "is.positive",
                  "vscore_RQUMW", "hessian_RQUMW","test.fun", "which.NA")) %dopar% {
    EST <- rep(NA, len_k)
    if (k < RF) {
      EST <- suppressWarnings(EST_RQUMW(y = sample[, k], X = X, W = W, Z = Z,
          tau = tau,g_mu = g_mu, method = method,ginv_mu = ginv_mu,g_gamma = g_gamma,
          ginv_gamma = ginv_gamma,g_lambda = g_lambda, ginv_lambda = ginv_lambda,
          start.theta = start.theta, applic = FALSE))
    }
    EST <- unname(as.numeric(EST))
    data.frame(t(EST), n_NA = any(is.na(EST)))
  },T)})
  foreach::registerDoSEQ()
  parallel::stopCluster(cl)
  #
  estimate_df <- do.call(rbind, estimate_list)
  n_NA_total <- sum(estimate_df$n_NA)
  output1 <- estimate_df[, 1:len_k]
  output1[output1 == 0] <- NA
  output1 <- which.NA(output1)
  if (nrow(output1) > RF) output1 <- output1[1:RF, ]
  #
  tabela <- tab_est(output1, theta)
  row.names(tabela) <- c(param_names_RQUMW(X, W, Z))
  print(tabela)
  #
  sec <- round(time[["elapsed"]])
  time_str <- sprintf("%02dh %02dm %02ds",sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
  return(list(output = output1,table = tabela,no_converg = n_NA_total,time = time_str))
}

## parameter names -------------------------------------------------------------

param_names_RQUMW <- function(X, W, Z)
{
  has_intercept <- function(M) {any(X[,1]==1) & var(M[, 1]) == 0}
  idx <- function(M) {
    p <- ncol(M)
    if (p == 0) return(integer(0))
    if (has_intercept(M)) c("(icpt)", seq_len(p - 1)) else seq_len(p)
  }
  c(paste0("μ", idx(X)),paste0("γ", idx(W)),paste0("λ", idx(Z)))
}


# Simulate RQ-UMW ==============================================================

#' Monte Carlo Simulation Study for the RQUMW Model
#'
#' Conducts a Monte Carlo simulation study for the regression quantile
#' UMW (RQUMW) model under different sample sizes and quantile levels.
#'
#' @param f A model formula that specifies the covariate structure of the RQUMW
#' model.
#' @param theta Numeric vector containing the true parameter values used
#'   in the data-generating process.
#' @param n Integer vector indicating the sample sizes to be considered.
#' @param re Integer indicating the number of Monte Carlo replications.
#' @param RF Integer indicating the number of bootstrap or resampling
#'   replications used in the estimation procedure.
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
#' @param method Optimization method passed to \code{\link[stats]{optim}}.
#' @param start.theta Optional numeric vector of initial values for the
#'   regression parameters. If \code{NULL}, default values are internally
#'   computed.
#' @param set_seed Integer. Optional random seed for reproducibility.
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{sample}: generated samples from the RQUMW model for each
#'     combination of \code{n} and \code{tau};
#'   \item \code{sim}: simulation results obtained from the estimation
#'     procedure;
#'   \item \code{cen_cov}: generated covariates used in the simulations;
#'   \item \code{link_mu}: the link function specification used for the
#'     location parameter.
#' }
#'
#' @export
simulate_RQUMW <- function(f = y ~ X1 + X2|1|1,theta=c(0.5,-0.6,0.2,1.5,2.3),method="BFGS",
                           n=50,re=100,RF=100,tau=0.5,link_mu="logit",save=F,start.theta = NULL,
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
      X<-cen_cov[[a1]][[a2]]$X;W<-cen_cov[[a1]][[a2]]$W;Z<-cen_cov[[a1]][[a2]]$Z
      funcWZ<-func_linkWZ(W=W,Z=Z)
      #
      nsamplerq[[a1]][[a2]] <- samplerq <- suppressWarnings(
        samples_RQUMW(n=nc,theta=theta,X=X,W=W,Z=Z,ginv_mu=ginv_mu,
                      ginv_gamma=funcWZ$ginv_gamma,ginv_lambda=funcWZ$ginv_lambda,
                      tau=tauc,re=re,n_cores = n_cores,set_seed=set_seed))
      sim_rumw[[a1]][[a2]] <- suppressWarnings(
        sim_est_RQUMW(sample = samplerq, X = X, W = W, Z = Z, n = nc,method = method,
                           theta = theta, tau = tauc, re = re, RF = RF,start.theta = start.theta,
                           n_cores =  n_cores,g_mu = g_mu, g_gamma = funcWZ$g_gamma,
                           g_lambda = funcWZ$g_lambda, ginv_mu = ginv_mu,
                           ginv_gamma = funcWZ$ginv_gamma, ginv_lambda = funcWZ$ginv_lambda))
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


