# Other Functions ==============================================================

## Parallelism -----------------------------------------------------------------

progresso<-function(iterations,sec="")
{
  pb <- progress::progress_bar$new(
    format = paste0(sec,":percent [:bar] :elapsed | Eta: :eta"),
    total = iterations,    # 100
    width = 60)
  #foreach:
  progress <- function(){pb$tick()}
  opts <- list(progress = progress)
  return(opts)
}

## Function to save simulation -------------------------------------------------

save_cen<-function(outputUMW,RF,n,write="sim")
{
  l1 <- nrow(outputUMW$sim$output)
  if(l1>=RF){
    save(outputUMW,file=paste0(write,"_R",RF,"_n",n,".RData"))
  }
  if(l1<RF){
    save(outputUMW,file=paste0(write,"_R",l1,"_n",n,".RData"))
  }
}

## sample replicate function ----------------------------------------------------


sample_UMW <- function(theta = c(0.7, 1.3, 0.5),n = 40,re = 100,set_seed = NULL,
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
                          .export = c("rUMW")) %dopar% {
    as.numeric(rUMW(n = n, theta = theta))
  },T)
  foreach::registerDoSEQ()
  parallel::stopCluster(cl)
  amostra <- as.data.frame(do.call(cbind, samples_list))
  return(amostra)
}

## UMW monte carlo simulation function -----------------------------------------

sim_est_UMW <- function(sample, theta = c(0.7, 1.3, 0.5), n = 40, re = 100, RF = 100,
                        method = "BFGS", start.theta = c(1, 1),
                        n_cores = (parallel::detectCores() - 1))
{
  require(foreach)
  opts <- progresso(iterations = re, sec = "[2/2]")
  estimate<-bigstatsr::FBM(nrow=re,ncol=3)
  n_NA<-bigstatsr::FBM(nrow=re,ncol=1)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  time <- system.time({try(foreach(k = 1:re,
                             .packages = c("foreach"),.options.snow = opts,
                             .export = c("Est_UMW", "which.NA", "test.fun", "is.positive",
                                         "hessian_UMW", "vscore_UMW", "llike_UMW")) %dopar% {
     output<-as.data.frame(estimate[])
     output[output==0]<-NA
     l_out<-nrow(which.NA(output))
     if(l_out<RF){
       estimate[k,]<-EST<-suppressWarnings(
         Est_UMW(x = sample[, k], start.theta = start.theta, method = method, applic = FALSE)
       )
       if(any(is.na(EST))){n_NA[j]<-1}
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
  if(nrow(output1)>=RF){output1<-output1[1:RF,]}
  #
  tabela<-tab_est(output1,theta)
  row.names(tabela) <- c("α", "γ", "λ")
  print(tabela)
  #
  sec <- round(time[["elapsed"]])
  time_str <- sprintf("%02dh %02dm %02ds", sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
  return(list(output = output1, table = tabela, no_converg = n_NA_total, time = time_str))
}

## summary table of the estimation ---------------------------------------------

tab_est<-function(coeff,theta){
  vmean <- colMeans(coeff)
  vbias <- vmean - theta
  vsd <- apply(coeff, 2, sd)
  vmse <- vbias^2 + vsd^2
  vac <- moments::skewness(coeff)
  vk <- moments::kurtosis(coeff)
  tabela <- data.frame(theta = theta, Mean = vmean, Bias = vbias,
                       MSE = vmse, AS = vac, K = vk)
  tabela <- round(tabela, 4)
}

## all positions of the vector are positive ------------------------------------

is.positive<-function(a)
{
  k<-length(a)
  tmp<-sum(a>0)
  return(k==tmp)
}

## Convergence test ------------------------------------------------------------

test.fun<-function(object)
{
  if(class(object)=="list"){
    if(object$convergence==0){
      parameters<-try(object$par,T)
      hess<-try(object$hessian,T)
      var.coef<-try(diag(solve(-hess)),T)
      if(is.numeric(parameters)==TRUE){
        if(is.numeric(var.coef)==TRUE){
          if(is.positive(var.coef)==TRUE){
            z<-c(parameters)
            return(z)
          }else{return(FALSE)}
        }else{return(FALSE)}
      }else{return(FALSE)}
    }else{return(FALSE)}
  }else{return(FALSE)}
}


## remove lines with NA --------------------------------------------------------

which.NA<-function(x)
{
  x<-data.frame(x)
  lines<-GLDEX::which.na(x[,1])
  tmp<-length(lines)
  if(tmp>0){
    y<-x[-lines,]
  }
  else{y<-x}
  return(y)
}


## Simulate UMW ================================================================

#' Simulation Study for the Unit-Modified Weibull(UMW) Distribution
#'
#' Conducts a Monte Carlo simulation study for the UMW distribution
#' under different sample sizes.
#'
#' @param theta Numeric vector of length 3 containing the true parameter
#'   values \code{(alpha, gamma, lambda)}, all strictly positive.
#' @param n Integer vector indicating the sample sizes to be considered.
#' @param re Integer indicating the number of Monte Carlo replications.
#' @param RF Integer number indicating the desired number of replicates.
#' @param save Logical; if \code{TRUE}, the simulation results are saved
#'   to disk.
#' @param n_cores Integer indicating the number of CPU cores to be used
#'   for parallel computation.
#' @param method Optimization method used by \code{\link[stats]{optim}}.
#'   Possible values are \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"} and \code{"SANN"}.
#' @param cen_name Character string used as a prefix for saved files.
#' @param start.theta Numeric vector of initial values for
#'   \eqn{(\gamma, \lambda)}.
#' @param set_seed Integer. Optional random seed for reproducibility.
#'
#' @return
#' A named list where each element corresponds to a sample size in
#' \code{n}. For each sample size, the list contains:
#' \itemize{
#'   \item \code{sample}: a numeric matrix containing the generated samples
#'   from the UMW distribution. Samples are stored column-wise, so that each
#'   column corresponds to one random sample.
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
#'   }
#' }
#'
#' @examples
#' theta = c(0.7,1.3,0.5)
#'
#' # Example without saving
#' simUMW<-simulate_UMW(theta,n = c(40,80),re = 100,RF = 100,method="BFGS",set_seed = 123)
#'
#' # Example saving
#' simUMW<-simulate_UMW(theta,n = c(40,80),re = 100,RF = 100,method="BFGS",set_seed = 123,
#' save = TRUE,cen_name = "sim")
#'
#' @export
simulate_UMW <- function(theta=c(0.7,1.3,0.5),n=40,re=100,RF=100,save=F,start.theta=c(1,1),
                         method="BFGS",n_cores = (parallel::detectCores()-1),
                         cen_name="sim",set_seed=NULL){
  if (!is.null(set_seed)){
    RNGkind("L'Ecuyer-CMRG")
    set.seed(set_seed)
  }
  outputUMW <- list()
  for (nc in n){
    cat("
n =",nc,"----
   ")
    outputUMW[[paste0("n", nc)]]$sample <- sample1 <- sample_UMW(theta = theta, n = nc,re = re,
                                                                 n_cores = n_cores,set_seed=set_seed)
    outputUMW[[paste0("n", nc)]]$sim <- suppressWarnings(sim_est_UMW(sample = sample1,theta = theta,
                                                 n = nc,re = re, RF = RF, n_cores = n_cores,start.theta=start.theta))
    # Save
    if(save==T){
      save_cen(outputUMW = outputUMW[[paste0("n", nc)]], RF = RF,  n = nc, write = paste0(cen_name,"_UMW"))
    }
  }
  invisible(outputUMW)
}


