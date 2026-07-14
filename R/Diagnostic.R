
# Informação matrix equality (IME) =============================================

#' Information Matrix Equality diagnostics for RQ-UMW models
#'
#' Computes case-deletion diagnostics based on the Information Matrix Equality
#' (IME) for RQ-UMW models, including the statistics
#' \eqn{D_t}, \eqn{S_{3,t}} and \eqn{S_{5,t}}.
#'
#' @param fit An object of class \code{RQUMW} returned by \code{fit_RQUMW()}.
#' @param M Integer vector indicating which diagnostics to compute:
#'   \code{1} for \eqn{D_t}, \code{2} for \eqn{S_{3,t}} and \code{3} for \eqn{S_{5,t}}.
#' @param base_size Base font size for diagnostic plots.
#' @param graphic Logical; if \code{TRUE}, diagnostic plots are produced.
#' @param print Logical; if \code{TRUE}, results are printed to the console.
#' @param act_parallel Logical; if \code{TRUE}, parallel computation is used.
#' @param g_num Logical; if \code{TRUE}, observation indices are displayed in plots.
#' @param n_cores Number of cores to be used in parallel computation.
#'
#' @return A list containing:
#' \itemize{
#'   \item diagnostic measures (\eqn{D_t}, \eqn{S_{3,t}}, \eqn{S_{5,t}});
#'   \item identified influential observations;
#'   \item information matrix components;
#'   \item diagnostic plots (if requested).
#' }
#'
#' @examples
#' library(UMW)
#'
#' set.seed(25)
#' X <- runif(100, 0, 1)
#' y <- 0.3 + 0.7*X + rnorm(100, 0, 0.05)
#' y <- pmin(pmax(y, 0.01), 0.99)
#' fit<-fit_RQUMW(f = y ~ X)
#'
#' diagn_IME_RQUMW(fit,M = c(1,2,3))
#'
#' diagn_IME_RQUMW(fit,M = c(2,3),g_num = FALSE)
#'
#' @export

diagn_IME_RQUMW<-function(fit,M=c(1,2,3),base_size=12,graphic = TRUE,print = TRUE,
                          act_parallel = TRUE,g_num = TRUE, single_obs = NULL,
                          influence = TRUE, n_cores = (parallel::detectCores()-1)){
  if (!inherits(fit, "RQ-UMW")) {
    stop("Object must be of class 'RQ-UMW'.")
  }
  if(!any(M %in% c(1,2,3))){
    stop("Return some of the methods: 1 = D_t, 2 = M3_t and/or 3 = M5_t")
  }
  if (!is.null(single_obs) && length(single_obs) != 1) {
    stop("Error: 'single_obs' must be NULL or contain exactly one observation.")
  }
  out<-c()
  n<-fit$n
  # m3 calculation ----------------------------------------------------------- #
  calc_m3 <- calc_m3_RQUMW(fit$pars,fit$y,fit$X,fit$Z,tau=fit$quantile,
                          g_mu=fit$link$g_mu,g_lambda=fit$link$g_lambda,
                          ginv_mu=fit$link$ginv_mu,ginv_lambda=fit$link$ginv_lambda)
  out$m3 <- calc_m3$m3
  out$A_n <- calc_m3$A_n
  out$B_n <- calc_m3$B_n
  out$C3_n <- calc_m3$C3_n
  # -------------------------------------------------------------------------- #
  if(!is.null(single_obs)){
    k <- n <- single_obs
    act_parallel <- FALSE
    influence <- FALSE
  }else{k<-1}
  require(foreach)
  opts <- progresso(iterations = n)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  `%op%` <- if (act_parallel) `%dopar%` else `%do%`
  time <- system.time({try(metrics_list <- foreach(i = k:n,.packages = c("foreach"),
      .options.snow = opts,.export = c("EST_RQUMW", "llike_RQUMW",
      "vscore_RQUMW","hessian_RQUMW", "test.fun", "which.NA","calc_m3_RQUMW", "vech")) %op% {
      D  <- m3 <- m5 <- NA
      y_i <- fit$y[-i]
      X_i <- fit$X[-i, , drop = FALSE]
      Z_i <- fit$Z[-i, , drop = FALSE]
      mod_i <- EST_RQUMW(y = y_i, X = X_i, Z = Z_i,tau = fit$quantile,
          applic = FALSE,start.theta = fit$pars, g_mu = fit$link$g_mu,
          g_lambda = fit$link$g_lambda,ginv_mu = fit$link$ginv_mu,
          ginv_lambda = fit$link$ginv_lambda,method = fit$method)
      res_i <- calc_m3_RQUMW(mod_i, y_i, X_i, Z_i,tau = fit$quantile,
          g_mu = fit$link$g_mu,g_lambda = fit$link$g_lambda,
          ginv_mu = fit$link$ginv_mu,ginv_lambda = fit$link$ginv_lambda)
      diff_theta <- mod_i - fit$pars
      if (any(M %in% c(1, 3))) {
        D <- (fit$n - 1) * as.numeric(t(diff_theta) %*% (-res_i$A_n) %*% diff_theta)
      }
      if (any(M %in% 2)) {m3 <- res_i$m3}
      if (any(M %in% 3)) {
        m5 <- ((fit$n - 1) / 2) * as.numeric(t(diff_theta) %*% (-res_i$A_n +
                                              res_i$B_n) %*% diff_theta)
      }
      data.frame(D = D, m3 = m3, m5 = m5)
    },TRUE)
  })
  foreach::registerDoSEQ()
  parallel::stopCluster(cl)
  metrics_m <- do.call(rbind, metrics_list)
  # D_i ---------------------------------------------------------------------- #
  if(any(M %in% c(1,3))){out$D_t<-as.vector(metrics_m[,1])}
    if(influence == TRUE && any(M %in% 1)){
      D_i_infl<-g_threshold_pos(M_i = out$D_t,leg_y = expression(D[t]),
                        graphic = graphic,res_num = g_num,base_size = base_size)
      out$obs_D_t<-D_i_infl$Obs_infl
      out$threshold_D_t <- D_i_infl$threshold
      if(graphic == TRUE){
        out$graphics$D_t <- D_i_infl$graphcs_M
      }
    }
  # I1 & I2: s3_i ------------------------------------------------------------ #
  if(any(M %in% 2)){
    out$m3_t <- as.vector(metrics_m[,2])
    s3 <- out$m3_t / out$m3
    out$s3_t <- ifelse(s3 < 1, 1/s3, s3) - 1
    if(influence == TRUE){
      s3_t_infl <- g_threshold_pos(M_i = out$s3_t,leg_y = bquote(S[3*","*t]),
                        graphic = graphic,res_num = g_num,base_size = base_size,
                        weight = function(n){log(2*n)})
      out$obs_s3_t <- s3_t_infl$Obs_infl
      out$threshold_s3_t <- s3_t_infl$threshold
      if(graphic == TRUE){
        out$graphics$s3_t<-s3_t_infl$graphcs_M
      }
    }
    # if(influence == TRUE){
    #   out$I_3<-interval_cook(sm = out$s3_t,Type = 1)
    #   if(graphic == TRUE){
    #     out$graphics$gI_3<-g_diagnostic_RQUMW(sm_i = out$s3_t,v = 1,
    #                        I1 = out$I_3$I1,leg_y = bquote(S[3*","*t]),n = fit$n,
    #                        y = fit$y,res_num = g_num,base_size = base_size)
    #   }
    # }
  }
  # I1 & I2: s5_i ------------------------------------------------------------ #
  if(any(M %in% 3)){
    out$Dm_t<-as.vector(metrics_m[,3])
    out$s5_t <- abs(out$Dm_t - out$D_t)
    if(influence == TRUE){
      s5_t_infl<-g_threshold_pos(M_i = out$s5_t,leg_y = bquote(S[5*","*t]),
                                graphic = graphic,res_num = g_num,base_size = base_size,
                                weight = function(n){log(n)^1.1})
      out$obs_s5_t<-s5_t_infl$Obs_infl
      out$threshold_s5_t <- s5_t_infl$threshold
      if(graphic == TRUE){
        out$graphics$s5_t<-s5_t_infl$graphcs_M
      }
    }
    # if(influence == TRUE){
    #   out$I_5<-interval_cook(sm = out$s5_t,Type = 2)
    #   if(graphic == TRUE){
    #     out$graphics$gI_5<- g_diagnostic_RQUMW(sm_i = out$s5_t,v = 0,
    #                         I1 = out$I_5$I1,leg_y = bquote(S[5*","*t]),n = fit$n,
    #                         y = fit$y,res_num = g_num,base_size = base_size)
    #   }
    # }
  }
  # Return ------------------------------------------------------------------- #
  if(influence == TRUE && print == TRUE){
    cat("Observation(s) considered as possible outlier: \n")
    if(any(M %in% 1)){cat("D_t:",out$obs_D_t,"\n")}
    if(any(M %in% 2)){cat("Measure 3:",out$obs_s3_t,"\n")}
    if(any(M %in% 3)){cat("Measure 5:",out$obs_s5_t,"\n")}
    cat("\n")
    if(graphic == TRUE){
      print(suppressWarnings(cowplot::plot_grid(plotlist = out$graphics,
                                                ncol = length(out$graphics))))
    }
  }
  invisible(out)
}


## Measure function 3 (m3): ----------------------------------------------------

calc_m3_RQUMW <- function(theta_hat, y, X, Z, tau, g_mu, g_gamma, g_lambda,
                         ginv_mu, ginv_gamma, ginv_lambda){
  n <- length(y)
  k <- length(theta_hat)
  S <- vscore_RQUMW(theta_hat, y, X, Z, tau=tau, g_mu = g_mu, ginv_mu = ginv_mu,
                    g_lambda = g_lambda,ginv_lambda = ginv_lambda, vsmatrix = T)
  B_n <- (t(S) %*% (S)) / n
  A_n <- hessian_RQUMW(theta_hat, y, X, Z, g_mu = g_mu, ginv_mu = ginv_mu,
                    g_lambda = g_lambda,ginv_lambda = ginv_lambda,tau=tau) / n
  P <- chol(-A_n)
  P_inv <- solve(P)
  C3_n <- P_inv %*% B_n %*% t(P_inv)
  m3 <- sqrt(sum(vech(C3_n - diag(k))^2))
  return(list(m3=m3, C3_n=C3_n, B_n=B_n, A_n=A_n))
}

## vech function: --------------------------------------------------------------

vech <- function(M){
  idx <- lower.tri(M, diag=TRUE)
  return(M[idx])
}

## diagnostic graphics function ------------------------------------------------

g_diagnostic_RQUMW<-function(sm_i,v,I1,n,y,res_num = T,leg_y = "S_ji",base_size=12,
                             color_l = "gray35",linetype_l = "dashed"){
  df<-data.frame(sm_i=sm_i,nk=1:n,y=y) |>
    dplyr::mutate(n_infl=dplyr::case_when(sm_i > I1[[1]] & sm_i < I1[[2]] ~ NA,TRUE ~ nk))
  ggplot2::ggplot(df, ggplot2::aes(x=nk,y=sm_i))+
    ggplot2::geom_linerange(ggplot2::aes(ymin = v, ymax = sm_i,color = is.na(n_infl)),
                            linewidth = 1,show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("FALSE" = "red3", "TRUE" = "black")) +
    ggplot2::geom_hline(yintercept = I1[[2]],color = color_l,linewidth = 0.3,
                        linetype = linetype_l) +
    ggplot2::geom_line(ggplot2::aes(x = nk, y = v),colour="gray40",linewidth=0.15,
                       linetype="solid") +
    ggplot2::geom_hline(yintercept = I1[[1]],color = color_l,linewidth = 0.3,
                        linetype = linetype_l) +
    ggplot2::labs(x = "Index", y = leg_y)+
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    if(res_num == T){ggplot2::geom_text(ggplot2::aes(label = n_infl),
                                        vjust = ifelse(sm_i < v, 1.35, -0.3),size=2.5)}
}

## Graf cook -------------------------------------------------------------------

g_threshold_pos<-function(M_i,leg_y = "y",res_num=T,base_size=12,graphic=T,
                          weight=function(n){log(n)},threshold=NULL,
                          color_l="gray35",linetype_l = "dashed"){
  out <- c()
  n <- length(M_i)
  if (is.null(threshold)) {threshold <- mean(M_i) + (weight(n) * sd(M_i))}
  df_M_i <- data.frame(Obs = 1:n,M_i = M_i) |>
    dplyr::mutate(n_infl=dplyr::case_when(M_i <= threshold ~ NA,TRUE ~ Obs))
  obs_infl <- as.numeric(df_M_i$n_infl[!is.na(df_M_i$n_infl)])
  if(length(obs_infl) == 0){(obs_infl<-NA_integer_)}else{obs_infl}
  out$Obs_infl <- obs_infl
  out$threshold <- threshold
  if(graphic == T){
    out$graphcs_M<-ggplot2::ggplot(df_M_i, ggplot2::aes(x = Obs, y = M_i)) +
      ggplot2::geom_line(ggplot2::aes(x = Obs, y = 0),colour="gray40",linewidth=0.15,
                         linetype="solid") +
      ggplot2::geom_linerange(ggplot2::aes(ymin = 0, ymax = M_i,color = is.na(n_infl)),
                              linewidth = 1,show.legend = FALSE) +
      ggplot2::scale_color_manual(values = c("FALSE" = "red3", "TRUE" = "black")) +
      ggplot2::labs(x = "Index", y = leg_y) +
      ggplot2::geom_hline(yintercept = threshold,color = color_l,linewidth = 0.3,
        linetype = linetype_l) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())+
      if(res_num == T){ggplot2::geom_text(ggplot2::aes(label = n_infl),
                                          vjust=-0.3,size=2.5)}
  }
  return(out)
}


## I1 & I2 ---------------------------------------------------------------------

interval_cook<-function(sm,Type,m_orig=NULL){
  # S3_i
  n <- length(sm)
  if(Type == 1){
    weight<-sqrt(n)/2
    I1 <- list(1 - weight*diff(quantile(sm, c(0.125,0.5))),
            1 + weight*diff(quantile(sm, c(0.5,0.875))))
  }
  # S5_i
  if(Type == 2){
    weight<-(n)/3
    I1 <- list(0 - weight*diff(quantile(sm, c(0.25,0.5))),
            0 + weight*diff(quantile(sm, c(0.5,0.75))))
  }
  # I_Max
  if(Type == 3){
    weight<-(log(n))^2/5
    I1 <- list((0 - weight*diff(quantile(sm,c(0.125,0.5)))),
            (0 + weight*diff(quantile(sm,c(0.5,0.875)))))
  }
  # pseudo_JaB
  if(Type == 4){
    n <- length(m_orig)
    q_sig<-1-(n/(n+1))^2
    I1<-quantile(sm,c(q_sig,1-(q_sig)))
    sm<-m_orig
  }
  atipicas_I1 <- which(sm < I1[[1]] | sm > I1[[2]])
  if (length(atipicas_I1) == 0){(atipicas_I1<-NA_integer_)}else{atipicas_I1}
  return(list(Obs_I1=atipicas_I1,I1=I1))
}


# Local influence RQ-UMW =======================================================

#' Local influence diagnostics for RQ-UMW models
#'
#' Computes local influence diagnostics for RQ-UMW models using case-weight
#' perturbation schemes, including the curvature measure \eqn{C_t}
#' and the maximum normal curvature direction \eqn{I_{max}}.
#'
#' @param fit An object of class \code{RQUMW} returned by \code{fit_RQUMW()}.
#' @param M Integer vector indicating which diagnostics to compute:
#'   \code{1} for \eqn{C_t} and \code{2} for \eqn{I_{max}}.
#' @param k Index of the parameter to be perturbed. Use \code{0} for all parameters.
#' @param g_num Logical; if \code{TRUE}, observation indices are displayed in plots.
#' @param base_size Base font size for diagnostic plots.
#' @param print Logical; if \code{TRUE}, results are printed to the console.
#' @param graphic Logical; if \code{TRUE}, diagnostic plots are produced.
#'
#' @return A list containing:
#' \itemize{
#'   \item diagnostic measures (\eqn{C_t}, \eqn{I_{max}})
#'   \item local influence measures;
#'   \item identified influential observations;
#'   \item diagnostic plots (if requested).
#' }
#'
#' @examples
#' library(UMW)
#'
#' set.seed(25)
#' X <- runif(100, 0, 1)
#' y <- 0.3 + 0.7*X + rnorm(100, 0, 0.05)
#' y <- pmin(pmax(y, 0.01), 0.99)
#' fit<-fit_RQUMW(f = y ~ X)
#'
#' UMW::diagn_LocI_RQUMW(fit,M = c(1,2),k = 0)
#'
#' UMW::diagn_LocI_RQUMW(fit,M = c(1),k = 2,g_num = F)
#'
#'
#' @export
diagn_LocI_RQUMW <- function(fit,M = c(1,2),k = 0,g_num = TRUE,base_size = 12,
                             print = TRUE,graphic = TRUE,influence = TRUE, single_obs = NULL) {
  if (!inherits(fit, "RQ-UMW")) {
    stop("Object must be of class 'RQ-UMW'.")
  }
  if(!any(M %in% c(1,2))){
    stop("Return some of the methods: 1 = C_t and/or 2 = I_max")
  }
  n_k<-length(fit$pars)
  if (!is.numeric(k) || k < 0 || k > n_k || k %% 1 != 0) {
    stop(paste0("'k' must be a numeric integer between 0 and ",n_k," (number of estimated parameters)."))
  }
  out<-c()
  n<-fit$n
  # case–weight -------------------------------------------------------------- #
  out$Delta <- Delta <- t(vscore_RQUMW(fit$pars,y = fit$y,X = fit$X,
                          Z = fit$Z,tau=fit$quantile, g_mu = fit$link$g_mu,
                          g_lambda = fit$link$g_lambda,ginv_mu = fit$link$ginv_mu,
                          ginv_lambda = fit$link$ginv_lambda, vsmatrix = T))
  # -------------------------------------------------------------------------- #
  I_inv <- solve(-fit$hessian)
  I_inv_neg<-matrix(0,nrow = n_k,ncol = n_k)
  if(k!=0){I_inv_neg[k,k]<-I_inv[k,k]}
  names_pars <- c("bold(theta)","gamma","lambda",
                  sapply(0:(n_k-3), function(i) paste0("beta[", i, "]")))
  param_symbols <- c("\u03B8","\u03B3","\u03BB",
                     sapply(0:(n_k-3), function(i) paste0("\u0392", i))
  )
  if(!is.null(single_obs)){
    h <- n <- single_obs
    act_parallel <- FALSE
    influence <- FALSE
  }else{h<-1}
  # C_t ---------------------------------------------------------------------- #
  if(any(M %in% c(1))){
    C_t<- c()
    for(t in h:n){
      C_t[t] <- 2 * abs(as.numeric(t(Delta[,t]) %*% (I_inv-I_inv_neg) %*% Delta[,t]))
    }
    out$C_t <- C_t
    if(influence == TRUE){
      C_t_infl<-g_threshold_pos(M_i = out$C_t,leg_y = expression(C[t]),graphic = graphic,
                                res_num = g_num,base_size = base_size)
      out$outliers$C_t <- C_t_infl$Obs_infl
      out$threshold$C_t <- C_t_infl$threshold
      if(graphic == TRUE){
        out$graphics$C_t<-C_t_infl$graphcs_M
      }
    }
  }
  # I_max -------------------------------------------------------------------- #
  if(any(M %in% c(2))){
    eigen_res <- eigen(-t(Delta) %*% (I_inv-I_inv_neg) %*% Delta)
    out$C_max <- max(as.numeric(eigen_res$values))
    out$I_max <- as.numeric(eigen_res$vectors[, which.max(eigen_res$values)])[h:n]
    if(influence == TRUE){
      out$I_max_lb<-interval_cook(sm = out$I_max,Type = 3)
      out$outliers$I_max <- out$I_max_lb$Obs_I1
      if(graphic == TRUE){
        out$graphics$I_max<-g_diagnostic_RQUMW(sm_i = out$I_max,v = 0,I1 = out$I_max_lb$I1,
                             leg_y = expression(I[max]),n = fit$n,y = fit$y,
                             res_num = g_num,base_size = base_size)
      }
    }
  }
  # Return ------------------------------------------------------------------- #
  if(influence == TRUE && print == TRUE){
    cat("Observation(s) considered as possible outlier for",param_symbols[k+1],":\n")
    if(any(M %in% 1)){cat("C_t:",out$outliers$C_t,"\n")}
    if(any(M %in% 2)){cat("I_max:",out$outliers$I_max,"\n")}
    cat("\n")
    if(graphic == TRUE){
      print(suppressWarnings(cowplot::plot_grid(plotlist = out$graphics,
                                                ncol = length(out$graphics))))    }
  }
  invisible(out)
}

# Influence measured distance ==================================================

#' Distance-based influence diagnostics for RQ-UMW models
#'
#' Computes influence diagnostics based on distance measures between
#' fitted models with and without each observation. The implemented
#' distances are the Frèchet distance and the Rao distance.
#'
#' @param fit An object of class \code{RQUMW} returned by \code{fit_RQUMW()}.
#' @param M Integer vector indicating which distances to compute:
#'   \code{1} for Frèchet distance and \code{2} for Rao distance.
#' @param g_num Logical; if \code{TRUE}, observation indices are displayed in plots.
#' @param base_size Base font size for diagnostic plots.
#' @param act_parallel Logical; if \code{TRUE}, parallel computation is used.
#' @param print Logical; if \code{TRUE}, results are printed to the console.
#' @param graphic Logical; if \code{TRUE}, diagnostic plots are produced.
#' @param n_cores Number of cores to be used in parallel computation.
#'
#' @return A list containing:
#' \itemize{
#'   \item diagnostic measures (\eqn{Frèchet}, \eqn{Rao})
#'   \item distance-based influence measures;
#'   \item identified influential observations;
#'   \item diagnostic plots (if requested).
#' }
#'
#' @examples
#' library(UMW)
#'
#' set.seed(25)
#' X <- runif(100, 0, 1)
#' y <- 0.3 + 0.7*X + rnorm(100, 0, 0.05)
#' y <- pmin(pmax(y, 0.01), 0.99)
#' fit<-fit_RQUMW(f = y ~ X)
#'
#' UMW::diagn_DIST_RUMW(fit,M = c(1,2))
#'
#' UMW::diagn_DIST_RUMW(fit,M = c(2),g_num = F)
#'
#' @export
diagn_DIST_RQUMW <- function(fit, M = c(1,2), g_num = TRUE, base_size = 12,
                            act_parallel = TRUE,influence = TRUE,print=TRUE,
                            graphic = TRUE,n_cores = (parallel::detectCores()-1),
                            single_obs=NULL){
  if (!inherits(fit, "RQ-UMW")) {
    stop("Object must be of class 'RQ-UMW'.")
  }
  if(!any(M %in% c(1,2))){
    stop("Return some of the methods: 1 = Frèchet distance and/or 2 = Rao distance")
  }
  if (!is.null(single_obs) && length(single_obs) != 1) {
    stop("Error: 'single_obs' must be NULL or contain exactly one observation.")
  }
  out<-c()
  n<-fit$n
  S1<-solve(-fit$hessian)
  if(!is.null(single_obs)){
    k <- n <- single_obs
    act_parallel <- FALSE
    influence <- FALSE
  }else{k<-1}
  opts<-progresso(iterations = n)
  require(foreach)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  `%op%` <- if (act_parallel) `%dopar%` else `%do%`
  metrics_list <- foreach(i = k:n,.packages = c("foreach"),.options.snow = opts,
    .combine = rbind,.export = c("EST_RQUMW", "llike_RQUMW","vscore_RQUMW",
    "hessian_RQUMW","test.fun", "which.NA","frechet_simple","rao_simple")) %op% {
    D <- R <- NA_real_
    y_i <- fit$y[-i]
    X_i <- as.matrix(fit$X[-i, , drop = FALSE])
    Z_i <- as.matrix(fit$Z[-i, , drop = FALSE])
    mod_i <- EST_RQUMW(y = y_i, X = X_i, Z = Z_i,tau = fit$quantile,
                       applic = TRUE,start.theta = fit$pars, g_mu = fit$link$g_mu,
                       g_lambda = fit$link$g_lambda,ginv_mu = fit$link$ginv_mu,
                       ginv_lambda = fit$link$ginv_lambda,method = fit$method)
    S2 <- solve(-mod_i$hessian)
    if (any(M %in% 1)) {
      D <- frechet_simple(mu1 = fit$pars, mu2 = mod_i$pars,S1  = S1,S2  = S2)
    }
    if (any(M %in% 2)) {
      R <- rao_simple(S1 = S1, S2 = S2)
    }
    data.frame(D = D, R = R)
  }
  foreach::registerDoSEQ()
  parallel::stopCluster(cl)
  metrics_m <- as.data.frame(metrics_list)
  # Frèchet distance --------------------------------------------------------- #
  if(any(M %in% c(1))){
    out$frechet<-as.vector(metrics_m[,1])
    if(influence == TRUE){
      frechet_infl<-g_threshold_pos(M_i = out$frechet,leg_y = "Frèchet distance",
                        graphic = graphic,res_num = g_num,base_size = base_size,
                        weight = function(n){log(n)^1.1})
      out$outliers$frechet <- frechet_infl$Obs_infl
      out$threshold$frechet <- frechet_infl$threshold
      if(graphic == TRUE){
        out$graphics$frechet<-frechet_infl$graphcs_M
      }
    }
  }
  # Rao distance ------------------------------------------------------------- #
  if(any(M %in% c(2))){
    out$rao<-as.vector(metrics_m[,2])
    if(influence == TRUE){
      rao_infl<-g_threshold_pos(M_i = out$rao,leg_y = "Rao distance",
                        graphic = graphic,res_num = g_num,base_size = base_size)
      out$outliers$rao <- rao_infl$Obs_infl
      out$threshold$rao <- rao_infl$threshold
      if(graphic == TRUE){
        out$graphics$rao<-rao_infl$graphcs_M
      }
    }
  }
  # Return ------------------------------------------------------------------- #
  if(influence == TRUE && print == TRUE){
    cat("Observation(s) considered as possible outlier:\n")
    if(any(M %in% 1)){cat("Frechet:",out$outliers$frechet,"\n")}
    if(any(M %in% 2)){cat("Rao:",out$outliers$rao,"\n")}
    cat("\n")
    if(graphic == T){
      print(suppressWarnings(cowplot::plot_grid(plotlist = out$graphics,
                                                ncol = length(out$graphics))))}
  }
  out1<-out
}

## The Frèchet distance as influence measure -----------------------------------

frechet_simple <- function(mu1,mu2,S1,S2) {
  dmu2 <- sum((mu1-mu2)^2)
  trace_term <- sum(diag(S1 + S2 - 2 * expm::sqrtm(expm::sqrtm(S1) %*% S2 %*% expm::sqrtm(S1))))
  d2 <- sqrt(dmu2 + trace_term)
  return(as.numeric(d2))
}

## The Rao distance as influence measure ---------------------------------------

rao_simple <- function(S1,S2) {
  M <- S2 %*% solve(S1)
  lambda <- eigen(M, only.values = TRUE)$values
  rho <- sqrt((1/2) * sum(log(lambda)^2))
  return(rho)
}

# JAB ==========================================================================

#' Jackknife-after-Bootstrap influence diagnostics for RQ-UMW models
#'
#' Computes Jackknife-after-Bootstrap (JaB) influence diagnostics for several
#' influence measures derived from RQ-UMW models.
#'
#' @param fit An object of class \code{RQUMW} returned by \code{fit_RQUMW()}.
#' @param M Integer vector indicating which influence measures to evaluate:
#'   \code{1} C_t, \code{2} I_max, \code{3} D_t, \code{4} S3_t,
#'   \code{5} S5_t, \code{6} Frèchet, \code{7} Rao.
#' @param B Number of bootstrap replications.
#' @param graphic Logical; if \code{TRUE}, diagnostic plots are produced.
#' @param set_seed Optional integer to set a random seed for reproducibility.
#' @param save Logical; if \code{TRUE}, results are saved to disk.
#' @param name_JaB Character string used as prefix for saved objects.
#' @param n_cores Number of cores to be used in parallel computation.
#'
#' @return A list containing, for each selected measure:
#' \itemize{
#'   \item original influence measure;
#'   \item bootstrap distributions (\eqn{Fstar});
#'   \item JaB influence values (\eqn{Fstar_i});
#'   \item identified influential observations;
#'   \item computation time.
#' }
#'
#' @examples
#' library(UMW)
#'
#' set.seed(25)
#' X <- runif(100, 0, 1)
#' y <- 0.3 + 0.7*X + rnorm(100, 0, 0.05)
#' y <- pmin(pmax(y, 0.01), 0.99)
#' fit<-fit_RQUMW(f = y ~ X)
#'
#' diagn_JaB_RUMW(fit,M = 1,B = 200,set_seed = 25)
#'
#' diagn_JaB_RUMW(fit,M = c(1:7),B = 200,set_seed = 25,save = T,name_JaB = "JaB")
#'
#'
#' @export
diagn_JaB_RQUMW <- function(fit, M = c(1,2,3,4,5,6,7), B = 200, graphic = T,
                            set_seed=NULL,save = F,name_JaB = "JaB", prob = NULL,
                            n_cores = (parallel::detectCores()-1)) {
  if(!any(M %in% 1:7)){stop("Select measurements from 1 to 7.")}
  out<-list()
  for (Mc in M) {
    met<-c("C_t","I_max","D_t","S3_t","S5_t","Frechet","Rao")
    cat(met[Mc],"----\n")
    M_orig <- func_model(fit, Mc)
    if (!is.null(set_seed)) {set.seed(set_seed)}
    indices_boot <- replicate(B, sample(1:fit$n, replace = TRUE))
    opts<-progresso(iterations = B)
    require(foreach)
    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    b_time <- system.time({
      M_boot_list <- foreach(b = 1:B,.packages = c("foreach"),.options.snow = opts,
                     .export = c("func_model","diagn_LocI_RQUMW","diagn_IME_RQUMW",
                     "diagn_DIST_RQUMW")) %dopar% {
        fit2 <- fit
        fit2$data <- fit$data[indices_boot[, b], ]
        res <- try(as.numeric(func_model(fit = fit2, M = Mc, orig = FALSE)),TRUE)
        if (inherits(res, "try-error")) {rep(NA_real_, fit$n)}
        else{res}
      }
    })
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
    M_boot <- as.data.frame(M_boot_list)
    M_i_boot <- Fstar <- F_i <-c()
    for (i in 1:fit$n) {
      keep <- which(!apply(indices_boot, 2, function(col) i %in% col))
      M_i_boot[[i]] <- as.numeric(M_boot[i, keep])
      Fstar[[i]]<-ecdf(M_i_boot[[i]])
      F_i[i]<-Fstar[[i]](M_orig[i])
    }
    if(Mc==2){F_i <- ifelse(M_orig<0,1-F_i,F_i)}
    max_len <- max(lengths(M_i_boot))
    M_i_boot_pad <- M_i_boot |> lapply(\(x) { length(x) <- max_len; x }) |>
      as.data.frame() |> (\(df) { names(df) <- paste0("V", seq_along(df)); df })()
    sec <- (round(b_time[["elapsed"]]))
    time<-(sprintf("%02dh %02dm %02ds",sec %/% 3600,(sec %% 3600) %/% 60,sec %% 60))
    out[[met[Mc]]]<-list(M_orig = M_orig,M_t_boot = M_i_boot_pad,time = time, B = B,
                         Fstar = Fstar,F_t = F_i)
    #
    if(graphic==T){
      leg_y<-c(expression(JaB * " - " * C[t]),expression(JaB * " - " * I[max]),
               expression(JaB * " - " * D[t]),bquote(JaB ~ "-" ~ S[3*","*t]),
               bquote(JaB ~ "-" ~ S[5*","*t]), "JaB - Frèchet distance",
               "JaB - Rao distance")
      if(is.null(prob)){prob <- (fit$n/(fit$n+1))^2}
      max_fi <- max(F_i)
      y_smin <- floor(max_fi / 0.05) * 0.05 - 0.1
      y_smax <- if (sum(F_i > prob) > 5) max_fi + 0.05 else y_smin
      g_M<-g_threshold_pos(M_i = F_i,leg_y = leg_y[Mc],
                           threshold = prob,base_size = 14)
      out[[met[Mc]]]$threshold <- prob
      out[[met[Mc]]]$Obs_infl <- g_M$Obs_infl
      out[[met[Mc]]]$graphic <- g_M$graphcs_M +
        ggplot2::coord_cartesian(ylim = c(y_smax,1))+
        ggplot2::scale_y_continuous(labels = function(x) sprintf("%.2f", x))
      withCallingHandlers(print(out[[met[Mc]]]$graphic),
                          warning = function(w) invokeRestart("muffleWarning"))
      cat("Possible influential observations:",out[[met[Mc]]]$Obs_infl,"\n")
    }
    if(save==T){
      outputJaB<-out
      save(outputJaB,file=paste0(name_JaB,"_B",B,"_",met[Mc],".RData"))
    }
  }
  invisible(out)
}

# pseudo-JaB ===================================================================

#' @export
diagn_pseudoJaB_RQUMW <- function(fit, M = c(1,2,3,4,5,6,7), B = 200, graphic = T,
                       set_seed=NULL,save = F,name_JaB = "JaB",prob=NULL,
                       n_cores = (parallel::detectCores()-1)) {
  if(!any(M%in%1:7)){stop("Select measurements from 1 to 7.")}
  out<-list()
  for (Mc in M) {
    met<-c("C_t","I_max","D_t","S3_t","S5_t","Frechet","Rao")
    cat(met[Mc],"----\n")
    M_orig <- func_model(fit, Mc)
    if (!is.null(set_seed)) {set.seed(set_seed)}
    B_data <- replicate(B, rRQUMW(n = fit$n,theta = fit$pars,X = fit$X,
              Z = fit$Z,tau = fit$quantile,ginv_mu = fit$link$ginv_mu,
              ginv_lambda = fit$link$ginv_lambda,set_seed = NULL))
    ib <- sample(1:fit$n,B, replace = TRUE)
    opts<-progresso(iterations = B)
    require(foreach)
    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    b_time <- system.time({
      dist_list <- foreach(b = 1:B,.packages = c("foreach"),.options.snow = opts,
         .export = c("func_model","diagn_LocI_RQUMW","diagn_IME_RQUMW",
                      "diagn_DIST_RQUMW")) %dopar% {
                        fit2 <- fit
                        fit2$data[,1] <- B_data[,b]
                        res <- try(as.numeric(func_model(fit = fit2, M = Mc,
                                   orig = FALSE,single_obs = ib[b])),TRUE)
                        if (inherits(res, "try-error")) NA_real_ else res
      }
    })
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
    dist_vec <- unlist(dist_list)
    dist_vec <- dist_vec[!is.na(dist_vec)]
    if(length(dist_vec) == 0){stop("All pseudo-JaB replications failed.")}
    #
    sec <- (round(b_time[["elapsed"]]))
    time<-(sprintf("%02dh %02dm %02ds",sec %/% 3600,(sec %% 3600) %/% 60,sec %% 60))
    out[[met[Mc]]]<-list(M_orig = M_orig,pJaB_values = dist_vec,time = time, B = B)
    #
    if(graphic==T){
      leg_y<-c(expression(pJaB * " - " * C[t]),expression(pJaB * " - " * I[max]),
               expression(pJaB * " - " * D[t]),bquote(pJaB ~ "-" ~ S[3*","*t]),
               bquote(pJaB ~ "-" ~ S[5*","*t]), "pJaB - Frèchet distance",
               "pJaB - Rao distance")
      if(Mc %in% c(1,3:7)){
        if(is.null(prob)){prob <- (fit$n/(fit$n+1))^2}
        threshold1<-quantile(dist_vec,prob)
        g_M<-g_threshold_pos(M_i = M_orig,leg_y = leg_y[Mc],
                             threshold = threshold1,base_size = 14)
        out[[met[Mc]]]$threshold <- threshold1
        out[[met[Mc]]]$Obs_infl <- g_M$Obs_infl
        out[[met[Mc]]]$graphic <- g_M$graphcs_M
      }else{
        Ib <- interval_cook(sm = dist_vec,Type = 4,m_orig = M_orig)
        out[[met[Mc]]]$threshold <- threshold1 <- Ib$I1
        out[[met[Mc]]]$Obs_infl <- Ib$Obs_I1
        V <- ifelse(Mc == 4, 1, 0)
        out[[met[Mc]]]$graphic <- g_diagnostic_RQUMW(sm_i = M_orig,v = V,
                            I1 = round(threshold1,3),leg_y = leg_y[Mc],n = fit$n,
                            y = fit$y,res_num = TRUE,base_size = 14)
      }
      withCallingHandlers(print(out[[met[Mc]]]$graphic),
                          warning = function(w) invokeRestart("muffleWarning"))
      cat("Possible influential observations:",out[[met[Mc]]]$Obs_infl,"\n")
    }
    if(save==T){
      outputJaB<-out
      save(outputJaB,file=paste0(name_JaB,"_B",B,"_",met[Mc],".RData"))
    }
  }
  invisible(out)
}

# selection of the measure for JaB ---------------------------------------------
func_model <- function(fit, M = c(1,2,3,4,5,6,7),orig=T,single_obs=NULL) {
  if(orig==TRUE){mod<-fit}else{
    mod <- fit_RQUMW(fit$formula,data = fit$data,tau = fit$quantile,
                     link_mu = fit$link$link_mu,method = fit$method,printmodel = F)
  }
  if(M == 1){# C_t
    M_i <- diagn_LocI_RQUMW(mod,M = c(1),k = 0,single_obs=single_obs,influence=FALSE)$C_t
  }
  if(M == 2){# I_max
    M_i <- diagn_LocI_RQUMW(mod,M = c(2),k = 0,single_obs=single_obs,influence=FALSE)$I_max
  }
  if(M == 3){# D_i
    M_i <- diagn_IME_RQUMW(mod,M = c(1),single_obs=single_obs,influence=FALSE)$D_t
  }
  if(M == 4){# S3_i
    M_i <- diagn_IME_RQUMW(mod,M = c(2),single_obs=single_obs,influence=FALSE)$s3_t
  }
  if(M == 5){# S5_i
    M_i <- diagn_IME_RQUMW(mod,M = c(3),single_obs=single_obs,influence=FALSE)$s5_t
  }
  if(M == 6){# Frechet
    M_i <- diagn_DIST_RQUMW(mod,M = c(1),single_obs=single_obs,influence=FALSE)$frechet
  }
  if(M == 7){# Rao
    M_i <- diagn_DIST_RQUMW(mod,M = c(2),single_obs=single_obs,influence=FALSE)$rao
  }
  return(M_i)
}
# ---------------------------------------------------------------------------- #



