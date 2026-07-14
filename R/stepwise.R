#' @export
stepwise_RQUMW <- function(data, response_var,
                           base_X = "1", base_Z = "1",
                           candidate_X=NULL, candidate_Z=NULL,
                           tau = 0.5, p_max = 0.1,link_mu = "probit",
                           intercept_X = TRUE, intercept_Z = TRUE, ...) {

  cat("=============================================================\n")
  cat(sprintf("Stepwise RQ-UMW (AIC and p-value <= %.2f)\n", p_max))
  cat("=============================================================\n\n")

  current_X <- base_X
  current_Z <- base_Z
  if(is.null(candidate_X)){candidate_X = setdiff(names(data), response_var)}
  if(is.null(candidate_Z)){candidate_Z = setdiff(names(data), response_var)}
  # FUNÇÃO DE FÓRMULA: Agora ela sabe se o bloco ainda está puramente na base inicial
  make_formula_str <- function(cx, cz, is_initial = FALSE) {

    # --- Bloco X ---
    # Se for o modelo inicial OU se o bloco X ainda não recebeu nenhuma variável candidata (continua sendo a base "1")
    if ((is_initial || all(cx == "1")) && ("1" %in% base_X)) {
      part_X <- "1"
    } else {
      if (!intercept_X) cx <- cx[cx != "1"]

      if (length(cx) == 0 || is.null(cx)) {
        part_X <- if (intercept_X) "1" else "0"
      } else {
        part_X <- paste(cx, collapse = " + ")
        if (!intercept_X) part_X <- paste(part_X, "- 1")
      }
    }

    # --- Bloco Z ---
    if ((is_initial || all(cz == "1")) && ("1" %in% base_Z)) {
      part_Z <- "1"
    } else {
      if (!intercept_Z) cz <- cz[cz != "1"]

      if (length(cz) == 0 || is.null(cz)) {
        part_Z <- if (intercept_Z) "1" else "0"
      } else {
        part_Z <- paste(cz, collapse = " + ")
        if (!intercept_Z) part_Z <- paste(part_Z, "- 1")
      }
    }

    paste0(response_var, " ~ ", part_X, " | ", part_Z)
  }

  clean_p_val <- function(p_str) {
    if (is.null(p_str) || is.na(p_str)) return(NULL)
    p_clean <- gsub("<", "", p_str)
    p_clean <- trimws(p_clean)
    return(as.numeric(p_clean))
  }

  get_p_value <- function(model_fit, var_name, prefix) {
    target_name <- paste0(prefix, ":", var_name)
    coef_table <- model_fit$coef

    if (target_name %in% rownames(coef_table)) {
      col_idx <- which(colnames(coef_table) == "Pr(>|z|)")
      if (length(col_idx) == 0) col_idx <- ncol(coef_table)
      p_raw <- coef_table[target_name, col_idx]
      return(clean_p_val(p_raw))
    }
    return(NULL)
  }

  # 1. AVALIA O MODELO BASE
  f_init <- make_formula_str(current_X, current_Z, is_initial = TRUE)
  fit_init <- tryCatch({
    fit_RQUMW(f = f_init, data = data, tau = tau,link_mu = link_mu, printmodel = FALSE, ...)
  }, error = function(e) NULL)

  if (is.null(fit_init)) stop("It was not possible to achieve model convergence with the provided base structure.")

  best_aic <- fit_init$metrics$aic
  cat(sprintf("Initial Base Model: (AIC: %.3f) \n %s \n\n", best_aic, f_init))

  converged_global <- FALSE
  ciclo <- 1

  while (!converged_global) {
    cat(sprintf("--------------- CYCLE %d ---------------\n", ciclo))
    aic_start_cycle <- best_aic

    # ------------------------------------------------------------------------
    # FASE A: Otimizar o bloco X
    # ------------------------------------------------------------------------
    cat("-> Optimizing block X...\n")
    remaining_X <- setdiff(candidate_X, current_X)
    best_new_X <- NULL
    chosen_X_aic <- best_aic

    if (length(remaining_X) > 0) {
      for (var in remaining_X) {
        test_X <- if (all(current_X == "1")) var else c(current_X, var)
        f_test <- make_formula_str(test_X, current_Z, is_initial = FALSE)

        fit_test <- tryCatch({
          fit_RQUMW(f = f_test, data = data, tau = tau,link_mu = link_mu, printmodel = FALSE, ...)
        }, error = function(e) NULL)

        if (!is.null(fit_test)) {
          test_aic <- fit_test$metrics$aic
          p_val <- get_p_value(fit_test, var, "β")

          if (test_aic < chosen_X_aic && !is.null(p_val) && p_val <= p_max) {
            chosen_X_aic <- test_aic
            best_new_X <- var
            saved_p_X <- p_val
          }
        }
      }

      if (!is.null(best_new_X)) {
        current_X <- if (all(current_X == "1")) best_new_X else c(current_X, best_new_X)
        best_aic <- chosen_X_aic
        cat(sprintf("   (+) Included variable '%s' (p-value: %.3f).\n New AIC: %.4f\n", best_new_X, saved_p_X, best_aic))
      } else {
        cat("   ( ) No additional variable from block X met the criteria.\n")
      }
    }

    # ------------------------------------------------------------------------
    # FASE B: Otimizar o bloco Z
    # ------------------------------------------------------------------------
    cat("-> Optimizing block Z...\n")
    remaining_Z <- setdiff(candidate_Z, current_Z)
    best_new_Z <- NULL
    chosen_Z_aic <- best_aic

    if (length(remaining_Z) > 0) {
      for (var in remaining_Z) {
        test_Z <- if (all(current_Z == "1")) var else c(current_Z, var)
        f_test <- make_formula_str(current_X, test_Z, is_initial = FALSE)

        fit_test <- tryCatch({
          fit_RQUMW(f = f_test, data = data, tau = tau, link_mu = link_mu, printmodel = FALSE, ...)
        }, error = function(e) NULL)

        if (!is.null(fit_test)) {
          test_aic <- fit_test$metrics$aic
          p_val <- get_p_value(fit_test, var, "δ")

          if (test_aic < chosen_Z_aic && !is.null(p_val) && p_val <= p_max) {
            chosen_Z_aic <- test_aic
            best_new_Z <- var
            saved_p_Z <- p_val
          }
        }
      }

      if (!is.null(best_new_Z)) {
        current_Z <- if (all(current_Z == "1")) best_new_Z else c(current_Z, best_new_Z)
        best_aic <- chosen_Z_aic
        cat(sprintf("   (+) Included variable '%s' (p-value: %.3f).\n New AIC: %.4f\n", best_new_Z, saved_p_Z, best_aic))
      } else {
        cat("   ( ) No additional variable from block Z met the criteria.\n")
      }
    }

    if (best_aic >= aic_start_cycle) {
      converged_global <- TRUE
      cat("\n=> Stopping criterion met. No variables added.\n")
    } else {
      ciclo <- ciclo + 1
    }
  }

  final_formula <- make_formula_str(current_X, current_Z, is_initial = FALSE)
  cat("\n=============================================================\n")
  cat("Final model selected:\n")
  cat(final_formula, "\n")

  final_fit <- fit_RQUMW(f = final_formula, data = data, tau = tau, link_mu = link_mu, printmodel = TRUE, ...)
  invisible(final_fit)
}

# fit<stepwise_RQUMW(data = d_test,response_var = "acuracy",tau = 0.5,link_mu = "logit",
#                   base_X = "1",
#                   base_Z = "1",
#                   candidate_X = c("dislex", "QI", "QI2", "intera", "intera2"),
#                   candidate_Z = c("dislex", "QI", "QI2", "intera", "intera2"),
#                   intercept_X = TRUE,intercept_Z = FALSE,p_max = 0.1)
# a$coef
# c("dislex", "QI2", "intera2")
# c("dislex", "QI", "QI2", "intera", "intera2")



