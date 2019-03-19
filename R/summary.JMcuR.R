#' Summary of a 'JMcuR' objects
#'
#' The function provides a summary of jointCureModel estimations.
#'
#' @param object an object inheriting from class 'JMcuR'
#' @param ... further arguments to be passed to or from other methods. They are ignored in this function.
#'
#' @return Returns NULL.
#'
#' @export
#'
#' @examples
summary.JMcuR <- function (object, ...)
  {
    if (!inherits(object, "JMcuR"))
      stop("Use only with 'JMcuR' objects.\n")

    #---- details of cure model
    if(object$jointCureModel=="MCmodel")
      cureModel <- "Mixture cure model (MCM)"
    if(object$jointCureModel=="JLCCmodel")
      cureModel <- "Joint latent class cure model (JLCCM)"
    if(object$jointCureModel=="FJCmodel")
      cureModel <- "Full joint cure model"
    if(object$param=="shared-RE")
      AssocStruc <- "Shared random effect(s)"
    if(object$param=="td-value")
      AssocStruc <- "Current value of longitudinal process"
    cat("#-- Statistical Model:", "\n")
        cat(paste("     -(Joint) cure model:", cureModel), "\n")
        if(object$jointCureModel=="FJCmodel")
          cat(paste("     -Association structure:", AssocStruc), "\n")
        if(object$jointCureModel %in%  c("FJCmodel","JLCCmodel")){
          cat(paste("     -Random effects assumed class-specific:", object$Sigma_d), "\n")
          cat(paste("     -Random effects assumed independent:", object$cov_prior=="inverse-gamma"), "\n")
        }
        cat(paste("     -Number of subjects:", object$control$n), "\n")
        cat(paste("     -Number of observations:", nrow(object$data)), "\n")
        cat("\n")

    #---- Parameter Estimations
    coefs <- object$coefficients
    strs <- object$StErr
    sds <- object$StDev
    CIs <- object$CIs
    # Incidence model
    incidence <- cbind("Value" = coefs$gamma,
                       "Std.Err" = strs$gamma, "Std.Dev" = sds$gamma,
                       "2.5%" = CIs$gamma[1, ], "97.5%" = CIs$gamma[2, ])
    cat("#-- Estimation of the incidence model:\n")
    prmatrix(incidence, na.print = "")
    cat("\n")
    # Latency model
    latency <- cbind("Value" = coefs$alpha,
                     "Std.Err" = strs$alpha, "Std.Dev" = sds$alpha,
                     "2.5%" = CIs$alpha[1, ], "97.5%" = CIs$alpha[2, ])
    latency <- rbind(latency,
                     shape = c(coefs$shape, strs$shape, sds$shape, CIs$shape[1 ], CIs$shape[2 ]))
    latency <- rbind(shape = latency["shape", ], latency[-nrow(latency), ])
    if(object$jointCureModel=="FJCmodel" && object$param=="td-value")
      row.names(latency)[grep("alphaL", row.names(latency), fixed = TRUE)] <- "Assoct."
    if(object$jointCureModel=="FJCmodel" && object$param=="shared-RE"){
      ii <- grep("alphaRE", row.names(latency), fixed = TRUE)
      row.names(latency)[ii] <- paste("Assoct.", colnames(object$x$U), sep = "")
    }
    cat("#-- Estimation of the latency model:\n")
    prmatrix(latency, na.print = "")
    cat("\n")
    cat("\n")

    #---- Longitudinal model
    if(object$jointCureModel!="MCmodel"){
      cat("#-- Estimation of the longitudinal model: \n")
      cat("\n")
      cat("Fixed effects for susceptible group (D=1): \n")
      coefs.beta.D1 <- cbind("Value" = coefs$beta.D1, "Std.Err" = strs$beta.D1,
                             "Std.Dev" = sds$beta.D1, "2.5%" = CIs$beta.D1[1, ],
                             "97.5%" = CIs$beta.D1[2, ])
      prmatrix(coefs.beta.D1, na.print = "")
      cat("\n")
      cat("Fixed effects for non susceptible group (D=0): \n")
      coefs.beta.D0 <- cbind("Value" = coefs$beta.D0, "Std.Err" = strs$beta.D0,
                             "Std.Dev" = sds$beta.D0, "2.5%" = CIs$beta.D0[1, ],
                             "97.5%" = CIs$beta.D0[2, ])
      prmatrix(coefs.beta.D0, na.print = "")
      cat("\n")
      if(ncol(object$x$U)>1){
        if(!object$Sigma_d){
          cat("Covariance matrix of the random-effects:\n")
          prmatrix(object$coefficients$Sigma2, na.print = "")
        }else{
          cat("Covariance matrix of the random-effects for susceptible group:\n")
          prmatrix(object$coefficients$Sigma2.D1, na.print = "")
          cat("\n")
          cat("Covariance matrix of the random-effects for non susceptible group:\n")
          prmatrix(object$coefficients$Sigma2.D0, na.print = "")
          cat("\n")
        }
        cat(paste("Residual standard error = ", object$coefficients$tau))
      }else{
        if(!object$Sigma_d){
          cat(paste("Variance of the random-effect for non susceptible group: ", object$coefficients$sigma2))
          cat("\n")
          cat(paste("Residual standard error = ", object$coefficients$tau))
        }else{
          cat(paste("Variance of the random-effect for susceptible group: ", object$coefficients$sigma2[1]))
          cat("\n")
          cat(paste("Variance of the random-effect for non susceptible group: ", object$coefficients$sigma2[2]))
          cat("\n")
          cat(paste("Residual standard error = ", object$coefficients$tau))
        }
      }

    }
  }
