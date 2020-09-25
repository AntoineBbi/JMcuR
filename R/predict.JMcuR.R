#' Predictions of the cure membership, the longitudinal or survival responses
#'
#' Predictions of the cure membership, the longitudinal or survival responses from the estimated joint longitudinal survival models with cure fraction
#'
#' @param object an object inheriting from class \code{JMcuR}
#' @param newdata a data frame in which to look for variables with which to predict (default is NULL, and the prediction is done on the data used to estimate the cure model)
#' @param type a character string indicating the type of predictions to compute, marginal or subject-specific.
#' @param state a (vector of) character string indicating the response ("cure","survival","longitudinal") to predict
#' @param level.interval a numeric scalar denoting the tolerance/confidence level for the credibility interval; the credible interval is return only for \cite{MCMC=TRUE}
#' @param Tsurv A numeric scalar denoting the time of interest associated with the prediction; 'Tsur' is NULL by default and then the censoring time is considered to predict the cure status.
#' @param Yt a numeric scalar denoting the value of the longitudinal measurement if type='profile'. yt is NULL by default
#' @param yVar a character string indicating the name of the variable in newdata that corresponds to the longitudinal outcome;
#' @param idVar a character string indicating the name of the variable in newdata that corresponds to the subject identifier; required when type = "Subject"
#' @param MCMC logical; if TRUE prediction procedure is done using Bayesian approach, and if FALSE (by default) the Frequentist approach is then considered for the prediction procedure
#' @param M a numerical scalar denoting the sampling size using Bayesian approach; is NULL by default and the lengh of MC chains in estimation step is then considered.
#' @param ... further arguments to be passed to or from other methods. They are ignored in this function.
#'
#' @return An object \code{curePredict} being a list with the following elements:
#'    \describe{
#'   \item{\code{out_pred}}{Dataframe with the (subject-specific) prediction results.}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters}
#'   \item{\code{data}}{the argument 'newdata' used}
#'   \item{\code{type}}{the argument 'type' used}
#'   \item{\code{state}}{the argument 'state' used}
#'    }
#'
#' @export
#'
#' @import mvtnorm joineR
#'
#' @author Antoine Barbieri and Catherine Legrand
#'
#' @seealso \code{\link{jointCureModel}}
#'
#' @examples
#'
#' ## For the exemple(s), use the data 'aids' from joineR package
#' data("aids", package = "joineR")
#'
#' ## estimation of the MCM for parameter initialisation
#' aids.id <- unique(aids[,c("id","time","death","drug","gender","prevOI","AZT")])
#' aids.id2 <- aids.id
#' aids.id2$drug <- as.numeric(aids.id$drug)-1
#' aids.id2$gender <- as.numeric(aids.id$gender)-1
#' aids.id2$prevOI <- as.numeric(aids.id$prevOI)-1
#' aids.id2$AZT <- as.numeric(aids.id$AZT)-1
#' smcure_out <- smcure(Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                      cureform=~ drug + gender + prevOI + AZT,
#'                      data = aids.id2,
#'                      model="ph")
#'
#' ## Estimation of the joint latent class cure model (JLCCM)
#' JLCCM <- jointCureModel(formFixed = CD4 ~ obstime + drug + gender + prevOI + AZT,
#'                         formRandom = ~ obstime,
#'                         timeVar= "obstime",
#'                         formLatency = Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                         formIncidence = ~ drug + gender + prevOI + AZT,
#'                         formID= ~ id,
#'                         IdVar = "id",
#'                         data = aids,
#'                         # model specifications
#'                         jointCureModel = "JLCCmodel",
#'                         survMod = "weibull-PH",
#'                         param = "shared-RE",
#'                         # prior options
#'                         n.iter = 1000,
#'                         n.burnin = 500,
#'                         Infprior_cure = TRUE,
#'                         smcure_out = smcure_out,
#'                         priorTau = 100,
#'                         # classification options
#'                         classif_trick = TRUE,
#'                         cov_prior = "inverse-gamma",
#'                         Sigma_d = TRUE)
#'
#' ## details of the estimated model
#' summary(JLCCM)
#'
#' ## prediction step with the JLCCM
#' # Frequentist approach
#' pred_JLCCM <- predict(object = JLCCM,
#'                       newdata = NULL,
#'                       type = "subject",
#'                       state = c("cure"),
#'                       level.intervalle=0.95,
#'                       Tsurv = NULL,
#'                       Yt = NULL,
#'                       yVar = "CD4",
#'                       idVar = "id",
#'                       MCMC = FALSE,
#'                       M = NULL)
#'
#' \dontrun{
#' # Bayesian approach
#' pred_JLCCM_MCMC <- predict(object = JLCCM,
#'                            newdata = NULL,
#'                            type = "subject",
#'                            state = c("cure"),
#'                            level.intervalle=0.95,
#'                            Tsurv=NULL,
#'                            Yt = NULL,
#'                            yVar = "CD4",
#'                            idVar = "id",
#'                            MCMC = TRUE,
#'                            M = 500)
#'}
#'
predict.JMcuR <- function(object,
                          newdata = NULL,
                          type = "subject",
                          state = "cure",
                          level.interval = 0.95,
                          Tsurv = NULL,
                          Yt = NULL,
                          yVar = "y",
                          idVar = "subject",
                          MCMC = FALSE,
                          M = NULL,
                          ...)
{

  # Extensions (see predict_class_tmp:
  # - predictions of the longitudinal outcomes at t*
  # - predictions of the probability to be event free at time t* given being event free at t<t*
  #
  # Modifications:
  # MCMC==TRUE: JLCCM given the matrix of variance of random effect specific by group !!
  #             FJCM too

  #--- initial conditions
  if (!inherits(object, "JMcuR"))
    stop("object does not belong to the 'JMcuR' class.\n")
  if(is.null(object$data) && is.null(newdata))
    stop("There is no data.\n")
  if(type=="profile" && is.null(newdata))
    stop("No profile prediction from dataset used to estimate the cure model.\n")
  if (object$survMod != "weibull-PH")
    stop("The survival model is not Weibul-PH.\n")
  if (!(object$jointCureModel %in%  c("FJCmodel","JLCCmodel","MCmodel")))
    stop("Only (joint) cure models belonging to 'FJCmodel', 'JLCCmodel', and 'MCmodel' are implemented for prediction step.\n")
  if (is.null(Yt) && type== "profile" && !(object$jointCureModel=="MCM"))
    stop("For profile prediction(s), it's necessary to give the associated 'Yt'.\n")
  if(!(is.null(Yt)) && state=="longitudinal")
    stop("To predict the longitudinal process at time 'Tsur', 'Yt' must be NULL and a set of longitudinal measurement must be avalaible.\n")
  if(object$jointCureModel=="MCmodel" && state=="longitudinal")
    stop("No longitudinal prediction can be done with MCmodel.\n")

  # Check if it is used...
  it_max=10


  #--- management data
  if(!is.null(newdata)) newdata <- cbind(newdata,Tsurv) else newdata <- object$data
  if(type=="profile")
    newdata <- cbind(newdata,Yt)
  # suprimer les 'time' qui sont supp?rieur ? Tsup
  if(nrow(newdata[all.vars(object$formID)])==1 && type=="subject")
    stop("Subject prediction need more one longitudinal measurement history. Otherwise use type='profile'.\n")
  if (is.null(newdata[[idVar]]))
    stop("'idVar' not in 'newdata.\n'")
  survMod <- object$survMod
  timeVar <- object$timeVar
  id <- as.numeric(unclass(newdata[[idVar]]))
  id <- id. <- match(id, unique(id))
  n.id <- length(unique(id))
  data.id <- newdata[!duplicated(id), ]
  idT <- data.id[[idVar]]
  idT <- match(idT, unique(idT))
  # incidence part (n.row=n.id)
  W <- model.matrix(object$formIncidence, model.frame(object$formIncidence, data = data.id[all.vars(object$formIncidence)]))
  # latency design vector for baseline regression
  Z <- model.matrix(delete.response(terms(object$formLatency)),
                    model.frame(delete.response(terms(object$formLatency)),
                                data = data.id[all.vars(delete.response(terms(object$formLatency)))]))

  if(!is.null(Tsurv)) Tsurv <- rep(Tsurv, n.id) else Tsurv <- as.vector(data.id[all.vars(object$formLatency)][, 1])

  #--- creation of objects
  codaFit <- do.call("rbind", object$codaFit)
  if(is.null(M))
    M <- nrow(codaFit)
  if (M > nrow(codaFit)) {
    warning("\n'M' cannot be set as greater than length Markov chain", nrow(codaFit))
    M <- nrow(codaFit)
    # res <- vector("list", M) # for the y
    success.rate <- matrix(FALSE, M, n.id)
  }
  codaFit <- codaFit[sample(nrow(codaFit), M), , drop = FALSE]
  n.sims <- M
  if(MCMC){
    postProba <- proba_uncured <- proba_cured <- S1_Tsurv <- matrix(NA, nrow = n.sims, ncol = n.id)
  }else{
    postProba <- proba_uncured <- proba_cured <- S1_Tsurv <- rep(NA, n.id)
  }
  out_pred <- as.data.frame(cbind(ID = data.id[, 1],
                                  delta = data.id[all.vars(object$formLatency)][, 2],
                                  time = Tsurv))


  #---------------------------------------------------------------------------------------------------#
  #--- Mixture cure model: Only latency and incidence model, without longitudinal process
  if(object$jointCureModel=="MCmodel"){

    if(MCMC){
      # Longer thant !MCMC, but the interest is about the CI on the quantities
      ind.thetas <- sapply(names(object$coefficients), grep, x = colnames(codaFit), fixed = TRUE)
      for(i in 1:n.id){
        #--- incidence probability
        proba_uncured[,i] <- 1/(1+exp(-codaFit[, ind.thetas$gamma]%*%(as.vector(W[i,]))))
        proba_cured[,i] <- 1-proba_uncured[,i]
        #--- latency
        S1_Tsurv[,i] <- exp(-exp(codaFit[, ind.thetas$alpha]%*%as.vector(Z[i,])) * (Tsurv[i]^codaFit[, ind.thetas$shape]))
        # posterior conditional probability of being cured (class membership)
        postProba[,i] <- (1-out_pred$delta[i]) * proba_cured[,i]/(proba_cured[,i]+proba_uncured[,i]*S1_Tsurv[,i])
      }
      # output
      out_pred <- cbind(out_pred, Proba.D0 = round(apply(postProba,2,median),3), Surv.D1 = round(apply(S1_Tsurv,2,median),3))
      out <- list(out_pred = out_pred)
      out$sims.list <- list(Latency.Tsurv = S1_Tsurv, Incidence.D1 = proba_uncured, Incidence.D0 = proba_cured, postProba.D0 = postProba)
    }else{
      # MCMC is FALSE, then no draw for the quantities, and then no CI
      #--- incidence probability
      proba_uncured <- 1/(1+exp(-as.matrix(W)%*%as.vector(object$coefficients$gamma)))
      proba_cured <- 1-proba_uncured
      #--- latency
      S1_Tsurv <- exp(-exp(as.matrix(Z)%*%as.vector(object$coefficients$alpha)) * (Tsurv^object$coefficients$shape))
      # posterior conditional probability of being cured (class membership)
      postProba <- (1-out_pred$delta) * proba_cured/(proba_cured+proba_uncured*S1_Tsurv)
      out_pred <- cbind(out_pred, Proba.D0 = round(postProba,3), Surv.D1 = round(S1_Tsurv,3))
      out <- list(out_pred = out_pred)
    }

  #---------------------------------------------------------------------------------------------------#
  #--- joint cure model including the longitudinal measurements
  }else{
    #--- Otherwise, prediction(s) including longitudinal process

    # longitudinal part (desing matrices, response variable)
    y <- model.response(model.frame(terms(object$formFixed),
                                    data = newdata))
    X <- model.matrix(delete.response(terms(object$formFixed)),
                      model.frame(delete.response(terms(object$formFixed)),
                                  data = newdata[all.vars(delete.response(terms(object$formFixed)))]))
    U <- model.matrix(object$formRandom, model.frame(object$formRandom, data = newdata[all.vars(object$formRandom)]))
    ncU <- ncol(U)
    ncX <- ncol(X)
    ncZ <- ncol(Z)
    data.tmp <- data.id
    data.tmp[, which(colnames(data.tmp)==object$timeVar)] <- Tsurv
    Xpred <- model.matrix(delete.response(terms(object$formFixed)),
                          model.frame(delete.response(terms(object$formFixed)),
                                      data = data.tmp[all.vars(delete.response(terms(object$formFixed)))]))
    Upred <- model.matrix(object$formRandom, model.frame(object$formRandom, data = data.tmp[all.vars(object$formRandom)]))
    # creation of objects
    if(MCMC){
      f1_Yt <- f0_Yt <- matrix(NA, nrow = n.sims, ncol = n.id)
    }else{
      f1_Yt <- f0_Yt <- rep(NA, n.id)
    }

    if("longitudinal" %in% state || object$jointCureModel=="FJCmodel"){
      Yhat_D0 <- Yhat_D1 <- rep(NA, n.id)
      RE_D0 <- RE_D1 <- matrix(NA, ncol=ncU, nrow = n.id)
      namesD0_uu <- namesD1_uu <- NULL
      for(uu in 1:ncU){
        namesD0_uu <- c(namesD0_uu, paste("RE",uu,".D0.",colnames(Upred)[uu],sep = ""))
        namesD1_uu <- c(namesD1_uu, paste("RE",uu,".D1.",colnames(Upred)[uu],sep = ""))
      }
      colnames(RE_D0) <- namesD0_uu
      colnames(RE_D1) <- namesD1_uu
    }

    #---------------------------------------------------------------------------------------------------#
    if(object$jointCureModel=="JLCCmodel"){

      codaFit <- cbind(codaFit, "tau" = sqrt(1 / codaFit[, "prec.tau2"]))
      codaFit <- codaFit[, -which(colnames(codaFit)=="prec.tau2")]
      colnames(codaFit)[which(grepl("beta\\[2,",colnames(codaFit)))] <- paste("beta.D0[",1:length(which(grepl("beta",colnames(codaFit)) & grepl("2,",colnames(codaFit)))),"]", sep = "")
      colnames(codaFit)[which(grepl("beta\\[1,",colnames(codaFit)))] <- paste("beta.D1[",1:length(which(grepl("beta",colnames(codaFit)) & grepl("1,",colnames(codaFit)))),"]", sep = "")
      # management of covariance matrix of random effects
      if(object$Sigma_d){
        colnames(codaFit)[which(grepl("prec.Sigma2\\[2,",colnames(codaFit)))] <- gsub("prec.Sigma2\\[2,", "prec.Sigma2.D0\\[",
                                                                                      colnames(codaFit)[which(grepl("prec.Sigma2\\[2,",colnames(codaFit)))])
        colnames(codaFit)[which(grepl("prec.Sigma2\\[1,",colnames(codaFit)))] <- gsub("prec.Sigma2\\[1,", "prec.Sigma2.D1\\[",
                                                                                      colnames(codaFit)[which(grepl("prec.Sigma2\\[1,",colnames(codaFit)))])
      }
      ind.thetas <- sapply(names(object$coefficients), grep, x = colnames(codaFit), fixed = TRUE)


      #--- MARGINAL prediction ---#
      if(type=="profile"){

        # management of variance of longitudinal process
        Uvar <- unique(Upred)
        if(!object$Sigma_d){
          var_RE <- NULL
          for(s in 1:n.sims){
            if(object$cov_prior=="wishart")
              mat_tmp  <- solve(matrix(codaFit[s,ind.thetas$Sigma2 ], ncU, ncU, TRUE))
            else
              mat_tmp  <- solve(diag(codaFit[s,ind.thetas$Sigma2 ]))
            var_RE <- c(var_RE, Uvar%*%mat_tmp%*%t(Uvar))
          }
          sd_LMM <- sqrt(var_RE + codaFit[, ind.thetas$tau]^2)
        }else{
          var_RE.D1 <- var_RE.D0 <- NULL
          for(s in 1:n.sims){
            if(object$cov_prior=="wishart"){
              mat_tmp.D1  <- solve(matrix(codaFit[s,ind.thetas$Sigma2.D1 ], ncU, ncU, TRUE))
              mat_tmp.D0  <- solve(matrix(codaFit[s,ind.thetas$Sigma2.D0 ], ncU, ncU, TRUE))
            }
            else{
              mat_tmp.D1  <- solve(diag(codaFit[s,ind.thetas$Sigma2.D1 ]))
              mat_tmp.D0  <- solve(diag(codaFit[s,ind.thetas$Sigma2.D0 ]))
            }
            var_RE.D1 <- c(var_RE, Uvar%*%mat_tmp.D1%*%t(Uvar))
            var_RE.D0 <- c(var_RE, Uvar%*%mat_tmp.D0%*%t(Uvar))
          }
          sd_LMM.D1 <- sqrt(var_RE.D1 + codaFit[, ind.thetas$tau]^2)
          sd_LMM.D0 <- sqrt(var_RE.D0 + codaFit[, ind.thetas$tau]^2)
        }

        if(length(Yt)==1)
          Yt <- rep(Yt, n.id)
        if(length(Yt)!=n.id)
          stop("'Yt' must be a scalar or a vector of n.id-length.")

        for(i in 1:n.id){
          #--- incidence probability
          proba_uncured[,i] <- 1/(1+exp(-codaFit[, ind.thetas$gamma]%*%(as.vector(W[i,]))))
          proba_cured[,i] <- 1-proba_uncured[,i]
          #--- latency contribution
          S1_Tsurv[,i] <- exp(-exp(codaFit[, ind.thetas$alpha]%*%as.vector(Z[i,])) * (Tsurv[i]^codaFit[, ind.thetas$shape]))
          #--- longitudinal contribution
          if(!object$Sigma_d){
            f1_Yt[,i] <- dnorm(Yt[i], mean = codaFit[, ind.thetas$beta.D1]%*%as.vector(Xpred[i,]), sd = sd_LMM)
            f0_Yt[,i] <- dnorm(Yt[i], mean = codaFit[, ind.thetas$beta.D0]%*%as.vector(Xpred[i,]), sd = sd_LMM)
          }else{
            f1_Yt[,i] <- dnorm(Yt[i], mean = codaFit[, ind.thetas$beta.D1]%*%as.vector(Xpred[i,]), sd = sd_LMM.D1)
            f0_Yt[,i] <- dnorm(Yt[i], mean = codaFit[, ind.thetas$beta.D0]%*%as.vector(Xpred[i,]), sd = sd_LMM.D0)
          }

          # output
          postProba[,i] <- (1-out_pred$delta[i]) * (f0_Yt[,i]*proba_cured[,i])/(f0_Yt[,i]*proba_cured[,i]+f1_Yt[,i]*proba_uncured[,i]*S1_Tsurv[,i])
        }
        # output
        out_pred <- cbind(out_pred, Yt, Proba.D0 = round(apply(postProba,2,median),3), Surv.D1 = round(apply(S1_Tsurv,2,median),3))
        out <- list(out_pred = out_pred)
        out$sims.list <- list(Latency_Tsurv = S1_Tsurv,
                              postProba.D0 = postProba,
                              Incidence.D1 = proba_uncured, Incidence.D0 = proba_cured,
                              f1_Yt = f1_Yt, f0_Yt = f0_Yt,
                              var_RE = var_RE,
                              sd_LMM = sd_LMM)
        if(!object$Sigma_d){
          out$sims.list$sd_LMM <- sd_LMM
        }else{
          out$sims.list$sd_LMM.D0 <- sd_LMM.D0
          out$sims.list$sd_LMM.D1 <- sd_LMM.D1
        }


      }else{
        #--- SUBJECT-SPECIFIC prediction  ---#

        mat0.RE <- mat1.RE <- array(NA, dim = c(n.id, n.sims, ncU))
        for(i in 1:n.id){
          #--- longitudinal contribution
          if(length(which(id==i))==1){
            U.id <- t(as.vector(U[which(id==i), ]))
            X.id <- t(as.vector(X[which(id==i), ]))
            y.id <- y[which(id==i)]
          }else{
            U.id <- U[which(id==i), ]
            X.id <- X[which(id==i), ]
            y.id <- y[which(id==i)]
          }
          if(MCMC){
            mat0.RE <- mat1.RE <- array(NA, dim = c(n.id, n.sims, ncU))
            #--- incidence probability
            proba_uncured[, i] <- 1/(1+exp(-codaFit[, ind.thetas$gamma]%*%(as.vector(W[i,]))))
            proba_cured[, i] <- 1-proba_uncured[,i]
            #--- latency contribution
            S1_Tsurv[,i] <- exp(-exp(codaFit[, ind.thetas$alpha]%*%as.vector(Z[i,])) * (Tsurv[i]^codaFit[, ind.thetas$shape]))
            for(s in 1:n.sims){
              if(!object$Sigma_d){
                tau2 <- codaFit[s, ind.thetas$tau]^2
                if(object$cov_prior=="wishart")
                  Sigma2 <- solve(matrix(codaFit[s,ind.thetas$Sigma2 ], ncU, ncU, TRUE))
                else
                  Sigma2 <- solve(diag(codaFit[s,ind.thetas$Sigma2 ]))
                if("longitudinal" %in% state){
                  ## conditional expectation of random effets following draw of parameters (cf. (G)LMM) given D=0 and D=1
                  mat0.RE[i,s,] <- (t(Sigma2%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                               +U.id%*%Sigma2%*%t(U.id))
                                      %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D0]))
                  )
                  mat1.RE[i,s,] <- (t(Sigma2%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                               +U.id%*%Sigma2%*%t(U.id))
                                      %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D1]))
                  )
                }
                ## marginal multivariate density
                f0_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D0]),
                                      sigma = (U.id)%*%as.matrix(Sigma2)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
                f1_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D1]),
                                      sigma = (U.id)%*%as.matrix(Sigma2)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))

              }else{ # class-specific covariance matrix of random effect
                tau2 <- codaFit[s, ind.thetas$tau]^2
                if(object$cov_prior=="wishart"){
                  Sigma2.D0 <- solve(matrix(codaFit[s,ind.thetas$Sigma2.D0 ], ncU, ncU, TRUE))
                  Sigma2.D1 <- solve(matrix(codaFit[s,ind.thetas$Sigma2.D1 ], ncU, ncU, TRUE))
                }
                else{
                  Sigma2.D0 <- solve(diag(codaFit[s,ind.thetas$Sigma2.D0 ]))
                  Sigma2.D1 <- solve(diag(codaFit[s,ind.thetas$Sigma2.D1 ]))
                }
                if("longitudinal" %in% state){
                  ## conditional espectation of random effets following draw of parameters (cf. (G)LMM) given D=0 and D=1
                  mat0.RE[i,s,] <- (t(Sigma2.D0%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                  +U.id%*%Sigma2.D0%*%t(U.id))
                                      %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D0]))
                  )
                  mat1.RE[i,s,] <- (t(Sigma2.D1%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                  +U.id%*%Sigma2.D1%*%t(U.id))
                                      %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D1]))
                  )
                }
                ## marginal multivariate density
                f0_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D0]),
                                      sigma = (U.id)%*%as.matrix(Sigma2.D0)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
                f1_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D1]),
                                      sigma = (U.id)%*%as.matrix(Sigma2.D1)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
              }

            }# end loop s

            ## posterior conditional probability of being cured
            postProba[, i] <- (1-out_pred$delta[i]) * (f0_Yt[, i]*proba_cured[, i])/(f0_Yt[, i]*proba_cured[, i]+f1_Yt[, i]*proba_uncured[, i]*S1_Tsurv[, i])

            ## output, longitudinal predictions (predict Y_i at time t given the group membership)
            if("longitudinal" %in% state){
              RE_D0[i,] <- apply(mat0.RE[i,,],2,mean)
              RE_D1[i,] <- apply(mat1.RE[i,,],2,mean)
              Yhat_D0[i] <- object$coefficients$beta.D0%*%as.vector(Xpred[i, ]) + RE_D0[i,]%*%as.vector(Upred[i, ])
              Yhat_D1[i] <- object$coefficients$beta.D1%*%as.vector(Xpred[i, ]) + RE_D1[i,]%*%as.vector(Upred[i, ])
            }

          }else{ #if MCMC is FALSE

            #--- incidence probability
            proba_uncured[i] <- 1/(1+exp(-object$coefficients$gamma%*%(as.vector(W[i,]))))
            proba_cured[i] <- 1-proba_uncured[i]
            #--- latency contribution
            S1_Tsurv[i] <- exp(-exp(object$coefficients$alpha%*%as.vector(Z[i,])) * (Tsurv[i]^object$coefficients$shape))

            #--- Variance components of longitudinal contribution
            tau2 <- object$coefficients$tau^2
            if("longitudinal" %in% state && length(which(id==i))!=1){
              ## conditional espectation of random effets following draw of parameters (cf. (G)LMM) given D=0 and D=1
              if(object$Sigma_d){
                RE_D0[i,] <- (t(object$coefficients$Sigma2.D0%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                                +U.id%*%object$coefficients$Sigma2.D0%*%t(U.id))
                                %*%(y.id-X.id%*%object$coefficients$beta.D0))
                )
                RE_D1[i,] <- (t(object$coefficients$Sigma2.D1%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                                +U.id%*%object$coefficients$Sigma2.D1%*%t(U.id))
                                %*%(y.id-X.id%*%object$coefficients$beta.D1))
                )
              }else{
                RE_D0[i,] <- (t(object$coefficients$Sigma2%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                             +U.id%*%object$coefficients$Sigma2%*%t(U.id))
                                %*%(y.id-X.id%*%object$coefficients$beta.D0))
                )
                RE_D1[i,] <- (t(object$coefficients$Sigma2%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                             +U.id%*%object$coefficients$Sigma2%*%t(U.id))
                                %*%(y.id-X.id%*%object$coefficients$beta.D1))
                )
              }
            }
            ## marginal multivariate density
            if(object$Sigma_d){
              f0_Yt[i] <- dmvnorm(y.id,
                                  mean = as.vector(X.id%*%object$coefficients$beta.D0),
                                  sigma = (U.id)%*%as.matrix(object$coefficients$Sigma2.D0)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
              f1_Yt[i] <- dmvnorm(y.id,
                                  mean = as.vector(X.id%*%object$coefficients$beta.D1),
                                  sigma = (U.id)%*%as.matrix(object$coefficients$Sigma2.D1)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
            }else{
              f0_Yt[i] <- dmvnorm(y.id,
                                  mean = as.vector(X.id%*%object$coefficients$beta.D0),
                                  sigma = (U.id)%*%as.matrix(object$coefficients$Sigma2)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
              f1_Yt[i] <- dmvnorm(y.id,
                                  mean = as.vector(X.id%*%object$coefficients$beta.D1),
                                  sigma = (U.id)%*%as.matrix(object$coefficients$Sigma2)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
            }


            ## posterior conditional probability of being cured
            postProba[i] <- (1-out_pred$delta[i]) * (f0_Yt[ i]*proba_cured[ i])/(f0_Yt[ i]*proba_cured[ i]+f1_Yt[ i]*proba_uncured[ i]*S1_Tsurv[ i])

            ## output, longitudinal preductions (predict Y_i at time t given the group membership)
            if("longitudinal" %in% state){
              Yhat_D0[i] <- object$coefficients$beta.D0%*%as.vector(Xpred[i, ]) + RE_D0[i,]%*%as.vector(Upred[i, ])
              Yhat_D1[i] <- object$coefficients$beta.D1%*%as.vector(Xpred[i, ]) + RE_D1[i,]%*%as.vector(Upred[i, ])
            }

          }# end !MCMC

        }#end loop for i


        # output
        if(MCMC){
          out_pred <- cbind(out_pred, Proba.D0 = round(apply(postProba,2,median),3), Surv.D1 = round(apply(S1_Tsurv,2,median),3))
        }else{
          out_pred <- cbind(out_pred, Proba.D0 = round(postProba,3), Surv.D1 = round(S1_Tsurv,3))
          mat0.RE <- mat1.RE <- NULL
        }

        if("longitudinal" %in% state)
          out_pred <- cbind(out_pred,
                            Yhat.D0 = Yhat_D0, RE_D0,
                            Yhat.D1 = Yhat_D1, RE_D1)
        out <- list(out_pred = out_pred)
        out$sims.list <- list(Latency_Tsurv = S1_Tsurv, postProba.D0 = postProba,
                              Incidence.D1 = proba_uncured, Incidence.D0 = proba_cured,
                              f1_Yt = f1_Yt, f0_Yt = f0_Yt,
                              mat0.RE = mat0.RE, mat1.RE = mat1.RE)

      }# end subject prediction



    #---  Full joint cure models --------------------------------------------------------------------------------------------#
    }else if(object$jointCureModel=="FJCmodel"){

      if(type=="profile"){
        #--- MARGINAL/PROFILE prediction  ---#
        stop("It would be preferable to consider 'JLCCmodel' with type='profile', otherwise use type='subject' if a set of longitudinal measurements is available by subject.\n")
        # see version september 2017 where a part is writing

      }else{
        #--- SUBJECT-SPECIFIC prediction  ---#

        param <- object$param

        if(!MCMC){
          #--- Estimation of the conditional probability of being cured (explicite solutions)

          # Take the estimated values of longitudinal parameter
          beta.D1 <- object$coefficients$beta.D1
          beta.D0 <- object$coefficients$beta.D0
          tau2 <- object$coefficients$tau^2
          if(object$Sigma_d && ncU>1){
            Sigma2.D0 <- object$coefficients$Sigma2.D0
            Sigma2 <- object$coefficients$Sigma2.D1
          }
          if(!object$Sigma_d && ncU>1){
            Sigma2 <- object$coefficients$Sigma2
            Sigma2.D0 <- object$coefficients$Sigma2
          }
          if(object$Sigma_d && ncU==1){
            Sigma2.D0 <- object$coefficients$sigma2[2]
            Sigma2 <- object$coefficients$sigma2[1]
          }
          if(!object$Sigma_d && ncU==1){
            Sigma2 <- object$coefficients$sigma2
            Sigma2.D0 <- object$coefficients$sigma2
          }
          # Take the estimated values of incidence parameters
          gamma <- object$coefficients$gamma
          # Take the estimated values of latency parameters
          shape <- object$coefficients$shape
          alpha <- object$coefficients$alpha[1:ncZ]
          if(param=="td-value"){
            alphaL <- object$coefficients$alpha[1+ncZ]
            alphaRE <- NULL
          }
          if(param=="shared-RE"){
            alphaL <- NULL
            alphaRE <- object$coefficients$alpha[(1+ncZ):(ncU+ncZ)]
          }
          list.thetas <- list(beta = beta.D1, beta.D0 = beta.D0,
                              tau2 = tau2, tau2.D0 = tau2,
                              Sigma2 = Sigma2, Sigma2.D0 = Sigma2.D0,
                              gamma = gamma,
                              alpha = alpha, alphaL = alphaL, alphaRE = alphaRE, shape = shape)
          list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
          thetas <- unlist(as.relistable(list.thetas))

          for(i in seq_len(n.id)){

            #--- longitudinal contribution
            if(length(which(id==i))==1){
              U.id <- t(as.vector(U[which(id==i), ]))
              X.id <- t(as.vector(X[which(id==i), ]))
              y.id <- y[which(id==i)]
            }else{
              U.id <- U[which(id==i), ]
              X.id <- X[which(id==i), ]
              y.id <- y[which(id==i)]
            }

            #--- incidence probability
            proba_uncured[i] <- 1/(1+exp(-list.thetas$gamma%*%(as.vector(W[i,]))))
            proba_cured[i] <- 1-proba_uncured[i]

            #--- marginal multivariate density
            f0_Yt[i] <- dmvnorm(y.id,
                                mean = as.vector(X.id%*%list.thetas$beta.D0),
                                sigma = (U.id)%*%as.matrix(list.thetas$Sigma2.D0)%*%t(U.id)+diag(list.thetas$tau2.D0,ncol=length(y.id),nrow=length(y.id)))
            f1_Yt[i] <- dmvnorm(y.id,
                                mean = as.vector(X.id%*%list.thetas$beta),
                                sigma = (U.id)%*%as.matrix(list.thetas$Sigma2)%*%t(U.id)+diag(list.thetas$tau2,ncol=length(y.id),nrow=length(y.id)))

            #--- Explicite normal distribution of random effect for (b_i|Y,theta)
            if("cure" %in% state){
              #--- Empiriacl bayes random effects (b_i|Y_i=y_i,D_i=d_i which is an explicite normal distribution)
              RE_D0[i,] <- (t(list.thetas$Sigma2.D0%*%t(U.id)%*%solve(diag(list.thetas$tau2.D0,ncol=length(y.id),nrow=length(y.id))
                                                                      +as.matrix(U.id)%*%list.thetas$Sigma2.D0%*%t(U.id))
                              %*%(y.id-X.id%*%list.thetas$beta.D0))
              )
              RE_D1[i,] <- (t(list.thetas$Sigma2%*%t(U.id)%*%solve(diag(list.thetas$tau2,ncol=length(y.id),nrow=length(y.id))
                                                                   +as.matrix(U.id)%*%list.thetas$Sigma2%*%t(U.id))
                              %*%(y.id-X.id%*%list.thetas$beta))
              )
            }

            #--- Latency estimation
            etaBaseline <- as.vector(list.thetas$alpha%*%Z[i,])
            # given shared random effects
            if(param=="shared-RE")
              S1_Tsurv[i] <- exp(-exp(etaBaseline + RE_D1[i,]%*%list.thetas$alphaRE)*Tsurv[i]^(list.thetas$shape))
            # given shared current value which need integral approximation
            if(param=="td-value"){
              #--- td-value case using Gauss_Kronrod quadrature
              if (!object$timeVar %in% names(object$data))
                stop("\n'timeVar' does not correspond to one of the columns in formulas")
              # maybe need gaussKronrod() if not source at the beginning
              wk <- gaussKronrod()$wk
              sk <- gaussKronrod()$sk
              K <- length(sk)
              P <- Tsurv[i]/2
              st <- outer(P, sk + 1)
              id.GK <- rep(i, each = K)
              data.id2 <- data.tmp[id.GK, ]
              data.id2[[timeVar]] <- c(t(st))
              mfX <- model.frame(object$formFixed, data = data.id2)
              mfU <- model.frame(object$formRandom, data = data.id2)
              Xs <- model.matrix(object$formFixed, mfX)
              Us <- model.matrix(object$formRandom, mfU)
              # if (one.RE)
              #   Us <- cbind(Us, rep(0, nrow(Us)))
              h_k <- as.vector(wk) * list.thetas$shape * as.vector(st^(shape - 1)) * exp(alphaL * as.vector(Xs%*%list.thetas$beta + Us%*%RE_D1[i,]))
              S1_Tsurv[i] <- exp(-exp(etaBaseline) * P * sum(h_k))
            }
            if("cure" %in% state){
              #--- posterior conditional probability of being cured
              postProba[i] <- (1-out_pred$delta[i]) * (f0_Yt[ i]*proba_cured[i])/(f0_Yt[ i]*proba_cured[i]+f1_Yt[i]*proba_uncured[i]*S1_Tsurv[i])
            }
          }

        }else{
          #---- if(MCMC)

          cat("It can take time! Have a coffee break ;)")

          # use the Markov Chain of parameters from the object of class JMcure
          codaFit <- cbind(codaFit, "tau" = sqrt(1 / codaFit[, "prec.tau2"]))
          codaFit <- codaFit[, -which(colnames(codaFit)=="prec.tau2")]
          colnames(codaFit)[which(grepl("beta\\[2,",colnames(codaFit)))] <- paste("beta.D0[",1:length(which(grepl("beta",colnames(codaFit)) & grepl("2,",colnames(codaFit)))),"]", sep = "")
          colnames(codaFit)[which(grepl("beta\\[1,",colnames(codaFit)))] <- paste("beta.D1[",1:length(which(grepl("beta",colnames(codaFit)) & grepl("1,",colnames(codaFit)))),"]", sep = "")
          # management of covariance matrix of random effects
          if(object$Sigma_d){
            colnames(codaFit)[which(grepl("prec.Sigma2\\[2,",colnames(codaFit)))] <- gsub("prec.Sigma2\\[2,", "prec.Sigma2.D0\\[",
                                                                                          colnames(codaFit)[which(grepl("prec.Sigma2\\[2,",colnames(codaFit)))])
            colnames(codaFit)[which(grepl("prec.Sigma2\\[1,",colnames(codaFit)))] <- gsub("prec.Sigma2\\[1,", "prec.Sigma2.D1\\[",
                                                                                          colnames(codaFit)[which(grepl("prec.Sigma2\\[1,",colnames(codaFit)))])
          }
          ind.thetas <- sapply(c(names(object$coefficients),"alphaL","alphaRE"), grep, x = colnames(codaFit), fixed = TRUE)
          # separate the different coefficient index
          if(param=="td-value"){
            ind.thetas$alpha <- setdiff(ind.thetas$alpha, ind.thetas$alphaL)
            ind.thetas$alphaRE <- NULL
          }
          if(param=="shared-RE"){
            ind.thetas$alpha <- setdiff(ind.thetas$alpha, ind.thetas$alphaRE)
            ind.thetas$alphaL <- NULL
          }

          # create the matrix of random effects
          mat0.RE <- mat1.RE <- array(NA, dim = c(n.id, n.sims, ncU))

          for(i in seq_len(n.id)){

            #--- longitudinal contribution
            if(length(which(id==i))==1){
              U.id <- t(as.vector(U[which(id==i), ]))
              X.id <- t(as.vector(X[which(id==i), ]))
              y.id <- y[which(id==i)]
            }else{
              U.id <- U[which(id==i), ]
              X.id <- X[which(id==i), ]
              y.id <- y[which(id==i)]
            }

            #--- incidence probability
            proba_uncured[, i] <- 1/(1+exp(-codaFit[, ind.thetas$gamma]%*%(as.vector(W[i,]))))
            proba_cured[, i] <- 1-proba_uncured[,i]

            #--- latency contribution at baseline
            etaBaseline <- as.vector(codaFit[, ind.thetas$alpha]%*%as.vector(Z[i,]))

            # To approximate the integrale which depends of the i-th subject
            if(param=="td-value"){
              # Initialization for the numerical intergration using Gauss-Kronrod quadrature
              if (!object$timeVar %in% names(object$data))
                stop("\n'timeVar' does not correspond to one of the columns in formulas")
              wk <- gaussKronrod()$wk
              sk <- gaussKronrod()$sk
              K <- length(sk)
              P <- Tsurv[i]/2
              st <- outer(P, sk + 1)
              id.GK <- rep(i, each = K)
              data.id2 <- data.tmp[id.GK, ]
              data.id2[[timeVar]] <- c(t(st))
              mfX <- model.frame(object$formFixed, data = data.id2)
              mfU <- model.frame(object$formRandom, data = data.id2)
              Xs <- model.matrix(object$formFixed, mfX)
              # if (one.RE)
              #   Us <- cbind(Us, rep(0, nrow(Us)))
              Us <- model.matrix(object$formRandom, mfU)
            }

            # beginning of the M draws
            for(s in 1:n.sims){

              # common matrix of random effect
              if(!object$Sigma_d){
                tau2 <- codaFit[s, ind.thetas$tau]^2
                if(object$cov_prior=="wishart")
                  Sigma2 <- solve(matrix(codaFit[s,ind.thetas$Sigma2 ], ncU, ncU, TRUE))
                else
                  Sigma2 <- solve(diag(codaFit[s,ind.thetas$Sigma2 ]))
                  ## conditional espectation of random effets following draw of parameters (cf. (G)LMM) given D=0 and D=1
                mat0.RE[i,s,] <- (t(Sigma2%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                             +U.id%*%Sigma2%*%t(U.id))
                                    %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D0]))
                )
                mat1.RE[i,s,] <- (t(Sigma2%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                             +U.id%*%Sigma2%*%t(U.id))
                                    %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D1]))
                )
                ## marginal multivariate density
                f0_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D0]),
                                      sigma = (U.id)%*%as.matrix(Sigma2)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
                f1_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D1]),
                                      sigma = (U.id)%*%as.matrix(Sigma2)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))

              }else{  # class-specific covariance matrix of random effect

                tau2 <- codaFit[s, ind.thetas$tau]^2
                if(object$cov_prior=="wishart"){
                  Sigma2.D0 <- solve(matrix(codaFit[s,ind.thetas$Sigma2.D0 ], ncU, ncU, TRUE))
                  Sigma2.D1 <- solve(matrix(codaFit[s,ind.thetas$Sigma2.D1 ], ncU, ncU, TRUE))
                }
                else{
                  Sigma2.D0 <- solve(diag(codaFit[s,ind.thetas$Sigma2.D0 ]))
                  Sigma2.D1 <- solve(diag(codaFit[s,ind.thetas$Sigma2.D1 ]))
                }
                ## conditional espectation of random effets following draw of parameters (cf. (G)LMM) given D=0 and D=1
                mat0.RE[i,s,] <- (t(Sigma2.D0%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                +U.id%*%Sigma2.D0%*%t(U.id))
                                    %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D0]))
                )
                mat1.RE[i,s,] <- (t(Sigma2.D1%*%t(U.id)%*%solve(diag(tau2,ncol=length(y.id),nrow=length(y.id))
                                                                +U.id%*%Sigma2.D1%*%t(U.id))
                                    %*%(y.id-X.id%*%codaFit[s,ind.thetas$beta.D1]))
                )
                ## marginal multivariate density
                f0_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D0]),
                                      sigma = (U.id)%*%as.matrix(Sigma2.D0)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
                f1_Yt[s,i] <- dmvnorm(y.id,
                                      mean = as.vector(X.id%*%codaFit[s,ind.thetas$beta.D1]),
                                      sigma = (U.id)%*%as.matrix(Sigma2.D1)%*%t(U.id)+diag(tau2,ncol=length(y.id),nrow=length(y.id)))
              }

              # latency contribution given shared current value which need integral approximation
              if(param=="td-value"){
                h_k <- as.vector(wk) * codaFit[s,ind.thetas$shape] * as.vector(st^(codaFit[s,ind.thetas$shape] - 1)) * exp(codaFit[s,ind.thetas$alphaL] * as.vector(Xs%*%codaFit[s,ind.thetas$beta.D1] + Us%*%mat1.RE[i, s, ]))
                S1_Tsurv[s, i] <- exp(-exp(etaBaseline[s]) * P * sum(h_k))
              }

            }# end loop s

            # latency contribution given shared random effects
            if(param=="shared-RE")
              S1_Tsurv[, i] <- exp(-exp(etaBaseline + as.vector(apply(codaFit[, ind.thetas$alphaRE]*mat1.RE[i, , ],1,sum)))*Tsurv[i]^(codaFit[, ind.thetas$shape]))

            # Posterior probability of being cured for censored subjects
            postProba[, i] <- (1-out_pred$delta[i]) * (f0_Yt[, i]*proba_cured[, i])/(f0_Yt[, i]*proba_cured[, i]+f1_Yt[, i]*proba_uncured[, i]*S1_Tsurv[, i])

          }# end loop i

      }# end if(MCMC)

        # output
        if(MCMC){
          out_pred <- cbind(out_pred, Proba.D0 = round(apply(postProba,2,median),3), Surv.D1 = round(apply(S1_Tsurv,2,median),3))
        }else{
          out_pred <- cbind(out_pred, Proba.D0 = round(postProba,3), Surv.D1 = round(S1_Tsurv,3))
          mat0.RE <- mat1.RE <- NULL
        }

        out <- list(out_pred = out_pred)
        out$sims.list <- list(Latency_Tsurv = S1_Tsurv,
                              postProba.D0 = postProba,
                              Incidence.D1 = proba_uncured, Incidence.D0 = proba_cured,
                              f1_Yt = f1_Yt, f0_Yt = f0_Yt,
                              mat0.RE = mat0.RE, mat1.RE = mat1.RE)

      }# end individual predictions with FJCMs

      #---------------------------------------------------------------------------------------------------#
    }else if(object$jointCureModel %in% c("JSECmodel","FJSECmodel")){ # individual prediction
      stop("Prediction is not available for 'JSECmodel' or 'FJSECmodel'.\n")
    }
  }

  # whatever the model
  # out$coefficients <- mean(postProba)
  # out$StDev <- sd(as.vector(postProba))
  # alpha <- 1-level.interval
  # out$CIs <- quantile(postProba, probs = c(alpha/2, 1-alpha/2))
  # out$Yt
  out$data <- newdata
  out$type <- type
  out$state <- state
  #out$object <- object
  class(out) <- "curePredict"
  out

}
