#' Estimation of a classical mixture cure model (survival model including a cure fraction)
#'
#' Function using JAGS to estimate the mixture cure model (MCM)[Farewell, 1982].
#'
#' @param formLatency survival formula as formula in survival package for latency submodel
#' @param formIncidence formula specifying covariate in incidence submodel
#' @param formID formula specifying the ID variable (e.g. = ~ subject), if NULL each row is a statistic unit
#' @param data dataset of observed variables
#' @param survMod form of survival submodel (only "weibull-PH" is available until now)
#' @param inf_prior if TRUE, the estimation is built considering informative prior from survival (coxph function) or smcure package. If FALSE, vague priors are used.
#' @param smcure_out object from smcure package estimating model and used for data-driven priors, otherwise is NULL by default
#' @param classif_trick if TRUE, subjects are cured if t_obs>max(t_obs(delta==1), otherwise D is only know for subjects having developed the event Taylor's rule about the cure status, D==0 if T_obs>max(t_obs[which(delta==1)])
#' @param n.chains the number of parallel chains for the model; default is 3
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is NULL
#' @param parallel if TRUE, the program is parallelized
#' @param C positive integer for the zero trick used to define the likelihood in JAGS, default is 1000
#' @param priorTau variance by default for vague prior distribution
#' @param save_jagsUI If TRUE (by default), the output of jagsUI package is return by the function
#' @param out_data Boolean such as TRUE if you want the data of different submodels in output or FALSE otherwise
#' @param verbose If set to FALSE, all text output in the console will be suppressed as the function runs (including most warnings).
#'
#' @return A \code{JMcuR} object which is a list with the following elements:
#'    \describe{
#'   \item{\code{mean}}{list of posterior mean for each parameter}
#'   \item{\code{median}}{list of posterior median for each parameter}
#'   \item{\code{modes}}{list of posterior mode for each parameter}
#'   \item{\code{StErr}}{list of standard error for each parameter}
#'   \item{\code{StDev}}{list of standard deviation for each parameter}
#'   \item{\code{Rhat}}{Gelman and Rubin diagnostic for all parameters}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95 for each parameters excepted for covariance parameters in covariance matrix of random effects. Otherwise, use save_jagsUI=TRUE to have the associated quantiles.}
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters and random effects}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{out_jagsUI}}{only if \code{save_jagsUI=TRUE} in argument: list including posterior mean, median, quantiles (2.5%, 25%, 50%, 75%, 97.5%), standart deviation for each parameter and each random effect.
#'   Moreover, this list also returns the MCMC draws, the Gelman and Rubin diagnostics (see output of jagsUI objects)}
#'  }
#'
#' @export
#'
#' @import jagsUI survival smcure
#'
#' @author Antoine Barbieri and Catherine Legrand
#'
#' @references Barbieri A and Legrand C. (2019). \emph{Joint longitudinal and time-to-event cure models for the assessment of being cured}. Statistical Methods in Medical Research. doi: 10.1177/0962280219853599.
#' @references Taylor JMG. (1995) \emph{Semi-Parametric Estimation in Failure Time Mixture Models}. Biometrics 1995; 51(3): 899-907.
#' @references Farewell VT. \emph{The use of mixture models for the analysis of survival data with long-term survivors}. Biometrics 1982; 38(4):1041–1046.
#'
#' @examples
#'
#' # For the exemple(s), use the data 'aids' from joineR package
#' data("aids", package = "joineR")
#'
#' ## Estimation of the mixture cure model (MCM)
#' MCM1 <- mCureModel(formLatency = Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                    formIncidence = ~ drug + gender + prevOI + AZT,
#'                    formID = ~ id,
#'                    data = aids,
#'                    # model specifications
#'                    survMod = "weibull-PH",
#'                    # prior options
#'                    n.iter = 1000,
#'                    n.burnin = 500,
#'                    smcure_out = NULL,
#'                    priorTau = 100,
#'                    # classification options
#'                    classif_trick = TRUE)
#'
#' # details of the estimated model
#' summary(MCM1)
#'
#' \dontrun{
#'
#' # estimation of the MCM with parameter initialisation
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
#' MCM2 <- mCureModel(formLatency = Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                    formIncidence = ~ drug + gender + prevOI + AZT,
#'                    data = aids.id,
#'                    # model specifications
#'                    survMod = "weibull-PH",
#'                    # prior options
#'                    smcure_out = smcure_out,
#'                    priorTau = 100,
#'                    # classification options
#'                    classif_trick = TRUE)
#'
#' # details of the estimated model
#' summary(MCM2)
#' }
#'
mCureModel <- function(formLatency,
                       formIncidence,
                       formID = NULL,
                       data,
                       survMod = "weibull-PH",
                       inf_prior = TRUE,
                       smcure_out = NULL,
                       classif_trick = TRUE,
                       n.chains = 3,
                       n.iter = 10000,
                       n.burnin = 5000,
                       n.thin = 1,
                       n.adapt = NULL,
                       C = 1000,
                       priorTau = 100,
                       save_jagsUI = TRUE,
                       out_data = T,
                       parallel = FALSE,
                       verbose = TRUE)
{

  # lag = 0
  # quiet = FALSE

  #### ---- cure part ----

  ## data
  if(!is.null(formID)){
    data_cure <- data[c(all.vars(formID),all.vars(formLatency),all.vars(formIncidence))]
    data_cure <- unique(data_cure)
  }else{
    data_cure <- data[c(all.vars(formLatency),all.vars(formIncidence))]
  }
  Time <- data_cure[all.vars(formLatency)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  delta <- data_cure[all.vars(formLatency)][, 2]   # vector of event indicator (delta)
  nTime <- length(Time)                            # number of subject having Time
  # design matrices
  mfZ <- model.frame(formLatency, data = data_cure)
  Z <- model.matrix(formLatency, mfZ)    # or: Z <- if (survMod == "weibull-PH") {if (is.null(Z)) cbind(rep(1, nT), rep(0, nT)) else cbind(1, Z)}
  mfW <- model.frame(formIncidence, data = data_cure)
  W <- model.matrix(formIncidence, mfW)  # or: W <- cbind(1,data_cure[all.vars(formIncidence)])
  zeros <- numeric(nTime)

  # partially latent variable D
  D <- delta
  D[which(D==0)] <- NA
  if(classif_trick){
    # time_threshold <-  max(range(boxplot(data_cure$t_obs[which(data_cure$delta==1)],plot=FALSE)$stats))
    time_threshold <- max(Time[which(delta==1)])
    D[which(delta==0 & Time>time_threshold)] <- 0
  }

  initial.values <- list(shape = 1)

  if(inf_prior && is.null(smcure_out)){
    cat("> Initiation of survival parameter values using 'survival' package. \n")
    tmp_model <- survival::coxph(formLatency,
                                 data = data_cure,
                                 x = TRUE)
    # prior parameters for latency model
    priorMean.alpha <- c(0, tmp_model$coefficients)
    priorTau.alpha <- diag(c(1/priorTau, 1/(priorTau*diag(tmp_model$var))))
    # prior parameters for incidence model
    priorMean.gamma <- as.numeric(rep(0,ncol(W)))
    priorTau.gamma <- diag(rep(1/priorTau,length(priorMean.gamma)))
    # initialisation values of survival parameters
    if(survMod=="weibull")
      initial.values$alpha <- c(0, tmp_model$coefficients)
    initial.values$gamma <- as.numeric(rep(0,ncol(W)))
  }

  if(class(smcure_out)=="smcure"){
    # prior parameters for latency model
    priorMean.alpha <- as.numeric(c(0,smcure_out$beta))
    priorTau.alpha <- diag(1/c(1,smcure_out$beta_var)/priorTau)
    # prior parameters for incidence model
    priorMean.gamma <- as.numeric(smcure_out$b)
    priorTau.gamma <- diag(1/smcure_out$b_var/priorTau)
    # initialisation values of survival parameters
    if(survMod=="weibull")
      initial.values$alpha <- c(0, tmp_model$coefficients)
    initial.values$gamma <- priorMean.gamma
  }
  if(!inf_prior){
    # prior parameters for latency model
    priorMean.alpha <- as.numeric(rep(0,ncol(Z)))
    priorTau.alpha <- diag(rep(1/priorTau,length(priorMean.alpha)))
    # prior parameters for incidence model
    priorMean.gamma <- as.numeric(rep(0,ncol(W)))
    priorTau.gamma <- diag(rep(1/priorTau,length(priorMean.gamma)))
    # initialisation values of survival parameters
    if(survMod=="weibull")
      initial.values$alpha <- c(0, tmp_model$coefficients)
    initial.values$gamma <- priorMean.gamma
  }

  # jags data
  jags.data <- list(C = C, zeros = numeric(nTime),
                    Time = Time, delta=delta,
                    N = nTime, ncZ = ncol(Z), ncW = ncol(W),
                    class = D,
                    Z = Z,
                    W = W,
                    priorMean.gamma = priorMean.gamma, priorTau.gamma = priorTau.gamma,
                    priorMean.alpha = priorMean.alpha, priorTau.alpha = priorTau.alpha,
                    priorA.shape = 1/priorTau, priorB.shape = 1/priorTau)

  param <- "shared-RE"
  one.RE <- TRUE

  #--- Model specification
  model <- switch(paste(survMod, sep = "/"),
                  # Weibull,
                  `weibull-PH` = MCM.Weib
  )

  #--- Parameter to save specification
  parms_to_save <- c("alpha", "gamma")
  parms_to_save <- switch(paste(survMod, sep = "/"),
                          ### Weibull
                          `weibull-PH` = c(parms_to_save, "shape","class")
                          )

  working.directory = getwd()
  write.model.jags(model=model,
                   intitled=file.path(working.directory,"JagsModel.txt"),
                   Data=jags.data,
                   jointCureModel="MCmodel",
                   param="shared-RE",
                   one.RE=TRUE)

  initial.values <- initial.values[names(initial.values) %in% parms_to_save]

  #---- Call to JAGS to estimate the model
  if(!verbose)
    cat("> JAGS is running, it can take time. Keep calm, have a cafee ;) \n")
  out_jags = jagsUI::jags(data = jags.data,
                          parameters.to.save = parms_to_save,
                          model.file = "JagsModel.txt",
                          n.chains = n.chains,
                          parallel = parallel,
                          n.adapt = n.adapt,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          verbose = verbose,
                          DIC = F#, codaOnly = c("class")
  )
  file.remove("JagsModel.txt")

  #--- Management of the outputs
  out <- list(data = data)
  out$control <- list(n.chains = n.chains,
                      parallel = parallel,
                      n.adapt = n.adapt,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin,
                      classif_trick = classif_trick,
                      formIncidence = formIncidence,
                      formLatency = formLatency,
                      formID = formID,
                      survMod = survMod,
                      n=nTime)

  # sims.list output
  out$sims.list <- out_jags$sims.list
  out$sims.list$class <- NULL

  # posterior classification of susceptible suject
  out$posterior_class_uncured <- out_jags$q05$class

  # median : Posterior median of parameters (if mean, you can use mean instead of q50)
  out$median <- out_jags$q50
  out$median$class <- NULL

  # mean : Posterior mean of parameters (if mean, you can use mean instead of q50)
  out$coefficients <- out_jags$mean
  out$coefficients$class <- NULL

  # modes of parameters
  out$modes <- lapply(out$sims.list, function(x) {
    m <- function(x) {
      d <- density(x, bw = "nrd", adjust = 3, n = 1000)
      d$x[which.max(d$y)]
    }
    if (is.matrix(x))
      as.array(apply(x, 2, m))
    else{
      if(is.array(x))
        apply(x, c(2,3), m)
      else m(x)
    }
  })

  # standard error of parameters
  out$StErr <- lapply(out$sims.list, function(x) {
    f <- function(x) {
      acf.x <- drop(acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
      acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
      ess <- length(x)/(1 + 2 * sum(acf.x))
      sqrt(var(x)/ess)
    }
    if (is.matrix(x))
      as.array(apply(x, 2, f))
    else{
      if(is.array(x))
        apply(x, c(2,3), f)
      else f(x)
    }
  })

  # standard deviation of parameters
  out$StDev <- out_jags$sd
  out$StDev$class <- NULL

  # Rhat : Gelman & Rubin diagnostic
  out$Rhat <- out_jags$Rhat
  out$Rhat$class <- NULL

  # names
  names(out$coefficients$alpha) <-
    names(out$median$alpha) <-
    names(out$Rhat$alpha) <-
    names(out$modes$alpha) <-
    names(out$StErr$alpha) <-
    names(out$StDev$alpha) <- colnames(Z)

  names(out$coefficients$gamma) <-
    names(out$median$gamma) <-
    names(out$Rhat$gamma) <-
    names(out$modes$gamma) <-
    names(out$StErr$gamma) <-
    names(out$StDev$gamma) <- colnames(W)

  # credible intervalles
  out$CIs$alpha <- cbind(as.vector(t(out_jags$q2.5$alpha)),
                        as.vector(t(out_jags$q97.5$alpha)))
  rownames(out$CIs$alpha) <- colnames(Z)
  colnames(out$CIs$alpha) <- c("2.5%", "97.5%")
  out$CIs$gamma <- cbind(as.vector(t(out_jags$q2.5$gamma)),
                         as.vector(t(out_jags$q97.5$gamma)))
  rownames(out$CIs$gamma) <- colnames(W)
  colnames(out$CIs$gamma) <- c("2.5%", "97.5%")

  out$CIs$shape <- c(out_jags$q2.5$shape,
                     out_jags$q97.5$shape)
  names(out$CIs$shape) <- c("2.5%", "97.5%")

  # save jags output if requires
  out$jointCureModel <- "MCmodel"
  if(save_jagsUI)
    out$out_jagsUI <- out_jags

  class(out) <- "JMcuR"
  out

}#end function
