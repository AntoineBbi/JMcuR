#' Estimation of the joint longitudinal time-to-event models including a cure fraction
#'
#' Function using JAGS to estimate the joint longitudinal survival models with cure fraction: mixutre cure model (MCmodel) [Farewell, 1982], joint latent class cure model (JLCCmodel),  full joint cure model (FJCmodel) [Law et al. 2002, YU et al. 2004,2008]
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without response variable
#' @param formLatency survival formula as formula in survival package for latency submodel
#' @param formIncidence formula specifying covariate in incidence submodel
#' @param formID formula specifying the ID variable (e.g. = ~ subject)
#' @param IdVar string specify the names of ID variable (Id subject)
#' @param data dataset of observed variables
#' @param timeVar string specify the names of time variable (time of repeated measurements)
#' @param jointCureModel joint cure model to estimate including: "MCmodel","JLCCmodel","FJCmodel"
#' @param survMod form of survival submodel (only "weibull-PH" is available until now)
#' @param param shared association for "FJCmodel": "shared-RE" (default) or "td-value"
#' @param Infprior_cure most of prior distributions are data-driven if 'TRUE', otherwise vague prior distribution are considered
#' @param smcure_out estimated model from smcure package, otherwise is NULL by default
#' @param classif_trick if TRUE, subjects are cured if t_obs>max(t_obs(delta==1), otherwise D is only know for subjects having developed the event Taylor's rule about the cure status, D==0 if T_obs>max(t_obs[which(delta==1)])
#' @param cov_prior ="wishart" prior distribtion for random effect covariance matrix or ="inverse-gamma" for independent prior distribtions with random effect variance matrix
#' @param Sigma_d if TRUE, the random effect matrix is class-specific, otherwise is the matrix is defined for the whole population
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter nteger specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is 5000
#' @param C positive integer for the zero trick used to define the likelihood in JAGS, default is 1
#' @param priorTau variance by default for vague prior distribution
#' @param out_data Boolean such as TRUE if you want the data of different submodels in output or FALSE otherwise
#'
#' @return A 'JMcuR' object which is a list with the following elements:
#'    \describe{
#'   \item{\code{Coefficients}}{list of posterior mean of each parameter}
#'   \item{\code{Modes}}{ list of posterior modes compute from the posterior samples of the parameters}
#'   \item{\code{Sd}}{list of Standard deviation of each parameters}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters}
#'    }
#'
#' @export
#'
#' @import smcure rjags lcmm coda
#'
#' @examples
#'
#' ## estimation of the MCM for parameter initialisation
#' aids.id <- unique(aids[,c("id","time","death","drug","gender","prevOI","AZT")])
#' table(aids.id$death)
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
#'                         jointCureModel = c("JLCCmodel"),
#'                         survMod = c("weibull-PH"),
#'                         param = c("shared-RE"),
#'                         # prior options
#'                         Infprior_cure = TRUE,
#'                         smcure_out = smcure_out,
#'                         priorTau = 100,
#'                         # classification options
#'                         classif_trick = TRUE,
#'                         cov_prior = c("inverse-gamma"),
#'                         Sigma_d = TRUE)
#'
#' ## details of the estimated model
#' summary(JLCCM)
#'
#' ## Estimation of the full joint cure model given shared curent value (FJCM_A1)
#' FJCM_A1 <- cureJMbayes(formFixed = CD4 ~ obstime,
#'                        formRandom = ~ obstime,
#'                        timeVar= "obstime",
#'                        formLatency = Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                        formIncidence = ~ drug + gender + prevOI + AZT,
#'                        formID= ~ id,
#'                        IdVar = "id",
#'                        data = aids,
#'                        # model specifications
#'                        jointCureModel = c("FJCmodel"),
#'                        survMod = c("weibull-PH"),
#'                        param = c("td-value"),
#'                        # prior options
#'                        Infprior_cure = TRUE,
#'                        smcure_out = smcure_out,
#'                        priorTau = 100,
#'                        # classification options
#'                        classif_trick = TRUE,
#'                        cov_prior = c("inverse-gamma"),
#'                        Sigma_d = TRUE,
#'                        # MCMC options
#'                        n.chains = 1, n.iter = 10000, n.burnin = 5000, n.thin = 1, n.adapt = 5000, quiet = FALSE, C = 1,
#'                        # out option
#'                        out_data=T)
#'
#' ## details of the estimated model
#' summary(FJCM_A1)
#'
#' ## Estimation of the full joint cure model given shared curent value (FJCM_A1)
#' FJCM_A2 <- cureJMbayes(formFixed = CD4 ~ obstime,
#'                        formRandom = ~ obstime,
#'                        timeVar= "obstime",
#'                        formLatency = Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                        formIncidence = ~ drug + gender + prevOI + AZT,
#'                        formID= ~ id,
#'                        IdVar = "id",
#'                        data = aids,
#'                        # model specifications
#'                        jointCureModel = c("FJCmodel"),
#'                        survMod = c("weibull-PH"),
#'                        param = c("shared-RE"),
#'                        # prior options
#'                        Infprior_cure = TRUE,
#'                        smcure_out = smcure_out,
#'                        priorTau = 100,
#'                        # classification options
#'                        classif_trick = TRUE,
#'                        cov_prior = c("inverse-gamma"),
#'                        Sigma_d = TRUE,
#'                        # MCMC options
#'                        n.chains = 1, n.iter = 10000, n.burnin = 5000, n.thin = 1, n.adapt = 5000, quiet = FALSE, C = 1,
#'                        # out option
#'                        out_data=T)
#'
#' ## details of the estimated model
#' summary(FJCM_A2)
jointCureModel <- function(formFixed,
                           formRandom,
                           formLatency,
                           formIncidence,
                           formID,
                           IdVar,
                           data,
                           timeVar,
                           jointCureModel = c("FJCmodel"),
                           survMod = c("weibull-PH"),
                           param = c("shared-RE"),
                           Infprior_cure=F,
                           smcure_out=NULL,
                           classif_trick=TRUE,
                           cov_prior=c("inverse-gamma"),
                           Sigma_d=FALSE,
                           n.chains = 1,
                           n.iter = 10000,
                           n.burnin = 5000,
                           n.thin = 1,
                           n.adapt = 5000,
                           C = 1,
                           priorTau = 100,
                           out_data=T)
  {


  #----- versions
  # May-2017    : first stable version;
  # Aug-2017    : addition of smcure object in arguments used in simulation programs;
  # 7-sept-2017 : addition of "shared-RE" as association structure;
  # 14-sept-2017: solving one random effect issue and last stable version (current version);
  # _-Novem-2017: vague priors without data-driven parameter priors
  #               + independence between random effects (inverse-gamma prior distributions while inverse-wishart)
  #               + option about D==0 if T_obs>max(t_obs[which(delta==1)])
  #               + calibration of markov chain initial values
  # May-2018    : Add dynamic prediction of the class membership using the draw of variable Ds
  #               + extension to ordinal longitudinal responses
  #               + add the specific-class variance of the error measurements


  #--------------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------------
  # -------- WARNINGs !!!
  # computed and put a dic in output (must be studied)
  # Add the specific-class variance of the error measurements
  # Add class membership prediction
  # -------- Improvements !!!
  # Add other latent strutures (Slope, AUC)
  # Add other survival models (B-spline)
  # Initial values: put random effect predictions
  # Add other GLMM: probability distribution function for longitudinal response variable and for random effets ? cumulative model or student
  #---------------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------------

  lag = 0
  quiet = FALSE

  #--- initial conditions
  if((jointCureModel %in%  c("JSECmodel","FJSECmodel") && Sigma_d)){
    warning("The 'JSECmodel', 'FJSECmodel' do not include class-specific random effect covariance, then 'Sigma_d==FALSE'\n")
    Sigma_d <- FALSE
  }
  if((jointCureModel %in%  c("JSECmodel","FJSECmodel") && cov_prior=="inverse-gamma")){
    warning("Only wishart prior distribution is implemended for 'JSECmodel', 'FJSECmodel'. \n")
    cov_prior <- "wishart"
  }



  #### ---- cure part ----

  ## data
  data_cure <- data[c(all.vars(formID),all.vars(formLatency),all.vars(formIncidence))]
  data_cure <- unique(data_cure)
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
  # time_threshold <-  max(range(boxplot(data_cure$t_obs[which(data_cure$delta==1)],plot=FALSE)$stats))
  time_threshold <- max(Time[which(delta==1)])
  if(classif_trick){
    D[which(delta==0 & Time>time_threshold)] <- 0
  }

  if(Infprior_cure && is.null(smcure_out)){
    ## Call function smcure
    smcure_out <- smcure::smcure(formLatency,
                                 cureform=formIncidence,
                                 data = data_cure,
                                 model="ph")
    # prior parameters for latency model
    priorMean.alpha <- as.numeric(c(0,smcure_out$beta))
    # priorTau.alpha <- diag(rep(1/priorTau,length(priorMean.alpha))) # if informative : 1/(as.numeric(smcure_out$beta_var)*100)
    priorTau.alpha <- diag(1/c(1,smcure_out$beta_var)/100)
    # prior parameters for incidence model
    priorMean.gamma <- as.numeric(smcure_out$b)
    priorTau.gamma <- diag(1/smcure_out$b_var/100)
  }
  if(Infprior_cure && !is.null(smcure_out)){
    # prior parameters for latency model
    priorMean.alpha <- as.numeric(c(0,smcure_out$beta))
    # priorTau.alpha <- diag(rep(1/priorTau,length(priorMean.alpha))) # if informative : 1/(as.numeric(smcure_out$beta_var)*100)
    priorTau.alpha <- diag(1/c(1,smcure_out$beta_var)/100)
    # prior parameters for incidence model
    priorMean.gamma <- as.numeric(smcure_out$b)
    priorTau.gamma <- diag(1/smcure_out$b_var/100)
  }
  if(!Infprior_cure){
    # prior parameters for latency model
    priorMean.alpha <- as.numeric(rep(0,ncol(Z)))
    priorTau.alpha <- diag(rep(1/priorTau,length(priorMean.alpha)))
    # prior parameters for incidence model
    priorMean.gamma <- as.numeric(rep(0,ncol(W)))
    priorTau.gamma <- diag(rep(1/priorTau,length(priorMean.gamma)))
  }

  initial.values <- list(gamma = priorMean.gamma,
                         alpha = priorMean.alpha,
                         shape = 3)


  if (jointCureModel %in% c("MCmodel")) {  # c("FJCmodel","JLCCmodel","JSECmodel","MCmodel")
    jags.data <- list(C = C, zeros = numeric(nTime),
                      Time = Time, delta=delta,
                      N = nTime, ncZ = ncol(Z), ncW = ncol(W),
                      class = D,
                      Z = Z,
                      W = W,
                      priorMean.gamma = priorMean.gamma, priorTau.gamma = priorTau.gamma,
                      priorMean.alpha = priorMean.alpha, priorTau.alpha = priorTau.alpha,
                      priorA.shape = 1/priorTau, priorB.shape = 1/priorTau)
    # initial.values <- list(gamma = priorMean.gamma, alpha= priorMean.alpha,
    #                        shape= 1)

    param <- "shared-RE"
    one.RE <- TRUE

  }else{

    #### ---- longitudinal part ----

    data_long <- data[unique(c(all.vars(formID),all.vars(formFixed),all.vars(formRandom)))]
    y.long <- data_long[all.vars(formFixed)][,1]
    mfX <- model.frame(formFixed, data = data_long)
    X <- model.matrix(formFixed, mfX)
    mfU <- model.frame(formRandom, data = data_long)
    U <- model.matrix(formRandom, mfU)
    # one.RE <- ncol(U) == 1
    # if (one.RE){
    #   U <- cbind(U, rep(0, nrow(U)))
    #   warning("'cureJMbayes' is more efficient when at least two random effects are considered.\n")
    # }
    one.RE <- formRandom=="~1"
    if (one.RE && cov_prior %in% c("wishart")){
      cov_prior <- "inverse-gamma"
      warning("'td-value' with only random intercept is only considered for cov_prior='inverse-gamma'.\n")
    }
    # if(Sigma_d && cov_prior %in% c("wishart")){
    #   Sigma_d <- FALSE
    #   warning(" class-specific covariance matrixis not implemented for cov_prior='wishart'.\n")
    # }
    id <- as.integer(data[all.vars(formID)][,1])
    offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

    if (length(unique(id)) != length(Time))
      stop("sample sizes in the longitudinal and event processes differ.\n")
    N <- length(Time)

    # prior parameters for longitudinal model
    if (jointCureModel %in% c("FJCmodel", "JLCCmodel")) {
      # who is cure ?
      mfW2 <- model.frame(formIncidence, data = data)
      W2 <- model.matrix(formIncidence, mfW2)  # or: W <- cbind(1,data_cure[all.vars(formIncidence)])
      D_hat <- round(data[all.vars(formLatency)][, 2]
                     + (1-data[all.vars(formLatency)][, 2])* 1/(1+exp(-as.matrix(W2)%*%as.vector(smcure_out$b))))
      # if(classif_trick){
        D_hat[which(data[all.vars(formLatency)][, 2]==0 & data[all.vars(formLatency)][, 1]>time_threshold)] <- 0
      # }
      data <- cbind(data,D_hat)
      # Estimation of priorMean by class given class 1 for D==1 and class 2 for D==0
      long_model0 <- lcmm::hlme(fixed = formFixed , random= formRandom ,subject = IdVar, data=data[which(data$D_hat==0),])
      priorMean.beta2 <- long_model0$best[1:ncol(X)]
      priorTau.beta2 <- diag(rep(1/priorTau,length(priorMean.beta2)))
      long_model1 <- lcmm::hlme(fixed = formFixed , random= formRandom ,subject = IdVar, data=data[which(data$D_hat==1),])
      priorMean.beta1 <- long_model1$best[1:ncol(X)]
      priorTau.beta1 <- diag(rep(1/priorTau,length(priorMean.beta1)))
      # if (one.RE)
      #   priorR.Sigma2[ncol(priorR.Sigma2),ncol(priorR.Sigma2)] <- 1
      sigma2 <- long_model1$best["stderr"]
      jags.data <- list(C = C, zeros = numeric(nTime),
                        Time = Time, delta=delta,
                        N = nTime, ncZ = ncol(Z), ncW = ncol(W),
                        Z = Z, W = W,
                        class = D,
                        priorMean.gamma = priorMean.gamma, priorTau.gamma = priorTau.gamma,
                        priorMean.alpha = priorMean.alpha, priorTau.alpha = priorTau.alpha,
                        priorA.shape = 1/priorTau, priorB.shape = 1/priorTau,
                        offset = offset,
                        ncX = ncol(X), # ncU =  ncol(U),
                        y = y.long, X = X, U = U,
                        priorA.tau = (1/sigma2)^2/10, priorB.tau = (1/sigma2)/10,
                        priorMean.beta2 = priorMean.beta2, priorTau.beta2 = priorTau.beta2,
                        priorMean.beta1 = priorMean.beta1, priorTau.beta1 = priorTau.beta1)
      # initRE <- rbind(long_model0$predRE,long_model1$predRE)
      # initRE <- initRE[order(initRE[,1],decreasing=F), ]
      # initial.values <- list(gamma = priorMean.gamma, alpha = priorMean.alpha, shape = 1, alphaL = 0,
      #                        beta = matrix(c(priorMean.beta1, priorMean.beta2), byrow = T, ncol=ncol(X), nrow = 2), xi = initRE,
      #                        prec.tau2 = 1, prec.Sigma2 = diag(1,nrow = ncol(U))
      #                        )
      initial.values$beta <- matrix(c(priorMean.beta1,priorMean.beta2),nrow=2, byrow = T)
      initial.values$prec.tau <- 1/sigma2
      ncU <- ncol(U)
      m1 <- m0 <- matrix(0, ncol = ncU, nrow = ncU)
      m1[lower.tri(m1, TRUE)] <- long_model1$best[ grepl( "varcov" , names(long_model1$best ), fixed = TRUE ) ]
      m0[lower.tri(m0, TRUE)] <- long_model0$best[ grepl( "varcov" , names(long_model1$best ), fixed = TRUE ) ]
      m1 <- m1 + t(m1)
      m0 <- m0 + t(m0)
      diag(m1) <- diag(m1)/2
      diag(m0) <- diag(m0)/2
      m0 <- solve(m0)
      m1 <- solve(m1)
      if(!Sigma_d){
        if(cov_prior=="inverse-gamma")
          initial.values$prec.Sigma2 <- (diag(m0)+diag(m1))/2
        else initial.values$prec.Sigma2 <- (m0+m1)/2
      }else{
        if(cov_prior=="inverse-gamma")
          initial.values$prec.Sigma2 <- as.matrix(rbind(diag(m1),diag(m0)))
        else{
          m <- array(0,dim = c(2,ncU,ncU))
          m[1,,] <- m1
          m[2,,] <- m0
          initial.values$prec.Sigma2 <- m
        }
      }
      initial.values$xi <- matrix(0, ncol=ncol(U), nrow = N)
    }
    # if (jointCureModel %in% c("JSECmodel")) {
    #   # utiliser lme4 sur tout le jeu de donn?e
    #   long_model <- lcmm::hlme(fixed = formFixed , random= formRandom ,subject = IdVar, data=data)
    #   priorMean.beta <- long_model$best[1:ncol(X)]
    #   priorTau.beta <- diag(rep(1/priorTau,length(priorMean.beta)))
    #   priorR.Sigma2 <- diag(rep(1/priorTau,ncol(U)))
    #   # if (one.RE)
    #   #   priorR.Sigma2[ncol(priorR.Sigma2),ncol(priorR.Sigma2)] <- 1
    #   sigma2 <- long_model$best["stderr"]
    #   jags.data <- list(C = C, zeros = numeric(nTime),
    #                     Time = Time, delta=delta,
    #                     N = nTime, ncZ = ncol(Z), ncW = ncol(W),
    #                     y = y.long, Z = Z, W = W,
    #                     priorMean.gamma = priorMean.gamma, priorTau.gamma = priorTau.gamma,
    #                     priorMean.alpha = priorMean.alpha, priorTau.alpha = priorTau.alpha,
    #                     priorA.shape = 1/priorTau, priorB.shape = 1/priorTau,
    #                     offset = offset,
    #                     ncX = ncol(X), # ncU =  ncol(U),
    #                     X = X, U = U,
    #                     priorA.tau = (1/sigma2)^2/10, priorB.tau = (1/sigma2)/10,
    #                     priorMean.beta = priorMean.beta, priorTau.beta = priorTau.beta)
    #   # initRE <- long_model$predRE
    #   # initial.values <- list(gamma = priorMean.gamma, alpha = priorMean.alpha, shape = 1, alphaL = 0,
    #   #                        beta = priorMean.beta, xi = initRE,
    #   #                        prec.tau2 = 1, prec.Sigma2 = diag(1,nrow = ncol(U))
    #   # )
    #   initial.values$beta <- priorMean.beta
    #   initial.values$prec.tau <- 1/sigma2
    #   ncU <- ncol(U)
    #   m <- matrix(0, ncol = ncU, nrow = ncU)
    #   m[lower.tri(m, TRUE)] <- long_model$best[ grepl( "varcov" , names(long_model$best ), fixed = TRUE ) ]
    #   m <- m + t(m)
    #   diag(m) <- diag(m)/2
    #   m <- solve(m)
    #   initial.values$prec.Sigma2 <- m
    #   initial.values$xi <- as.matrix(long_model$predRE[,-1])
    # }
    # if (jointCureModel %in% c("FJSECmodel")) {
    #   long_model <- lcmm::hlme(fixed = formFixed , random= formRandom ,subject = IdVar, data=data)
    #   priorMean.beta <- long_model$best[1:ncol(X)]
    #   priorTau.beta <- diag(rep(1/priorTau,length(priorMean.beta)))
    #   priorR.Sigma2 <- diag(rep(1/priorTau,ncol(U)))
    #   # if (one.RE)
    #   #   priorR.Sigma2[ncol(priorR.Sigma2),ncol(priorR.Sigma2)] <- 1
    #   sigma2 <- long_model$best["stderr"]
    #   initRE <- long_model$predRE
    #   # prior parameters for incidence model
    #   if(one.RE){
    #     priorMean.gammaRE <- 0
    #     priorTau.gammaRE <- 1/priorTau
    #   }else{
    #     priorMean.gammaRE <- as.numeric(rep(0,ncol(U)))
    #     priorTau.gammaRE <- diag(rep(1/priorTau,length(priorMean.gammaRE)))
    #   }
    #   jags.data <- list(C = C, zeros = numeric(nTime),
    #                     Time = Time, delta=delta,
    #                     N = nTime, ncZ = ncol(Z), ncW = ncol(W),
    #                     y = y.long, Z = Z, W = W,
    #                     priorMean.gamma = priorMean.gamma, priorTau.gamma = priorTau.gamma,
    #                     priorMean.alpha = priorMean.alpha, priorTau.alpha = priorTau.alpha,
    #                     priorA.shape = 1/priorTau, priorB.shape = 1/priorTau,
    #                     offset = offset,
    #                     ncX = ncol(X), # ncU =  ncol(U),
    #                     X = X, U = U,
    #                     priorA.tau = (1/sigma2)^2/10, priorB.tau = (1/sigma2)/10,
    #                     priorMean.beta = priorMean.beta, priorTau.beta = priorTau.beta,
    #                     priorMean.gammaRE = priorMean.gammaRE, priorTau.gammaRE = priorTau.gammaRE)
    #   initial.values <- list(gamma = c(priorMean.gamma, priorMean.gammaRE),
    #                          alpha = priorMean.alpha,
    #                          shape = 3)
    #   initial.values$beta <- priorMean.beta
    #   initial.values$prec.tau <- 1/sigma2
    #   ncU <- ncol(U)
    #   m <- matrix(0, ncol = ncU, nrow = ncU)
    #   m[lower.tri(m, TRUE)] <- long_model$best[ grepl( "varcov" , names(long_model$best ), fixed = TRUE ) ]
    #   m <- m + t(m)
    #   diag(m) <- diag(m)/2
    #   m <- solve(m)
    #   initial.values$prec.Sigma2 <- m
    #   initial.values$xi <- as.matrix(long_model$predRE[,-1])
    # }

    #--- td-value case
    data.id <- data_long[!duplicated(id), ]
    if (!timeVar %in% names(data_long))
      stop("\n'timeVar' does not correspond to one of the columns in formulas")
    if (jointCureModel %in% c("FJCmodel", "JSECmodel","FJSECmodel") && param %in% c("td-value")) {
      data.id[[timeVar]] <- pmax(Time - lag, 0)
      mfX.id <- model.frame(formFixed, data = data.id)
      Xtime <- model.matrix(formFixed, mfX.id)
      mfU.id <- model.frame(formRandom, data = data.id)
      Utime <- model.matrix(formRandom, mfU.id)
      # if (one.RE)
      #   Utime <- cbind(Utime, rep(0, nrow(Utime)))
      jags.data <- c(jags.data, list(Xtime = Xtime, Utime = Utime))
    }


    #---- approxitmation of the intergral via the Gaussian quadrature (Gauss Kronrod rule)

    if (jointCureModel %in% c("FJCmodel", "JSECmodel","FJSECmodel") && param %in% c("td-value")) {

      gaussKronrod <-
        function (k = 15) {
          sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
                  0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                  -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                  0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
          wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                    0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                    0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                    0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
          wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
                   0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
          if (k == 7)
            list(sk = sk[1:7], wk = wk7)
          else
            list(sk = sk, wk = wk15)
        }

      wk <- gaussKronrod()$wk
      sk <- gaussKronrod()$sk
      K <- length(sk)
      P <- Time/2
      st <- outer(P, sk + 1)
      id.GK <- rep(seq_along(Time), each = K)
      data.id2 <- data.id[id.GK, ]
      data.id2[[timeVar]] <- c(t(st))
      mfX <- model.frame(formFixed, data = data.id2)
      mfU <- model.frame(formRandom, data = data.id2)
      Xs <- model.matrix(formFixed, mfX)
      Us <- model.matrix(formRandom, mfU)
      # if (one.RE)
      #   Us <- cbind(Us, rep(0, nrow(Us)))

      jags.data <- c(jags.data, list(K = K, P = P, st = st, wk = wk, Xs = Xs, Us = Us))
    }

    #--- management of prior distriution parameter associated with random effects
    # covariance matrix
    if(cov_prior=="inverse-gamma"){
      jags.data <- c(jags.data, list(priorA.Sigma2 = 1/priorTau,
                                     priorB.Sigma2 = 1/priorTau,
                                     ncU = ncol(U)))
    }else{ #wishart distribution
      jags.data <- c(jags.data, list(mu0 = rep(0, ncol(U)),
                                     priorR.Sigma2 = diag(rep(1/priorTau, ncol(U))),
                                     priorK.Sigma2 = ncol(U),
                                     ncU =  ncol(U)))
    }
    # association parameters
    if(jointCureModel %in% c("FJCmodel", "JSECmodel","FJSECmodel")){
          jags.data <- c(jags.data, list(priorTau.alphaA = 1/priorTau))
    }

  }

  #--- Model specification
  model <- switch(paste(jointCureModel, survMod, param, cov_prior, Sigma_d, sep = "/"),
                  ### Weibull, td-value wishart
                  # `MCmodel/weibull-PH/td-value/wishart/FALSE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/td-value/wishart/FALSE` = JLCCM.Weib,
                  `FJCmodel/weibull-PH/td-value/wishart/FALSE` = FJCM.Weib.td.value,
                  `JSECmodel/weibull-PH/td-value/wishart/FALSE` = JSECM.Weib.td.value,
                  `FJSECmodel/weibull-PH/td-value/wishart/FALSE` = FJSECM.Weib.td.value,
                  ### Weibull, shared-RE wishart
                  `MCmodel/weibull-PH/shared-RE/wishart/FALSE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/shared-RE/wishart/FALSE` = JLCCM.Weib,
                  `FJCmodel/weibull-PH/shared-RE/wishart/FALSE` = FJCM.Weib.sharedRE,
                  `JSECmodel/weibull-PH/shared-RE/wishart/FALSE` = JSECM.Weib.sharedRE,
                  `FJSECmodel/weibull-PH/shared-RE/wishart/FALSE` = FJSECM.Weib.sharedRE,
                  # ### Weibull, td-value wishart class-specific
                  `MCmodel/weibull-PH/td-value/wishart/TRUE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/td-value/wishart/TRUE` = JLCCM.Weib.d,
                  `FJCmodel/weibull-PH/td-value/wishart/TRUE` = FJCM.Weib.td.value.d,
                  # ### Weibull, shared-RE wishart class-specific
                  `MCmodel/weibull-PH/shared-RE/wishart/TRUE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/shared-RE/wishart/TRUE` = JLCCM.Weib.d,
                  `FJCmodel/weibull-PH/shared-RE/wishart/TRUE` = FJCM.Weib.sharedRE.d,
                  ### Weibull, td-value inverse-gamma
                  `MCmodel/weibull-PH/td-value/inverse-gamma/FALSE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/td-value/inverse-gamma/FALSE` = JLCCM.Weib.IG,
                  `FJCmodel/weibull-PH/td-value/inverse-gamma/FALSE` = FJCM.Weib.td.value.IG,
                  ### Weibull, shared-RE inverse-gamma
                  `MCmodel/weibull-PH/shared-RE/inverse-gamma/FALSE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/shared-RE/inverse-gamma/FALSE` = JLCCM.Weib.IG,
                  `FJCmodel/weibull-PH/shared-RE/inverse-gamma/FALSE` = FJCM.Weib.sharedRE.IG,
                  ### Weibull, td-value inverse-gamma
                  `MCmodel/weibull-PH/td-value/inverse-gamma/TRUE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/td-value/inverse-gamma/TRUE` = JLCCM.Weib.IGd,
                  `FJCmodel/weibull-PH/td-value/inverse-gamma/TRUE` = FJCM.Weib.td.value.IGd,
                  ### Weibull, shared-RE inverse-gamma
                  `MCmodel/weibull-PH/shared-RE/inverse-gamma/TRUE` = MCM.Weib,
                  `JLCCmodel/weibull-PH/shared-RE/inverse-gamma/TRUE` = JLCCM.Weib.IGd,
                  `FJCmodel/weibull-PH/shared-RE/inverse-gamma/TRUE` = FJCM.Weib.sharedRE.IGd
                  )

  #--- Parameter to save specification
  parms_to_save <- c("alpha", "gamma")
  parms_to_save <- switch(paste(jointCureModel, survMod, param, sep = "/"),
                          ### Weibull, td-value
                          `MCmodel/weibull-PH/td-value` = c(parms_to_save, "shape","class"),
                          `JLCCmodel/weibull-PH/td-value` = c(parms_to_save, "shape","beta","prec.Sigma2","prec.tau2","xi"),
                          `FJCmodel/weibull-PH/td-value` = c(parms_to_save, "shape","beta","prec.Sigma2","prec.tau2","alphaL","xi"),
                          `JSECmodel/weibull-PH/td-value` = c(parms_to_save, "shape","beta","prec.Sigma2","prec.tau2","alphaL","xi"),
                          `FJSECmodel/weibull-PH/td-value` = c(parms_to_save, "shape","gammaRE","beta","prec.Sigma2","prec.tau2","alphaL","xi"),
                          ### Weibull, shared-RE
                          `MCmodel/weibull-PH/shared-RE` = c(parms_to_save, "shape","class"),
                          `JLCCmodel/weibull-PH/shared-RE` = c(parms_to_save, "shape","beta","prec.Sigma2","prec.tau2","xi"),
                          `FJCmodel/weibull-PH/shared-RE` = c(parms_to_save, "shape","beta","prec.Sigma2","prec.tau2","alphaRE","xi"),
                          `JSECmodel/weibull-PH/shared-RE` = c(parms_to_save, "shape","beta","prec.Sigma2","prec.tau2","alphaRE","xi"),
                          `FJSECmodel/weibull-PH/shared-RE` = c(parms_to_save, "shape","gammaRE","beta","prec.Sigma2","prec.tau2","alphaRE","xi")
                          )

  working.directory = getwd()
  write.model.jags(model=model,
                   intitled=file.path(working.directory,"JagsModel.txt"),
                   Data=jags.data,
                   jointCureModel=jointCureModel,
                   param=param,
                   one.RE=one.RE)

  if(jointCureModel %in% c("FJCmodel", "JSECmodel","FJSECmodel") && param %in% c("td-value"))
    initial.values$alphaL <- c(0)
  if(jointCureModel %in% c("FJCmodel", "JSECmodel","FJSECmodel") && param %in% c("shared-RE"))
    initial.values$alphaRE <- rep(0,ncol(U))
  initial.values <- initial.values[names(initial.values) %in% parms_to_save]

  #---- Call to JAGS to estimate the model

  if (!require("rjags"))
    stop("'rjags' is required.\n")
  JMjags.model <- rjags::jags.model(file = "JagsModel.txt",
                                    data = jags.data,
                                    inits = list(initial.values),
                                    n.chains = n.chains,
                                    n.adapt = n.adapt,
                                    quiet = quiet)
  rjags::update(JMjags.model, n.burnin)
  fit <- rjags::coda.samples(JMjags.model,
                             variable.names=parms_to_save,
                             n.iter = n.iter - n.burnin,
                             thin = n.thin)

  #--- Management of the outputs

  # file.remove(file.path(con$working.directory, "JagsModel.txt"))
  codaFit <- as.mcmc.list(fit)

  # draws of variables
  if(jointCureModel %in% c("FJCmodel","JLCCmodel","JSECmodel","FJSECmodel")){
    # draws of the class membership variable
    Bs <- do.call(rbind, codaFit)
    ind.ds <- grep("class[", colnames(Bs), fixed = TRUE)
    var_class <- Bs[, ind.ds, drop = FALSE]
    # draws of random effects
    ind.bs <- grep("xi[", colnames(Bs), fixed = TRUE)
    ranef <- Bs[, ind.bs, drop = FALSE]
    if(!one.RE){
      ord.col <- sapply(strsplit(colnames(ranef), "[", fixed = TRUE),"[", 2)
      ord.col <- sapply(strsplit(ord.col, ",", fixed = TRUE), "[",1)
      ord.col <- order(as.numeric(ord.col))
      ranef <- ranef[, ord.col, drop = FALSE]
    }
    postMeans <- matrix(colMeans(ranef), ncol = ncol(U), byrow = TRUE)
    dimnames(postMeans) <- list(levels(factor(id)), colnames(U))
    postVars <- vector("list", nrow(postMeans))
    ind.var <- matrix(seq_len(ncol(ranef)), ncol = ncol(U), byrow = TRUE)
    for (i in seq_along(postVars)){
      postVars[[i]] <- var(ranef[, ind.var[i, ]])
    }
    if(!one.RE){
      postVars[] <- lapply(postVars,
                           function(m) {
                             dimnames(m) <- list(colnames(postMeans), colnames(postMeans))
                             m
                           }
      )
    }
    names(postVars) <- rownames(postMeans)
    # if (one.RE) {
    #   postMeans <- postMeans[, 1, drop = FALSE]
    #   postVars <- lapply(postVars, function(m) m[1, 1, drop = FALSE])
    # }
    n.sims <- nrow(Bs)
    sims.list <- vector("list", length(parms_to_save))
    names(sims.list) <- parms_to_save
    for (p in seq_along(parms_to_save)) {
      ii <- grep(paste("^", parms_to_save[p], sep = ""), colnames(Bs))
      sims.list[[p]] <- Bs[, ii]
    }
  # h <- function(b, theta) {
  #   theta <- relist(theta, skeleton = list.theta)
  #   beta <- theta$beta
  #   Sigma <- theta$Sigma
  #   tau <- theta$tau
  #   gamma <- theta$gamma
  #   alpha <- theta$alpha
  #   shape <- theta$shape
  #   n <- length(y$Time)
  #   # posterior class membership following gamma parameters, delta and survival function...
  #   # consideration of index output ?
  #   mu.y <- c(x$X %*% beta) + rowSums(x$U * xi[id, , drop = FALSE])
  #   log.p.y.xi <- dnorm(y$y, mu.y, tau, log = TRUE)
  #   log.p.y.xi <- tapply(log.p.y.xi, id, sum)
  #   log.p.xi <- dmvnorm(xi, rep(0, ncol(xi)), Sigma, log = TRUE)
  #   eta.t <- c(x$W %*% rep(alpha, length.out = ncol(x$W)))
  #   id.GK <- rep(seq_along(log(y$Time)), each = 15)
  #   wk.long <- rep(wk, n)
  #   if  (jointCureModel %in% c("JSECmodel")) {  # c("FJCmodel","JLCCmodel","JSECmodel","MCmodel")
  #     Y <- c(x$Xtime %*% beta) + rowSums(x$Utime * xi)
  #     Ys <- c(x$Xs %*% beta) + rowSums(x$Us * xi[id.GK, , drop = FALSE])
  #   }
  #   longSurv <- switch(param, `td-value` = alphas * Y, `td-extra` = Dalphas *
  #                        Yderiv, `td-both` = alphas * Y + Dalphas * Yderiv,
  #                      `shared-RE` = c(b %*% alphas))
  #   longSurv.s <- switch(param, `td-value` = alphas * Ys,
  #                        `td-extra` = Dalphas * Ys.deriv, `td-both` = alphas *
  #                          Ys + Dalphas * Ys.deriv, `shared-RE` = c(b %*%
  #                                                                     alphas)[id.GK])
  #   if (survMod == "weibull-PH") {
  #     log.hazard <- log(sigma.t) + (sigma.t - 1) * y$logT +
  #       eta.t + longSurv
  #     log.survival <- -exp(eta.t) * P * tapply(wk.long *
  #                                                exp(log(sigma.t) + (sigma.t - 1) * log(c(t(st))) +
  #                                                      longSurv.s), id.GK, sum)
  #   }
  #   else {
  #     log.hazard <- c(x$W2 %*% Bs.gammas) + eta.t + longSurv
  #     log.survival <- -exp(eta.t) * P * tapply(wk.long *
  #                                                exp(c(W2s %*% Bs.gammas) + longSurv.s), id.GK,
  #                                              sum)
  #   }
  #   log.p.t.b <- event * log.hazard + log.survival
  #   -2 * sum(log.p.y.b + log.p.t.b + log.p.b, na.rm = TRUE)
  # }
    #--- out
    out <- list(codaFit = lapply(codaFit, function(x) x[, -ind.bs, drop = FALSE]),
                postMeans = postMeans,
                postVars = postVars)

    out$sims.list.xi <- sims.list$xi
    sims.list$xi <- NULL

  }else{ # for MCM
    # draws of the class membership variable
    Bs <- do.call(rbind, codaFit)
    ind.ds <- grep("class[", colnames(Bs), fixed = TRUE)
    Ds <- Bs[, ind.ds, drop = FALSE]
    n.sims <- nrow(Bs)
    sims.list <- vector("list", length(parms_to_save))
    names(sims.list) <- parms_to_save
    for (p in seq_along(parms_to_save)) {
      ii <- grep(paste("^", parms_to_save[p], sep = ""), colnames(Bs))
      sims.list[[p]] <- Bs[, ii]
    }
    # output
    out <- list(codaFit = lapply(codaFit, function(x) x[, -ind.ds, drop = FALSE]))

    out$sims.list.Ds <- sims.list$class
    sims.list$class <- NULL
    }

  class(out$codaFit) <- "mcmc.list"

  if(jointCureModel %in% c("MCmodel")){
    sims.list <- list(gamma = sims.list$gamma, alpha = sims.list$alpha, shape = sims.list$shape)
    colnames(sims.list$gamma) <- colnames(W)
    colnames(sims.list$alpha) <- c(colnames(Z))
    out$sims.list <- sims.list
    out$coefficients <- lapply(sims.list, function(x) if (is.matrix(x))
      colMeans(x)
      else mean(x))
    out$StErr <- lapply(sims.list, function(x) {
      f <- function(x) {
        acf.x <- drop(acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
        acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
        ess <- length(x)/(1 + 2 * sum(acf.x))
        sqrt(var(x)/ess)
      }
      if (is.matrix(x))
        apply(x, 2, f)
      else f(x)
    })
    out$StDev <- lapply(sims.list, function(x) if (is.matrix(x))
      apply(x, 2, sd)
      else sd(x))
    out$CIs <- lapply(sims.list, function(x) if (is.matrix(x))
      apply(x, 2, quantile, probs = c(0.025, 0.975))
      else quantile(x, probs = c(0.025, 0.975)))
    out$vcov <- var(do.call(cbind, sims.list))

    out$sims.list <- sims.list
    rm(Ds, Bs, codaFit, sims.list)

  }else{ # for joint cure models
    if(!one.RE){
      ncU <- ncol(U)
      indSigma2 <- cbind(rep(1:ncU, each = ncU), rep(1:ncU, ncU))
      if( !Sigma_d && cov_prior=="wishart"){
        Sigma2 <- t(sapply(seq_len(n.sims),
                           function(i) {
                             m <- matrix(0, ncU, ncU)
                             m[indSigma2] <- sims.list$prec.Sigma2[i, ]
                             d <- solve(m)
                             d[lower.tri(d, TRUE)]
                           }))
        sims.list <- list(beta = sims.list$beta,
                          tau = sqrt(1/sims.list$prec.tau2),
                          gamma = sims.list$gamma,
                          alpha = sims.list$alpha,
                          shape = sims.list$shape,
                          Sigma2 = Sigma2)
      }
      if(Sigma_d && cov_prior=="wishart"){
        sims.list$prec.Sigma2.D1 <- sims.list$prec.Sigma2[ , grepl( "prec.Sigma2[1," , colnames( sims.list$prec.Sigma2 ),fixed = TRUE ) ]
        sims.list$prec.Sigma2.D0 <- sims.list$prec.Sigma2[ , grepl( "prec.Sigma2[2," , colnames( sims.list$prec.Sigma2 ),fixed = TRUE ) ]
        Sigma2.D1 <- t(sapply(seq_len(n.sims),
                           function(i) {
                             m <- matrix(0, ncU, ncU)
                             m[indSigma2] <- sims.list$prec.Sigma2.D1[i, ]
                             d <- solve(m)
                             d[lower.tri(d, TRUE)]
                           }))
        Sigma2.D0 <- t(sapply(seq_len(n.sims),
                              function(i) {
                                m <- matrix(0, ncU, ncU)
                                m[indSigma2] <- sims.list$prec.Sigma2.D0[i, ]
                                d <- solve(m)
                                d[lower.tri(d, TRUE)]
                              }))
        sims.list <- list(beta = sims.list$beta,
                          tau = sqrt(1/sims.list$prec.tau2),
                          gamma = sims.list$gamma,
                          alpha = sims.list$alpha,
                          shape = sims.list$shape,
                          Sigma2.D1 = Sigma2.D1,
                          Sigma2.D0 = Sigma2.D0)
      }
      if(cov_prior=="inverse-gamma" ){
        sims.list <- list(beta = sims.list$beta,
                          tau = sqrt(1/sims.list$prec.tau2),
                          gamma = sims.list$gamma,
                          alpha = sims.list$alpha,
                          shape = sims.list$shape,
                          Sigma2 = 1/sims.list$prec.Sigma2)
      }

    }else{
      sims.list <- list(beta = sims.list$beta,
                        tau = sqrt(1/sims.list$prec.tau2),
                        gamma = sims.list$gamma,
                        alpha = sims.list$alpha,
                        shape = sims.list$shape,
                        sigma2 = 1/sims.list$prec.Sigma2)
    }

    #---- colnames of sims.list's elements
    if (jointCureModel %in% c("FJSECmodel")) {
      colum_RE <- c((ncol(W)+1):length(colnames(sims.list$gamma)))
      # sims.list$gammaRE <- sims.list$gamma[,colum_RE]
      # sims.list$gamma <- sims.list$gamma[,-colum_RE]
      # colnames(sims.list$gammaRE) <- colnames(U)
      colnames(sims.list$gamma)[(1:ncol(W))] <- colnames(W)
    }else{
      colnames(sims.list$gamma) <- colnames(W)
    }
    # colnames(sims.list$alpha)[1:length(colnames(Z))] <- if (jointCureModel %in% c("FJCmodel","JSECmodel","FJSECmodel")) {c(colnames(Z),"Assoc.")}else{c(colnames(Z))}
    colnames(sims.list$alpha)[1:length(colnames(Z))] <- colnames(Z)
    if (jointCureModel %in% c("JSECmodel","FJSECmodel")) {
      colnames(sims.list$beta) <- colnames(X)
    }else{
      sims.list$beta.D1 <- sims.list$beta[ , grepl( "beta[1" , colnames( sims.list$beta ),fixed = TRUE ) ]
      sims.list$beta.D0 <- sims.list$beta[ , grepl( "beta[2" , colnames( sims.list$beta ),fixed = TRUE ) ]
      colnames(sims.list$beta.D1) <- colnames(sims.list$beta.D0) <- colnames(X)
      sims.list$beta <- NULL
    }
    if(!one.RE && !Sigma_d && cov_prior=="wishart"){
      tmp_mat <- matrix(1, ncol = ncU, nrow = ncU)
      colnames(sims.list$Sigma2) <- paste("Sigma2[", row(tmp_mat)[lower.tri(tmp_mat,TRUE)], ", ", col(tmp_mat)[lower.tri(tmp_mat, TRUE)], "]", sep = "")
    }
    if(!one.RE && Sigma_d && cov_prior=="wishart"){
      tmp_mat <- matrix(1, ncol = ncU, nrow = ncU)
      # colnames(sims.list$Sigma2.D1) <- paste("Sigma2.D1[", row(tmp_mat)[lower.tri(tmp_mat,TRUE)], ", ", col(tmp_mat)[lower.tri(tmp_mat, TRUE)], "]", sep = "")
      # colnames(sims.list$Sigma2.D0) <- paste("Sigma2.D0[", row(tmp_mat)[lower.tri(tmp_mat,TRUE)], ", ", col(tmp_mat)[lower.tri(tmp_mat, TRUE)], "]", sep = "")
      colnames(sims.list$Sigma2.D1) <- colnames(sims.list$Sigma2.D0) <- paste("Sigma2[", row(tmp_mat)[lower.tri(tmp_mat,TRUE)], ", ", col(tmp_mat)[lower.tri(tmp_mat, TRUE)], "]", sep = "")
    }
    if(!one.RE && !Sigma_d && cov_prior=="inverse-gamma"){
      colnames(sims.list$Sigma2) <- gsub("prec.", "", colnames(sims.list$Sigma2))
    }
    if(!one.RE && Sigma_d && cov_prior=="inverse-gamma"){
      colnames(sims.list$Sigma2) <- gsub("prec.", "", colnames(sims.list$Sigma2))
      sims.list$Sigma2.D1 <- sims.list$Sigma2[ , grepl( "Sigma2[1," , colnames( sims.list$Sigma2 ),fixed = TRUE ) ]
      sims.list$Sigma2.D0 <- sims.list$Sigma2[ , grepl( "Sigma2[2," , colnames( sims.list$Sigma2 ),fixed = TRUE ) ]
      colnames(sims.list$Sigma2.D1) <- gsub("1,", "", colnames(sims.list$Sigma2.D1))
      colnames(sims.list$Sigma2.D0) <- gsub("2,", "", colnames(sims.list$Sigma2.D0))
      sims.list$Sigma2 <- NULL
    }
    sims.list <- sims.list[!sapply(sims.list, is.null)]
    out$modes <- lapply(sims.list, function(x) {
      m <- function(x) {
        d <- density(x, bw = "nrd", adjust = 3, n = 1000)
        d$x[which.max(d$y)]
      }
      if (is.matrix(x))
        apply(x, 2, m)
      else m(x)
    })
    out$coefficients <- lapply(sims.list, function(x) if (is.matrix(x))
      colMeans(x)
      else mean(x))
    out$StErr <- lapply(sims.list, function(x) {
      f <- function(x) {
        acf.x <- drop(acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
        acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
        ess <- length(x)/(1 + 2 * sum(acf.x))
        sqrt(var(x)/ess)
      }
      if (is.matrix(x))
        apply(x, 2, f)
      else f(x)
    })
    out$StDev <- lapply(sims.list, function(x) if (is.matrix(x))
      apply(x, 2, sd)
      else sd(x))
    out$CIs <- lapply(sims.list, function(x) if (is.matrix(x))
      apply(x, 2, quantile, probs = c(0.025, 0.975))
      else quantile(x, probs = c(0.025, 0.975)))
    out$vcov <- var(do.call(cbind, sims.list))
    if(!one.RE){
      if(cov_prior=="wishart" && !Sigma_d){
        Sigma2 <- matrix(0, ncol = ncU, nrow = ncU)
        Sigma2[lower.tri(Sigma2, TRUE)] <- out$coefficients$Sigma2
        Sigma2 <- Sigma2 + t(Sigma2)
        diag(Sigma2) <- diag(Sigma2)/2
        out$coefficients$Sigma2 <- Sigma2
        dimnames(out$coefficients$Sigma2) <- list(colnames(U),colnames(U))
        Sigma2 <- matrix(0, ncol = ncU, nrow = ncU)
        Sigma2[lower.tri(Sigma2, TRUE)] <- out$modes$Sigma2
        Sigma2 <- Sigma2 + t(Sigma2)
        diag(Sigma2) <- diag(Sigma2)/2
        out$modes$Sigma2 <- Sigma2
        dimnames(out$modes$Sigma2) <- list(colnames(U),colnames(U))
      }
      if(cov_prior=="wishart" && Sigma_d){
        # D=0
        Sigma2 <- matrix(0, ncol = ncU, nrow = ncU)
        Sigma2[lower.tri(Sigma2, TRUE)] <- out$coefficients$Sigma2.D0
        Sigma2 <- Sigma2 + t(Sigma2)
        diag(Sigma2) <- diag(Sigma2)/2
        out$coefficients$Sigma2.D0 <- Sigma2
        dimnames(out$coefficients$Sigma2.D0) <- list(colnames(U),colnames(U))
        Sigma2 <- matrix(0, ncol = ncU, nrow = ncU)
        Sigma2[lower.tri(Sigma2, TRUE)] <- out$modes$Sigma2.D0
        Sigma2 <- Sigma2 + t(Sigma2)
        diag(Sigma2) <- diag(Sigma2)/2
        out$modes$Sigma2.D0 <- Sigma2
        dimnames(out$modes$Sigma2.D0) <- list(colnames(U),colnames(U))
        # D=1
        Sigma2 <- matrix(0, ncol = ncU, nrow = ncU)
        Sigma2[lower.tri(Sigma2, TRUE)] <- out$coefficients$Sigma2.D1
        Sigma2 <- Sigma2 + t(Sigma2)
        diag(Sigma2) <- diag(Sigma2)/2
        out$coefficients$Sigma2.D1 <- Sigma2
        dimnames(out$coefficients$Sigma2.D1) <- list(colnames(U),colnames(U))
        Sigma2 <- matrix(0, ncol = ncU, nrow = ncU)
        Sigma2[lower.tri(Sigma2, TRUE)] <- out$modes$Sigma2.D1
        Sigma2 <- Sigma2 + t(Sigma2)
        diag(Sigma2) <- diag(Sigma2)/2
        out$modes$Sigma2.D1 <- Sigma2
        dimnames(out$modes$Sigma2.D1) <- list(colnames(U),colnames(U))
      }
      if(cov_prior=="inverse-gamma" && !Sigma_d){
        out$coefficients$Sigma2 <- diag(out$coefficients$Sigma2)
        out$modes$Sigma2 <- diag(out$modes$Sigma2)
        dimnames(out$modes$Sigma2) <- dimnames(out$coefficients$Sigma2) <- list(colnames(U),colnames(U))
      }
      if(cov_prior=="inverse-gamma" && Sigma_d){
        out$coefficients$Sigma2.D1 <- diag(out$coefficients$Sigma2.D1)
        out$coefficients$Sigma2.D0 <- diag(out$coefficients$Sigma2.D0)
        out$modes$Sigma2.D1 <- diag(out$modes$Sigma2.D1)
        out$modes$Sigma2.D0 <- diag(out$modes$Sigma2.D0)
        dimnames(out$modes$Sigma2.D1) <- dimnames(out$coefficients$Sigma2.D1) <- dimnames(out$modes$Sigma2.D0) <- dimnames(out$coefficients$Sigma2.D0) <- list(colnames(U),colnames(U))
      }
    }

    # if (one.RE) {
    #   out$coefficients$Sigma2 <- out$coefficients$Sigma2[1, 1, drop = FALSE]
    #   dimnames(out$coefficients$Sigma2) <- list(colnames(U)[1],
    #                                        colnames(U)[1])
    #   if (param == "shared-RE") {
    #     out$coefficients$alpha <- out$coefficients$alpha[-length(out$coefficients$alpha)]
    #     #names(out$coefficients$alpha[length(out$coefficients$alpha)]) <- colnames(U)[1]
    #   }
    #   x$U <- x$U[, 1, drop = FALSE]
    #   if (!is.null(x$Utime))
    #     x$Utime <- x$Utime[, 1, drop = FALSE]
    # }

    #
    # # list.theta <- out$coefficients
    # # theta <- unlist(as.relistable(list.theta))
    # # --- Dic computation here
    # # Sigma2.hat <- h(postMeans, theta)
    # # M <- nrow(sims.list$beta)
    # # hat.Ds <- numeric(M)
    # # for (m in seq_len(M)) {
    # #   postMeans.m <- matrix(ranef[m, ], ncol = ncol(Z), byrow = TRUE)
    # #   thetas.m <- lapply(sims.list, function(x) if (is.matrix(x))
    # #     x[m, ]
    # #     else x[m])
    # #   D <- matrix(0, nrow(VC), ncol(VC))
    # #   D[lower.tri(D, TRUE)] <- thetas.m$D
    # #   D <- D + t(D)
    # #   diag(D) <- diag(D)/2
    # #   thetas.m$D <- D
    # #   thetas.m <- unlist(as.relistable(thetas.m))
    # #   hat.Ds[m] <- h(postMeans.m, thetas.m)
    # # }
    # # hat.D <- mean(hat.Ds, na.rm = TRUE)
    # # out$pD <- hat.D - D.hat
    # # out$DIC <- out$pD + hat.D
    #
    out$sims.list <- sims.list
    rm(ranef, Bs, codaFit, sims.list)
  }

  # --- Out
  if(out_data){
    y <- list(Time = Time, event = delta, zeros = zeros, lag = lag)
    y <- switch(jointCureModel,
                `FJCmodel` = c(y, list(y = y.long, offset = offset)),
                `JLCCmodel` = c(y, list(y = y.long, offset = offset)),
                `JSECmodel` = c(y, list(y = y.long, offset = offset)),
                `FJSECmodel` = c(y, list(y = y.long, offset = offset)),
                `MCmodel` = y)
    x <- list(Z = Z, W = W)
    if(jointCureModel %in% c("FJCmodel","JLCCmodel","JSECmodel","FJSECmodel")){
      x <- switch(paste(jointCureModel, survMod, param, sep = "/"),
                  ### Weibull, td-value
                  `FJCmodel/weibull-PH/td-value` = c(x, list(X = X, U = U, Xtime = Xtime, Utime = Utime)),
                  `JLCCmodel/weibull-PH/td-value` = c(x, list(X = X, U = U)),
                  `JSECmodel/weibull-PH/td-value` = c(x, list(X = X, U = U, Xtime = Xtime, Utime = Utime)),
                  `FJSECmodel/weibull-PH/td-value` = c(x, list(X = X, U = U, Xtime = Xtime, Utime = Utime)),
                  # `MCmodel/weibull-PH/td-value` = x,
                  ### Weibull, shared-RE
                  # `MCmodel/weibull-PH/shared-RE` = x,
                  `JLCCmodel/weibull-PH/shared-RE` = c(x, list(X = X, U = U)),
                  `FJCmodel/weibull-PH/shared-RE` = c(x, list(X = X, U = U)),
                  `JSECmodel/weibull-PH/shared-RE` = c(x, list(X = X, U = U)),
                  `FJSECmodel/weibull-PH/shared-RE` = c(x, list(X = X, U = U))
      )
    }

    out$x <- x
    out$y <- y
    if(jointCureModel %in% c("FJCmodel","JLCCmodel","JSECmodel","FJSECmodel")){
      out$id <- factor(id)
      names(out$id) <- all.vars(formID)
      out$data.id <- data.id
    }
    out$data_cure
    out$times <- data[[timeVar]]
    out$data <- data
    out$smcure <- c(smcure_out$b,smcure_out$beta)
  }

  out$survMod <- survMod
  out$formID <- formID
  out$formRandom <- formRandom
  out$formFixed <- formFixed
  out$formIncidence <- formIncidence
  out$formLatency <- formLatency
  out$survMod <- survMod
  out$timeVar <- timeVar
  out$control <- list(program = "JAGS", n.chains = n.chains, n.iter = n.iter,
                      n.burnin = n.burnin, n.thin = n.thin, n.adapt = n.adapt,
                      C = C, working.directory = getwd(), quiet = FALSE, lag = lag, n = length(Time))
  out$jointCureModel <- jointCureModel
  out$param <- param
  out$cov_prior <- cov_prior
  out$Sigma_d <- Sigma_d
  out$classif_trick <- classif_trick
  class(out) <- "JMcuR"
  out

}#end function
