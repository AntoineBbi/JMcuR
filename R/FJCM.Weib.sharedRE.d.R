FJCM.Weib.sharedRE.d <-
  function () {
    # Likelihood using trick of zeros
    for (i in 1:N) {
      # Longitudinal Part
      for(j in offset[i]:(offset[i+1]-1)){
        y[j] ~ dnorm(mu[j], prec.tau2)
        mu[j] <- inprod(beta[index[i], 1:ncX], X[j, 1:ncX]) + inprod(xi[i, 1:ncU], U[j, 1:ncU])
      }
      # Incidence part
      logit(pi[i]) <- inprod(gamma[1:ncW], W[i, 1:ncW])
      class[i] ~ dbern(pi[i])
      index[i] <- 1*equals(class[i],1) + equals(class[i],0)*2
      # Latency part
      etaBaseline[i] <- inprod(alpha[1:ncZ], Z[i, 1:ncZ]) + inprod(alphaRE[1:ncU], xi[i, 1:ncU])
      log_S1[i] <- -exp(etaBaseline[i]) * pow(Time[i], shape)
      log_h1[i] <- log(shape) + (shape-1)*log(Time[i]) + etaBaseline[i]
      logL[i] <- class[i]*log(pi[i]) + class[i]*delta[i]*log_h1[i] + class[i]*log_S1[i] + (1-delta[i])*(1-class[i])*log(1-pi[i])
      # log_f1[i] <- log(shape) + (shape-1)*log(Time[i]) + etaBaseline[i] - exp(etaBaseline[i]) *pow(Time[i], shape)
      # logL[i] <- class[i]*log(pi[i]) + class[i]*delta[i]*log_f1[i] + (1-delta[i])*class[i]*log_S1[i] + (1-delta[i])*(1-class[i])*log(1-pi[i])
      zeros[i] ~ dpois(mlogL[i])
      mlogL[i] <- -logL[i] + C
      # Random Effects Part : Hyper-parameterization of random effects
      xi[i,1:ncU] ~ dmnorm(mu0[], prec.Sigma2[index[i], , ])
    }
    # Longitudinal priors
    prec.tau2 ~ dgamma(priorA.tau,priorB.tau)
    beta[1,1:ncX] ~ dmnorm(priorMean.beta1[], priorTau.beta1[, ])       # class 1 (uncured subjects, D=1)
    beta[2,1:ncX] ~ dmnorm(priorMean.beta2[], priorTau.beta2[, ])       # class 2 (cured subjects, D=0)
    prec.Sigma2[1, 1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[,], priorK.Sigma2)  # Random Effects Part class 1 (uncured subjects, D=1)
    prec.Sigma2[2, 1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[,], priorK.Sigma2)  # Random Effects Part class 2 (cured subjects, D=0)
    # Latency priors
    alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
    # alphaRE[1:ncU] ~ dmnorm(priorMean.alphaRE[], priorTau.alphaRE[, ])
    for(rr in 1:ncU){
      alphaRE[rr] ~ dnorm(0, priorTau.alphaA)
    }
    shape ~ dgamma(priorA.shape,priorB.shape)
    # Incidence priors
    gamma[1:ncW] ~ dmnorm(priorMean.gamma[], priorTau.gamma[, ])
  }