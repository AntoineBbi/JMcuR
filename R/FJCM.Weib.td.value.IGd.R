FJCM.Weib.td.value.IGd <-
  function () {
    # Likelihood using trick of zeros
    for (i in 1:N) {
      # Longitudinal Part
      for(j in offset[i]:(offset[i+1]-1)){
        y[j] ~ dnorm(mu[j], prec.tau2)
        mu[j] <- inprod(beta[index[i], 1:ncX], X[j, 1:ncX]) + inprod(xi[i, 1:ncU], U[j, 1:ncU])
      }
      # Random Effects Part : Hyper-parameterization of random effects
      for(r in 1:ncU){
        xi[i,r] ~ dnorm(0, prec.Sigma2[index[i], r])
      }
      # Incidence part
      logit(pi[i]) <- inprod(gamma[1:ncW], W[i, 1:ncW])
      class[i] ~ dbern(pi[i])
      index[i] <- 1*equals(class[i],1) + equals(class[i],0)*2 
      # index[i] <- delta[i] + (1-delta[i]) * (equals( class[i], 1)*1 + equals(class[i], 0)*2)
      # Latency part
      etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ])
      shareY[i] <- inprod(beta[index[i], 1:ncX], Xtime[i, 1:ncX]) + inprod(xi[i, 1:ncU], Utime[i, 1:ncU])
      log_h1[i] <- log(shape) + (shape - 1) * log(Time[i]) + etaBaseline[i] + alphaL * shareY[i]
      for (k in 1:K) {
        shareY.s[i, k] <- inprod(beta[index[i], 1:ncX], Xs[K * (i - 1) + k, 1:ncX]) + inprod(xi[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])
        SurvLong[i, k] <- wk[k] * shape * pow(st[i, k], shape - 1) * exp(alphaL * shareY.s[i, k])
      }
      log_S1[i] <- (-exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ]))
      logL[i] <- class[i]*log(pi[i]) + class[i]*delta[i]*log_h1[i] + class[i]*log_S1[i] + (1-delta[i])*(1-class[i])*log(1-pi[i])
      mlogL[i] <- -logL[i] + C
      zeros[i] ~ dpois(mlogL[i])
    }
    # Longitudinal priors
    prec.tau2 ~ dgamma(priorA.tau,priorB.tau)
    beta[1,1:ncX] ~ dmnorm(priorMean.beta1[], priorTau.beta1[, ])          # class 1 (uncured subjects, D=1)
    beta[2,1:ncX] ~ dmnorm(priorMean.beta2[], priorTau.beta2[, ])          # class 2 (cured subjects, D=0)
    for(rr in 1:ncU){
      prec.Sigma2[1, rr] ~ dgamma(priorA.Sigma2,priorB.Sigma2)
      prec.Sigma2[2, rr] ~ dgamma(priorA.Sigma2,priorB.Sigma2)
    }
    # Latency priors
    alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
    shape ~ dgamma(priorA.shape,priorB.shape)
    alphaL ~ dnorm(0, priorTau.alphaA)
    # Incidence priors
    gamma[1:ncW] ~ dmnorm(priorMean.gamma[], priorTau.gamma[, ])
  }
