MCM.Weib <-
  function () {
    # Likelihood using trick of zeros
    for (i in 1:N) {
      # Incidence part
      logit(pi[i]) <- inprod(gamma[1:ncW], W[i, 1:ncW])
      class[i] ~ dbern(pi[i])
      # Latency part
      etaBaseline[i] <- inprod(alpha[1:ncZ], Z[i, 1:ncZ])
      log_S1[i] <- -exp(etaBaseline[i]) * pow(Time[i], shape)
      log_h1[i] <- log(shape) + (shape - 1) * log(Time[i]) + etaBaseline[i]
      logL[i] <- class[i] * log(pi[i]) + class[i] * delta[i] * log_h1[i] + class[i] * log_S1[i] + (1 - delta[i]) * (1 - class[i]) * log(1 - pi[i])
      zeros[i] ~ dpois(mlogL[i])
      mlogL[i] <- -logL[i] + C
    }
    # Latency priors
    alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
    shape ~ dgamma(priorA.shape,priorB.shape)
    # Incidence priors
    gamma[1:ncW] ~ dmnorm(priorMean.gamma[], priorTau.gamma[, ])
  }