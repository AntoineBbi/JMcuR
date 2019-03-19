MCM.Weib <-
  function () {
    # Likelihood using trick of zeros
    for (i in 1:N) {
      # Incidence part
      logit(pi[i]) <- inprod(gamma[1:ncW], W[i, 1:ncW])
      # Latency part
      etaBaseline[i] <- inprod(alpha[1:ncZ], Z[i, 1:ncZ])
      logL[i] <- delta[i]*log(pi[i]) + delta[i] * logf_1[i] + (1-delta[i]) * log(1 - pi[i] + pi[i] * S_1[i])
      S_1[i] <- exp(-exp(etaBaseline[i]) * pow(Time[i], shape))
      logf_1[i] <- log(shape) + (shape-1) * log(Time[i]) + etaBaseline[i] - exp(etaBaseline[i]) *pow(Time[i], shape)
      zeros[i] ~ dpois(mlogL[i])
      mlogL[i] <- -logL[i] + C
    }
    # Latency priors
    alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
    shape ~ dgamma(priorA.shape,priorB.shape)
    # Incidence priors
    gamma[1:ncW] ~ dmnorm(priorMean.gamma[], priorTau.gamma[, ])
  }