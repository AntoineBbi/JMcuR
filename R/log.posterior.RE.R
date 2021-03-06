log.posterior.RE <-
  function (b, y, Mats, survMod, ii) {
    id.i <- id %in% ii
    idT.i <- idT %in% ii
    X.i <- X[id.i, , drop = FALSE]
    U.i <- U[id.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% beta.new) + rowSums(U.i * rep(b, each = nrow(U.i)))
    logY <- dnorm(y[id.i], mu.y, tau.new, TRUE)
    log.p.yb <- sum(logY)
    log.p.b <- dmvnorm(b, rep(0, ncol(U)), Sigma2.new, TRUE)
    st <- Mats[[ii]]$st
    wk <- Mats[[ii]]$wk
    P <- Mats[[ii]]$P
    Xs <- Mats[[ii]]$Xs
    Us <- Mats[[ii]]$Us
    # Xs.deriv <- Mats[[ii]]$Xs.deriv
    # Zs.deriv <- Mats[[ii]]$Zs.deriv
    Ws.intF.vl <- Mats[[ii]]$Ws.intF.vl
    # Ws.intF.sl <- Mats[[ii]]$Ws.intF.sl
    #ind <- Mats[[ii]]$ind
    if (param %in% c("td-value", "td-both"))
      Ys <- as.vector(Xs %*% beta.new + rowSums(Us * rep(b, each = nrow(Us))))
    # if (param %in% c("td-extra", "td-both"))
    #   Ys.deriv <- as.vector(Xs.deriv %*% betas.new[indFixed]) + 
    #   rowSums(Zs.deriv * rep(b[indRandom], each = nrow(Zs)))
    tt <- switch(param,
                 "td-value" = c(Ws.intF.vl %*% assoc.alpha.new) * Ys,
                 # "td-extra" = c(Ws.intF.sl %*% Dalpha.new) * Ys.deriv,
                 # "td-both" = c(Ws.intF.vl %*% assoc.alpha.new) * Ys + c(Ws.intF.sl %*% Dalpha.new) * Ys.deriv,
                 "shared-RE" = rep(sum(b * assoc.alpha.new), length(st)))
    eta.tz <- if (!is.null(W)) {
      as.vector(Z[ii, , drop = FALSE] %*% alpha.new)
    } else 0
    log.survival <- if (survMod == "weibull-PH") {
      Vi <- exp(log(shape.new) + (shape.new - 1) * log(st) + tt) 
      - exp(eta.tz) * P * sum(wk * Vi)
    } #  else if (survMod == "spline-PH") {
    #   kn <- object$control$knots
    #   W2s <- splineDesign(unlist(kn, use.names = FALSE), st, 
    #                       ord = object$control$ordSpline, outer.ok = TRUE)
    #   Vi <- exp(c(W2s %*% Bs.gammas.new) + tt)
    #   idT <- rep(seq_along(P), each = 15)
    #   - sum(exp(eta.tw) * P * tapply(wk * Vi, idT, sum))
    # }
    if (all(st == 0))
      log.survival <- 1
    log.p.yb + log.survival + log.p.b
  }
