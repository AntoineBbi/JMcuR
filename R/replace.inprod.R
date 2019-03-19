replace.inprod <- 
  function (body.model, Data, jointCureModel, param, one.RE) {
    mt <- deparse(body.model, width.cutoff = 200L)
    # long
    if(jointCureModel %in% c("FJCmodel", "JLCCmodel")){
      ncX <- Data$ncX
      rplc <- paste(paste("beta[index[i], ",1:ncX, "] * X[j, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[index[i], 1:ncX], X[j, 1:ncX])", rplc, mt, fixed = TRUE)
    }
    if(jointCureModel %in% c("JSLSCmodel","FJSECmodel")){
      ncX <- Data$ncX
      rplc <- paste(paste("beta[",1:ncX, "] * X[j, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[1:ncX], X[j, 1:ncX])", rplc, mt, fixed = TRUE)
    }
    if(jointCureModel %in% c("FJCmodel", "JLCCmodel", "JSLSCmodel","FJSECmodel") && (!one.RE)){
      ncU <- Data$ncU
      rplc <- paste(paste("xi[i, ", 1:ncU, "] * U[j, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(xi[i, 1:ncU], U[j, 1:ncU])", rplc, mt, fixed = TRUE)
    }
    # Incidence
    ncW <- Data$ncW
    rplc <- paste(paste("gamma[", 1:ncW, "] * W[i, ", 1:ncW, "]", sep = ""), collapse = " + ")
    mt <- gsub("inprod(gamma[1:ncW], W[i, 1:ncW])", rplc, mt, fixed = TRUE)
    if(jointCureModel %in% c("FJSECmodel") && (!one.RE)){
      ncU <- Data$ncU
      rplc <- paste(paste("gammaRE[", 1:ncU, "] * xi[i, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(gammaRE[1:ncU], xi[i, 1:ncU])", rplc, mt, fixed = TRUE)
    }
    # Latency
    ncZ <- Data$ncZ
    rplc <- paste(paste("alpha[", 1:ncZ, "] * Z[i, ", 1:ncZ, "]", sep = ""), collapse = " + ")
    mt <- gsub("inprod(alpha[1:ncZ], Z[i, 1:ncZ])", rplc, mt, fixed = TRUE)
    # share latent structure
    if (jointCureModel %in% c("JSLSCmodel","FJSECmodel") && param %in% c("td-value")) {
      rplc <- paste(paste("beta[", 1:ncX, "] * Xtime[i, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[1:ncX], Xtime[i, 1:ncX])", rplc, mt, fixed = TRUE)
      rplc <- paste(paste("xi[i, ", 1:ncU, "] * Utime[i, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(xi[i, 1:ncU], Utime[i, 1:ncU])", rplc, mt, fixed = TRUE)
      #
      rplc <- paste(paste("beta[", 1:ncX, "] * Xs[K * (i - 1) + k, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[1:ncX], Xs[K * (i - 1) + k, 1:ncX])", rplc, mt, fixed = TRUE)
      rplc <- paste(paste("xi[i,", 1:ncU, "] * Us[K * (i - 1) + k, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(xi[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])", rplc, mt, fixed = TRUE)
    }
    if (jointCureModel %in% c("FJCmodel") && param %in% c("td-value")) {
      rplc <- paste(paste("beta[index[i], ", 1:ncX, "] * Xtime[i, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[index[i], 1:ncX], Xtime[i, 1:ncX])", rplc, mt, fixed = TRUE)
      rplc <- paste(paste("xi[i, ", 1:ncU, "] * Utime[i, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(xi[i, 1:ncU], Utime[i, 1:ncU])", rplc, mt, fixed = TRUE)
      #
      rplc <- paste(paste("beta[index[i], ", 1:ncX, "] * Xs[K * (i - 1) + k, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[index[i], 1:ncX], Xs[K * (i - 1) + k, 1:ncX])", rplc, mt, fixed = TRUE)
      rplc <- paste(paste("xi[i, ", 1:ncU, "] * Us[K * (i - 1) + k, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(xi[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])", rplc, mt, fixed = TRUE)
    }
    if (jointCureModel %in% c("FJCmodel","JSLSCmodel","FJSECmodel") && param %in% c("shared-RE") && (!one.RE)) {
      ncU <- Data$ncU
      rplc <- paste(paste("alphaRE[", 1:ncU, "] * xi[i,  ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(alphaRE[1:ncU], xi[i, 1:ncU])", rplc, mt, fixed = TRUE)
    }
    c("model", mt)
  }
