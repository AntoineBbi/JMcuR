ModelMats <-
function (time, ii) {
    id.GK <- rep(ii, each = 15)
    wk <- gaussKronrod()$wk
    sk <- gaussKronrod()$sk
    P <- time / 2
    st <- P * (sk + 1)
    data.id2 <- data.id[id.GK, ]
    data.id2[[timeVar]] <- pmax(st - lag, 0)
    out <- list(st = st, wk = rep(wk, length(P)), P = P)
    # if (param %in% c("td-value", "td-both")) {
    if (param %in% c("td-value")) {
        out$Xs <- model.matrix(delete.response(terms(object$formFixed)), 
                               model.frame(delete.response(terms(object$formFixed)), 
                                           data = data.id2[all.vars(delete.response(terms(object$formFixed)))]))
        out$Us <- model.matrix(object$formRandom, model.frame(object$formRandom, data = data.id2[all.vars(object$formRandom)]))
        out$Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
    }
    # if (param %in% c("td-extra", "td-both")) {
    #     mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
    #     mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
    #     out$Xs.deriv <- model.matrix(extraForm$fixed, mfX.deriv)
    #     out$Zs.deriv <- model.matrix(extraForm$random, mfZ.deriv)
    #     out$Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
    # }
    out
}
