write.model.jags <-
  function (model, intitled, Data, jointCureModel, param, one.RE) {
    model <- replace.inprod(body(model), Data, jointCureModel, param, one.RE)
    writeLines(model, intitled)
  }
