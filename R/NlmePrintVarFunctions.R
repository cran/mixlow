`NlmePrintVarFunctions` <-
function()   {
  ## Convienience function to print summary of var.function indices

  writeLines ("\nSee nlme() help files for information on specific variance functions.\n")  

  writeLines ("\n  ---- variance functions for single drugs ----\n")
  writeLines ("var.function 1:  sigma")
  writeLines ("var.function 2:  sigma*E[response]")
  writeLines ("var.function 3:  sigma*E[response]^beta")
  writeLines ("var.function 4:  sigma*(beta1 + E[response])")


  writeLines ("\n  ---- variance functions for multiple drugs ----\n")
  writeLines ("var.function 1:  sigma")
  writeLines ("var.function 2:  sigma*alpha, where alpha is drug-dependent")
  writeLines ("var.function 3:  sigma*E[response]")
  writeLines ("var.function 4:  sigma*E[response]^beta, where beta is drug-dependent\n")

  }

