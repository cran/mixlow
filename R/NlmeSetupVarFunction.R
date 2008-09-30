`NlmeSetupVarFunction` <-
function(ord, dat1)
  {
  ## sets up the variance functions that are used when calling NLME
  var.function.m1 = NULL
  var.function.m2 = NULL
  var.function.m3 = NULL
  var.function.m4 = NULL
  var.function.s1 = NULL
  var.function.s2 = NULL
  var.function.s3 = NULL
  var.function.s4 = NULL
  
  if (length(ord) > 1)
    {
    flist = list(0)
    for(ii in 1:length(levels(dat1$drg)))
      {
      strr = paste("flist$'",ii,"' = 0.5", sep="",collapse="")
      eval(parse(text=strr))
      }
    flist = flist[-1]

    var.function.m1 <- NULL
    var.function.m2 <- varIdent(form = ~1|drg)
    var.function.m3 <- varPower(form = ~(fitted(.))|drg, fixed= flist)
    var.function.m4 <- varPower(form = ~(fitted(.))|drg)

    # set fixed effects list
    fixedEffects=list(0)
    fixedEffects[[1]] <- as.formula(g~drg-1)
    fixedEffects[[2]] <- as.formula(p~drg-1)
    fixedEffects[[3]] <- as.formula(u~1)
    }

  # if only one drug
  if (length(ord) == 1)
    {
    var.function.s1 <- NULL
    var.function.s2 <- varPower(form = ~fitted(.), fixed= 0.5)
    var.function.s3 <- varPower(form = ~fitted(.))
    var.function.s4 <- varConstPower(form = ~fitted(.), fixed= list(power = 0.5) )

    # set fixed effects list
    fixedEffects=list(0)
    fixedEffects[[1]] <- as.formula(g~1)
    fixedEffects[[2]] <- as.formula(p~1)
    fixedEffects[[3]] <- as.formula(u~1)
    }


  # set up random effects list
  randomEffects = list(0)
  randomEffects$tray[[1]] = as.formula(u ~ 1)
  randomEffects = randomEffects[-1]
  
  return (list(fixedEffects=fixedEffects, var.function.m1=var.function.m1, var.function.m2=var.function.m2, 
    var.function.m3=var.function.m3, var.function.m4=var.function.m4, var.function.s1=var.function.s1, 
    var.function.s2=var.function.s2, var.function.s3=var.function.s3, var.function.s4=var.function.s4, randomEffects=randomEffects))
  }

