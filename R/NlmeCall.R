`NlmeCall` <-
function(formu, dat1, fixedEffects, starting, randomEffects, varst, method, chk, L1, analysis, verbose)
  {
  ## wrapper to call the nlme function
  var.function.m1 = L1$var.function.m1
  var.function.m2 = L1$var.function.m2
  var.function.m3 = L1$var.function.m3
  var.function.m4 = L1$var.function.m4
  var.function.s1 = L1$var.function.s1
  var.function.s2 = L1$var.function.s2
  var.function.s3 = L1$var.function.s3
  var.function.s4 = L1$var.function.s4
  
  
  
  try(
  withCallingHandlers(
    {
    mv1 <- nlme(formu,
      data= dat1,
      fixed= as.list(fixedEffects),
      start= starting,
      random= randomEffects,
      weights= eval(parse(text=paste("var.function.", substr(analysis[1],1,1), varst, sep=""))),
      method = method)
    anova(mv1)
    }, warning = NlmeIfError), silent=TRUE)

  if(exists("mv1"))
    {
    # check lambda parameter for failure
    a = summary(mv1)
    a = a$tTable
    if (chk=="with.lambda")
      {
      if (a[4,5] > .05 | a[4,1] < 0)
        {
        if (verbose==TRUE) print (a)
        rm(mv1)
        if (verbose==TRUE) writeLines("\n ***** failed (lambda p> .05 or lambda < 0) *****\n")
        }
      }
    }
  
  if(exists("mv1")) return(mv1)  
  }

