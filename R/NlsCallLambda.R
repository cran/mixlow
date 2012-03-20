`NlsCallLambda` <-
function(numRandomTries, p.ran, g.ran, x, y, u0, tr, lambda0, verbose) 
  {
  ## calls the NLS function when the model includes a lambda parameter
  ##library(nlme)
  hits = data.frame(g=numeric(), p=numeric(), bic=numeric())
  flg = 0
  cnt = 0
  
  # loop thought NLS model for each random set of p and g
  for (j in seq(1,numRandomTries+1)) {
    try(
    {
    p1 = p.ran[j]
    g1 = g.ran[j]
    model.lambda <- nls(y~  exp(u1) *ifelse(x==0, 1, 1/(1+(exp(log(x) - log(exp(p1))))^exp(g1))) + exp(u1)*lambda, start= c(g1 = g1, 
      p1=p1, u1=u0, lambda=lambda0))
    cnt = cnt +1
    flg = 1
    bic = BIC(model.lambda)
    hits= rbind(hits, c(g1,p1,bic))
    }, silent=TRUE)
    }
  names(hits) = c("param.g.random", "param.p.random","BIC")
  if (length(hits$BIC) == 0) {
    stop("** The nls algorithm could not converge (for param.lambda >=0), try using different starting values")
    }
  
  # print results
  if (verbose==TRUE) writeLines(paste("\n  tray= ",tr, ",  nls starting params and BIC for different random starting values (param.lambda>=0):",sep=""))
  if (verbose==TRUE) print (format(hits, digits = 4, width = 8))
  
  # obtain best starting values from top-scoring random-start models
  minn = min(hits$BIC)
  hits2 = hits[abs(hits$BIC-minn) < abs(minn/20),]
  if (verbose==TRUE) writeLines ("\n mean parameters where BIC is near minimum")
  #hits2.mean = mean(hits2)
  hits2.mean = colMeans(hits2)
  if (verbose==TRUE) print (hits2.mean)
  
  # call NLS model using best starting values
  g1 = as.vector(hits2.mean[1])
  p1 = as.vector(hits2.mean[2])
  model.lambda <- nls(y~  exp(u1) *ifelse(x==0, 1, 1/(1+(exp(log(x) - log(exp(p1))))^exp(g1))) + exp(u1)*lambda, start= c(g1 = g1, p1=p1, 
    u1=u0, lambda=lambda0))
  
  bic.lambda = BIC(model.lambda)
  sum.model.lambda = summary(model.lambda)
  if (verbose==TRUE) writeLines("\n-- lambda>=0 model --")
  if (verbose==TRUE) print (sum.model.lambda$parameters)
  if (verbose==TRUE) writeLines("\n")
  
  return(list(sum.model.lambda=sum.model.lambda, model.lambda=model.lambda, bic.lambda=bic.lambda))
  }

