`doNls` <-
function(mixlowData, parameterDefaults, lambdaThreshold=0.05, numRandomTries=10, verbose=FALSE)  {
  ## outer loop that calls the functions that call the NLS function for each analysis desired

  if (!inherits(mixlowData, "mixlowData")) 
    stop("use only with \"mixlowData\" objects")

  nls.control(maxiter = 500)
  
  trays = parameterDefaults$tray
  crData = mixlowData$concentrationResponse
  drugRatios = mixlowData$drugRatios
  
  v = numeric(length(trays))
  a = character(length(trays))
  nlsEstimates = data.frame(tray=a, drug=a, cell=a, g=v, p=v,u=v,lambda=v)
  nlsGraphing = vector(mode="list", length = length(trays))
  trayCnt = 0
  
  # ----------------- loop through each tray and call NLS model -----------------
  for (tr in trays) {
    if (verbose==TRUE) writeLines(paste("\n---- NLS assessment of tray ",tr," ----\n",sep=""))
    trayCnt = trayCnt + 1
    drug = drugRatios$drug[drugRatios$tray==tr]
    cell = drugRatios$cell[drugRatios$tray==tr]
    Units = drugRatios$Units[drugRatios$tray==tr]
    
    # setup data for NLS
    dat2 = crData[crData$tray==tr & crData$label == "rx",]
    cntrl = mean(dat2$adj_resp[dat2$conc==0])
    dat2$scaled = dat2$adj_resp/cntrl
      
    # setup xx to be used for graphing predictions
    x2 = dat2$adj_conc
    y = dat2$adj_resp
    x = dat2$conc
    xx = exp(seq(log(min(x2)),1.1*log(max(x2)),length.out = 500))
    
    # collect response values to find response closest to 1/2 average control-well response
    tmp = as.data.frame(cbind(x= dat2$adj_conc[dat2$tray==tr], y = dat2$adj_resp[dat2$tray==tr]))
    tmp = aggregate(tmp$y, list(X=tmp$x) ,mean)
    names(tmp) <- c("x","y")
    # default: choose IC50s by concentration nearest to 1/2 the control response
    idx = which.min(abs(cntrl/2 - tmp$y))
    
    # default g = 1, lambda = 0  --- use log of g and p
    p0 = log(as.numeric(as.vector(tmp$x[idx])))
    g0 = log(1)
    lambda0 = NaN
    p0.std =  abs(p0)/2
    g0.std = 1
    u0 = log(cntrl)
    
    # check for user-supplied parameter values for this tray
    if (is.nan(parameterDefaults[trayCnt,]$param.g) ==FALSE) g0 = parameterDefaults[trayCnt,]$param.g
    if (is.nan(parameterDefaults[trayCnt,]$param.p) ==FALSE) p0 = parameterDefaults[trayCnt,]$param.p
    if (is.nan(parameterDefaults[trayCnt,]$param.lambda) ==FALSE) lambda0 = parameterDefaults[trayCnt,]$param.lambda
    if (is.nan(parameterDefaults[trayCnt,]$param.g.std) ==FALSE) g0.std = parameterDefaults[trayCnt,]$param.g.std
    if (is.nan(parameterDefaults[trayCnt,]$param.p.std) ==FALSE) p0.std = parameterDefaults[trayCnt,]$param.p.std
    
    # make list of random starting values for p and g
    p.ran = c(p0, rnorm(numRandomTries,p0,p0.std))
    g.ran = c(g0, abs(rnorm(numRandomTries,g0, g0.std)))
    
    # first estimate NLS model with lambda=0
    L1 = NlsCallNoLambda (numRandomTries, p.ran, g.ran, x, y, u0, tr, verbose) 
    sum.model.no.lambda = L1$sum.model.no.lambda
    bic.no.lambda = L1$bic.no.lambda
    model.no.lambda = L1$model.no.lambda
    
    #defaults (use no.lambda model unless lambda model is superior)
    msg = "use model with lambda=0"
    sum.model.nls = sum.model.no.lambda  
    sum.model.nls$parameters = rbind(sum.model.nls$parameters, c(0,0,0,0))
    pred = predict(model.no.lambda,list(x=xx))
    pred.orig = predict(model.no.lambda)
  
    # now estimate NLS model with lambda != 0
    if (is.nan(lambda0)==FALSE) 
      {
      L2 = NlsCallLambda (numRandomTries, p.ran, g.ran, x, y, u0, tr, lambda0, verbose) 
      sum.model.lambda = L2$sum.model.lambda
      bic.lambda = L2$bic.lambda
      model.lambda = L2$model.lambda
      
      # choose best model between lambda and no.lambda
      if ((sum.model.lambda$parameters[4,1] > as.vector(lambdaThreshold)) & bic.lambda < bic.no.lambda )  
        {
        sum.model.nls = sum.model.lambda
        pred = predict(model.lambda,list(x=xx))
        pred.orig = predict(model.lambda)
        msg = "use model with lambda >0"
        }
      }
    if (verbose==TRUE) writeLines(msg)
    
    # save results in matrix nlsEstimates
    g = sum.model.nls$parameters[1,1]
    p = sum.model.nls$parameters[2,1]
    u = sum.model.nls$parameters[3,1]
    lambda = sum.model.nls$parameters[4,1]
    levels(nlsEstimates$tray) = c(levels(nlsEstimates$tray),as.vector(tr))
    nlsEstimates[trayCnt,1] = tr
    levels(nlsEstimates$drug) = c(levels(nlsEstimates$drug),as.vector(drug))
    nlsEstimates[trayCnt,2] = drug
    levels(nlsEstimates$cell) = c(levels(nlsEstimates$cell),as.vector(cell))
    nlsEstimates[trayCnt,3] = cell
    nlsEstimates[trayCnt,4] = g
    nlsEstimates[trayCnt,5] = p
    nlsEstimates[trayCnt,6] = u
    nlsEstimates[trayCnt,7] = lambda
    
    # save results for graphing
    nlsGraphing[[trayCnt]] = list(tray= tr, drug=drug, xx=xx, pred=pred, y=y, x2=x2)
    }

  names(nlsEstimates) <- c("tray","drug","cell","g","p","u","lambda")
  if (verbose==TRUE) writeLines("\n---- Final results of NLS assessment ----")
  if (verbose==TRUE) print (nlsEstimates)
  if (verbose==TRUE) writeLines("\n")

  returnList = list(nlsEstimates=nlsEstimates, nlsGraphing=nlsGraphing)
  class(returnList) <- "nlsData"
  return (returnList) 
  }

