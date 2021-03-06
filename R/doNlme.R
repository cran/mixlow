`doNlme` <-
function(mixlowData, nlsData, drugs=getDrugs(mixlowData), analysis="multiple", 
  varFunction=1, method="ML", verbose=FALSE) {
  ## outer loop that sets up calls to NLME

  if (!inherits(nlsData, "nlsData")) 
    stop("use only with \"nlsData\" objects")  
  
  if (analysis=="single") {
    if (is.list(varFunction)==FALSE)        stop("VarFunction must be a list of vectors named by drug")
    if (length(names(varFunction)) == 0)    stop("VarFunction must be a list of vectors named by drug") 
    }
  varFunction0 = varFunction
  
  nls.estimates = nlsData$nlsEstimates
  data00 = mixlowData$concentrationResponse
  drugRatios = mixlowData$drugRatios
  
  # check for number of drugs
  newDrugVector = drugs[1]   # default, for single drug
  
  if (length(drugs) !=1)  {
    # length will be 1 if a single drug is being analyzed alone.  Else, the drug set will contain a mixture.

    # order drugs so mixture is last
    newDrugVector = numeric(0)
    mixFlag = 0
    for (i in seq(1,length(drugs))) {
      dr = drugs[i]
      ratio = drugRatios[,dr]

      
      if (all(ratio==0)) {
        if (mixFlag == 1)
          stop("Drugs for analysis can contain only one mixture")
        mixFlag = 1
        mix = dr
        }
      if (any(ratio>0))
        newDrugVector = c(newDrugVector,dr)
      }
    
    if (mixFlag != 0) {
      #stop("Drugs must contain at least two drugs and one mixture.  No mixture present.")
      newDrugVector = c(newDrugVector,mix)
      }

    }    

  drugs0 = newDrugVector
  

       
  # remove entries for blanks
  data00 <- data00[data00$label=="rx",]

  nlmeResults = list(0)
  nlmeGraph = list(0)
  nlmeModels = list(0)
  
  numberOfAnalysis = 1
  if (analysis == "single")
    numberOfAnalysis = length(drugs0)
    
  
  
  # loop over each set of drugs to be modeled ----------------------------------------------------------------------
  for (setNum in seq(1,numberOfAnalysis)) {
    
    drugs = drugs0
    if (numberOfAnalysis > 1) {
      drugs = drugs0[setNum]
      varFunction = varFunction0[[drugs]]
      }

    #dr1 = drugs[1]
    cell = unique(as.vector(drugRatios$cell[drugRatios$drug==drugs[1]]))
        
    useTrays = drugRatios[(drugRatios$drug %in% drugs),]
        
    # add drug and cell data to data0
    data0 = merge(data00, useTrays, by= "tray")
    data0$tray = as.factor(data0$tray)
    
    # collect data for drugs in desired set
    inlist <- data0$drug %in% drugs
    dat1 <- data0[inlist & data0$cell== cell, ]
    
    # set order for assessing drugs
    ord = drugs
    ord2 = paste(ord,"|",sep="",collapse="")
    ord2 = substr(ord2,1,nchar(ord2)-1)
    
    if (verbose==TRUE) writeLines("######################################################################\n")
    if (verbose==TRUE) writeLines(paste("\n **** NLME analysis, drugs for round ", setNum, " are: ", ord2, "  ****",sep=""))
    
    # if only one tray, then duplicate data
    new = 20000
    for (dr in ord)
      {
      utray = unique(dat1$tray[dat1$drug==dr])
      if (length(utray) == 1)
        {
        if (verbose==TRUE) writeLines("\n********** duplicating data, only one tray present ***************\n")
        tmp <- dat1[dat1$drug==dr,]
        tray = unique(tmp$tray)
        # assumes that no tray number is called new20000+  
        tmp$tray <- paste("new", new,sep="")
        new = new + 1
        dat1 <- rbind(dat1,tmp)
        }
      }

    # order drugs
    dat1 <- dat1[order(dat1$cell,dat1$drug, dat1$tray, dat1$conc),]
    dat1$drg <- rep(1,length(dat1$adj_resp))
    for (i in 1:length(ord)) dat1$drg[dat1$drug==ord[i]] <- i
    dat1$drg = as.factor(dat1$drg)
    
    # make parameter list
    tmp1 = unique(dat1[,c("drug","tray")])
    combo = merge(nls.estimates,tmp1, by= "tray")
    tmp5 = cbind(as.numeric(as.vector(combo$g)),as.numeric(as.vector(combo$p)),as.numeric(as.vector(combo$u)),as.numeric(as.vector(combo$lambda)))
    if (is.null(combo$drug)==FALSE) {param = aggregate(tmp5, list(drug=combo$drug) ,mean)}
    if (is.null(combo$drug) == TRUE) {param = aggregate(tmp5, list(drug=combo$drug.x) ,mean)}
    names(param) = c("drug","g","p","u","lambda")
    paramList = as.list(param)
    

    
    ord3 = rep(0,length(ord))
    for (i in 1:length(ord)) ord3[ord == paramList$drug[i]] <- i
    #for (i in 1:length(ord)) {ord3[paramList$drug==ord[i]] <- i
    

    paramList$drug = paramList$drug[ord3]
    paramList$g = paramList$g[ord3]
    paramList$p = paramList$p[ord3]
    paramList$u = paramList$u[ord3]
    paramList$lambda = paramList$lambda[ord3]
    
    
    if (length(paramList$drug) != length(ord)) {stop("** Mismatch in length of drugs and parameter estimates")}
    if (any(paramList$drug != ord)) {
      stop("** Drugs/param list in wrong order")}
    
    paramList$u = mean(as.numeric(as.vector(combo$u)))
    paramList$p = paramList$p # p is in log form
    paramList$g = paramList$g # g is in log form
    
    # setup call to NLME and then call NLME, return NLME results for best model
    meth = method
    

    best2 <- NlmeSetupCall(ord, paramList, dat1, method, varFunction, cell, analysis, verbose)
    
    if(best2$best==0) {
      if (verbose==TRUE) print('no models were suitable')
      next
      }

    # collect results
    temp <- as.list(c(nam =best2$best, method=meth, drug= ord2, setNum=setNum, cell= cell, lik= logLik(best2$mbest),
        bic = BIC(best2$mbest), sig = best2$mbest$sigma, modelstruct = unlist(best2$mbest$modelStruct),
        rc = best2$rcdrug, r2= best2$r2drug, se=unlist(summary(best2$mbest)$tTable[,2]),
        coeff = unlist(best2$mbest$coefficients), df = unlist(best2$mbest$fixDF)[1], covv=as.list(best2$mbest$varFix)   ))
    nlmeResults[[setNum]] <- temp
    
    # collect results for graphing
    nlmeGraph[[setNum]] = list(pred0=best2$pred0, dat1=dat1, ord=ord, residu=best2$residu, best=best2$best)
    
    # alter formula so that anova can be called on models
    #best2$mbest$call$model = y~1
    nlmeModels[[setNum]] = best2$mbest

    }
  returnList = list(nlmeResults=nlmeResults, nlmeGraph=nlmeGraph, nlmeModels=nlmeModels)
  class(returnList) <- c("nlmeData")
  return (returnList)

  }

