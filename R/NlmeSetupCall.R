`NlmeSetupCall` <-
function(ord, paramList, dat1, method, varFunction, cell, analysis, verbose)
  {
  ## collects parameters and sets up a call to the NLME function
  library(nlme)
  
  # make drug string
  ord2 = paste(ord,"+",sep="",collapse="")
  ord2 = substr(ord2,1,nchar(ord2)-1)
  
  # keep track of lowest BIC
  lowBIC <- 1e12

  # if parameter lambda>0 for any drug, then use lambda model, otherwise use no.lambda model
  fFlag = any(paramList$lambda > 0)
  if (fFlag == FALSE) paramList$lambda = NULL
  
  # set up parameter list
  paramList$drug = NULL
  
  # setup error (variance) functions, and fixed-effects and random effects lists 
  L1 = NlmeSetupVarFunction (ord, dat1)
  fixedEffects = L1$fixedEffects
  randomEffects = L1$randomEffects

  # loop through NLME call using lambda model, then with no.lambda model unless fFlag == FALSE 
  for (chk in c("with.lambda", "no.lambda"))
    {
    if (chk == "with.lambda" & fFlag == FALSE) next
    if (verbose==TRUE) writeLines(paste("\n============  using: ",chk," ============\n",sep=""))
    
    # make formula and adjust parameter and fixed-effects lists, and make list for starting values, all for no.lambda model
    formu = as.formula(adj_resp ~ ifelse(conc>0, (exp(u)* 1/(1+(exp(log(conc)-p))^exp(g)) ), exp(u)  ))
    param2 = paramList
    param2$lambda = NULL
    starting = unlist(param2)
    fixedEffects[[4]] = NULL
    if (chk == "with.lambda")
      {
      # make formula etc. if lambda model is used
      starting = unlist(paramList)
      fixedEffects[[4]] <- as.formula(lambda~drg-1)
      if (length(ord) == 1) fixedEffects[[4]] <- as.formula(lambda~1)
      formu = as.formula(adj_resp ~ ifelse(conc>0, ((1-lambda)*exp(u)* 1/(1+(exp(log(conc)-p))^exp(g)) + exp(u)*lambda), 
        ((1-lambda)*exp(u) + exp(u)*lambda)  ))
      }
    
                                                  
    # ------------- loop through NLME call for each error function desired 
    if (is.list(varFunction))
      varFunction = as.numeric(unlist(varFunction))
      
    for (varst in varFunction)
      {
      # run model
      if (verbose==TRUE) writeLines(paste("\n-------------------- variance function= ",varst,"--------------\n",sep=""))
      
      # remove object mv1 if present
      objs = objects()
      grp = grep("mv1",objs)
      grlen = length(grp)
      if (grlen>0) rm(mv1)
      rm(grp)

      if (verbose==TRUE) writeLines ("starting")
      if (verbose==TRUE) print(starting)
      
      
      # call NLME, return mv1=NULL if unsuccessful
      mv1 = NlmeCall(formu, dat1, fixedEffects, starting, randomEffects, varst, method, chk, L1, analysis, verbose)
      
      # save mv1 with new name and print summary of results if NLME call is successful
      if(length(mv1) > 0)
        {
        if (chk == "with.lambda") strr = paste("model_varFun",varst,"_lambda",sep="",collapse="")
        if (chk == "no.lambda") strr = paste("model_varFun",varst,"_no.lambda",sep="",collapse="")
        strr2 = paste(strr," = mv1")
        eval(parse(text=strr2))
        if (verbose==TRUE) writeLines(c("\n************************",strr,"*********************************\n"))
        if (verbose==TRUE) print(summary(mv1))
        if (verbose==TRUE) writeLines("\n")
        if (verbose==TRUE) print("coef")
        if (verbose==TRUE) print(mv1$coefficients)
        # check F values
        ano <- anova(mv1)
        if (verbose==TRUE) print(ano)
        }
      
      # obtain new low BIC if NLME call is successful
      if(length(mv1) > 0)
        {
        if (BIC(mv1) < lowBIC)
          {
          if (verbose==TRUE) writeLines(paste("\n new BIC = ",BIC(mv1),"\n",sep=""))
          lowBIC <- BIC(mv1)
          }
        }
      }
    }

  # remove mv1 object if it exists 
  objs = objects()
  grp = grep("mv1",objs)
  grlen = length(grp)
  if (grlen>0) rm(mv1)
  rm(grp)
  
  if (verbose==TRUE) writeLines("\n\n********")
  if (verbose==TRUE) writeLines(paste("---------- finished all models, drugs= ", ord2, "  ----------\n",sep = ""))
  
  # find best model by searching for saved model objects, then comparing BIC values to best BIC
  good <- ls(pattern="^model_varFun")
  if (length(good)<1)
    {
    if (verbose==TRUE) print("no models were suitable")
    best = 0
    mbest = 0
    if (verbose==TRUE) writeLines(" ** best= NULL **  ")
    return (list(best=best, mbest=mbest, rcdrug=0, r2drug=0 ))
    }
  {if (length(good)>1) {
    if (verbose==TRUE) writeLines("\n   BIC criterion:")
    for (g in good){
      if (eval(parse(text=paste("BIC(",g,")",sep=""))) == lowBIC) best = g
      if (verbose==TRUE) writeLines(paste("    model= ",g,"     BIC= ",eval(parse(text=paste("BIC(",g,")",sep=""))),sep=""))
      }
          
    if (verbose==TRUE) writeLines(paste("\n----- best= ", best," -----\n", sep="",collapse=""))
    }
  else
    {
    # length(good)==1
    best <- good          
    if (verbose==TRUE) writeLines(paste("\n----- best= ",best, " -----\n", sep="",collapse=""))
    } }
  
  # create object mbest for best model
  mbest = eval(parse(text=best))  

  if (verbose==TRUE) print(summary(mbest))
  if (verbose==TRUE) writeLines("\n")
  if (verbose==TRUE) writeLines("-- coeffficients --\n")
  if (verbose==TRUE) print(mbest$coefficients)
  if (verbose==TRUE) writeLines("\n-- covariance --\n")
  if (verbose==TRUE) print(mbest$varFix)
  if (verbose==TRUE) writeLines("\n")
  
  ano <- anova(mbest)
  if (verbose==TRUE) print(ano)

  if (verbose==TRUE) writeLines("\n")

  # calculate R^2 values and rc goodness-of-fit statistic
  L3 = NlmeCalculateR2(mbest, dat1, verbose)
  
  pred0 = mbest$fitted
  residu = resid(mbest,type="p")

  if (verbose==TRUE) writeLines("********\n")

  return (list(best=best, mbest=mbest, rcdrug=L3$rcdrug, r2drug=L3$r2drug, pred0=pred0, residu=residu ))
  }

