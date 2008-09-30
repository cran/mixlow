`doLoewe` <-
function(mixlowData, nlmeData, verbose=FALSE)
  {
  ## calculates the drug interaction Loewe index  

  if (!inherits(nlmeData, "nlmeData")) 
    stop("use only with \"nlmeData\" objects")

  # remove entries for blanks
  dat0 <- mixlowData$concentrationResponse
  dat0 = dat0[dat0$label=="rx",]
  drugRatios = mixlowData$drugRatios
  
  nlmeResults = nlmeData$nlmeResults
  
  # add drug and cell data to dat0
  dat0 = merge(dat0, drugRatios, by= "tray")

  # collect data for Loewe analysis
  mixtureList = LoeweCollectData(nlmeResults, drugRatios)
  
  # do each mixture set -----------------------------------------------------------
  loeweResults = list()
  drugs = mixtureList$drugs
  fractions = mixtureList$fractions
  locations = mixtureList$locations
  mix = mixtureList$mixture
  ndrugs = length(drugs) 

  if ( (mixtureList$analysis == "multiple") & (length(locations) > 1) ) 
    stop("Analysis=multiple but more than one NLME result is being used")
  if ( (mixtureList$analysis == "single") & (length(locations) != ndrugs+1) ) 
    stop("Analysis=single but wrong number of NLME results are being used")
  if ( is.null(locations)) 
    stop("Cannot find needed set of NLME results.  Check single vs. multiple analysis.")
  
  # get covariance matrix
  L3 = LoeweGetCov(nlmeResults, drugs, locations, ndrugs)
  degf = L3$degf
  covv = L3$covv
  crit = L3$crit 
  g0 = L3$g0
  p0 = L3$p0
  lam0 = L3$lam0
  u0 = L3$u0
  
  if (verbose==TRUE) writeLines("\n-----------------------------------------------------------------------------")
  if (verbose==TRUE) writeLines(paste("    Mixture= ", mix, "\n", sep=""))
  if (verbose==TRUE) writeLines("    drugs= ")
  if (verbose==TRUE) print (drugs)
  if (verbose==TRUE) writeLines("\n    fractions= ")
  if (verbose==TRUE) print (fractions)
  if (verbose==TRUE) writeLines("\n    NLME covariance (ordered by g1,g2,g3...,p1,p2,p3...)= ")
  if (verbose==TRUE) print (covv)    
  if (verbose==TRUE) writeLines("\n") 
  
  
  # loop through algorithm, write results to ciL  ----------------------------------------------
  gamm = g0
  psi = p0 
  uu = mean(u0)
  
  iii = seq(.02, 0.98, 0.01)
  ciL = matrix(nrow= length(iii), ncol=5)
  
  cnt = 0
  # todo: allow for more than one mixture
  for (phi in iii)
    {
    k= phi/(1-phi)
    # *************** Log version ********************************
    cnt = cnt+1
    Dgamm = rep(0.0, length(drugs)+1 )
    Dpsi =  rep(0.0, length(drugs)+1 )
  
    #gamm = c(.12, .21, .32)
    #psi = c(.34, .43, .54)
    #fractions = c(.5, .5)
    #k=.6
  
    # first do diff(logL,gamm)
    for (r in seq(1, length(drugs)+1))
      {
      dgamm = rep(0.0, length(drugs)+1 )
      dgamm[r] = 1.0
      numer = 0
      for (j in seq(1, length(drugs)))
        {
        numer = numer + fractions[j]*(-1/exp(gamm[ndrugs+1])*dgamm[ndrugs+1]*log(k) + 1/exp(gamm[j])*dgamm[j]*log(k))*
        exp(log( k^(1/exp(gamm[ndrugs+1])) ) + psi[ndrugs+1]- log( k^(1/exp(gamm[j])) )-psi[j])
        }
      denom = 0
      for (j in seq(1, length(drugs)))
        {
        denom = denom + fractions[j]*exp(log( k^(1/exp(gamm[ndrugs+1])) ) + psi[ndrugs+1]-log( k^(1/exp(gamm[j])) ) -psi[j])
        }
      Dgamm[r] = numer/denom
      }
  
  
    # second do diff(logL,psi)
    for (r in seq(1, length(drugs)+1))
      {
      dpsi = rep(0.0, length(drugs)+1 )
      dpsi[r] = 1.0
      numer = 0
      for (j in seq(1, length(drugs)))
        {
        numer = numer + fractions[j]*(dpsi[ndrugs+1]-dpsi[j])*exp(log( k^(1/exp(gamm[ndrugs+1])) )+ psi[ndrugs+1]-log( k^(1/exp(gamm[j]))) -psi[j])
        }
      denom = 0
      for (j in seq(1, length(drugs)))
        {
        denom = denom + fractions[j]*exp(log( k^(1/exp(gamm[ndrugs+1])) ) +psi[ndrugs+1]-log( k^(1/exp(gamm[j])) )-psi[j])
        }
      Dpsi[r] = numer/denom
      }
  
    # put DV and Dw together
    Dr = c(Dgamm, Dpsi)
    seL = sqrt(t(Dr) %*% covv %*% Dr)
  
    # calculate loewe
    summ = 0
    for (j in seq(1, length(drugs)))
      {
      summ = summ + fractions[j]*exp(log( k^(1/exp(gamm[ndrugs+1])) )+ psi[ndrugs+1]-log( k^(1/exp(gamm[j])) ) -psi[j])
      }
    L = summ
    Llog = log(L)
  
  
    ciL[cnt,] = c(phi, L, seL, exp(Llog-crit*seL), exp(Llog+seL*crit))
    }
  ciL = as.data.frame(ciL)
  names(ciL) = c("Fraction.Affected","Log.Index", "SE.Log.Index", "Lower.CI","Upper.CI")
  
  if (verbose==TRUE) writeLines(paste("\n ------ Results for analysis, mixture= ", mix, " -----------------", sep=""))
  if (verbose==TRUE) print (ciL)      
  
  # calculate synergism/antagonism scores over interval
  L4 = LoeweCalculateScore(ciL)
  
  # collect results
  Units = unique(as.vector(drugRatios$Units))
  loeweResults = list(drugs=drugs, mix=mix, ciL = ciL, covv = covv, g0=g0, p0=p0, uu=uu, dat0=dat0, 
    score.interval=L4$score.interval, syner=L4$syner, syner.total=L4$syner2, antag=L4$antag, antag.total=L4$antag2, Units = Units )
  
  class(loeweResults) <- "loeweData"
  return (loeweResults)
  }

