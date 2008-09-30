`LoeweGetCov` <-
function(nlme.results, drugs, loc, ndrugs )
  {
  ## collects the NLME parameter covariance matrix for use in Loewe analysis
  # all drugs analyzed together?
  if (length(loc) == 1) {
    # all drugs were analyzed together
    uu = nlme.results[[loc]]$coeff.fixed.u
    degf = nlme.results[[loc]]$df.X.g.drg1
    g0 = numeric()
    p0 = numeric()
    lam0 = numeric()
    u0 = nlme.results[[loc]]$coeff.fixed.u
    for (i in seq(1,length(drugs)+1)) {
      g0 = c(g0, eval(parse(text= paste("nlme.results[[",loc,"]]$coeff.fixed.g.drg",i,sep=""))))
      p0 = c(p0, eval(parse(text= paste("nlme.results[[",loc,"]]$coeff.fixed.p.drg",i,sep=""))))
      lam0 = c(lam0, eval(parse(text= paste("nlme.results[[",loc,"]]$coeff.fixed.lambda.drg",i,sep=""))))
        }
    # values for f are not needed in Loewe analysis
    if (length(lam0)==0) ncov = (2*(ndrugs+1)+1)^2
    if (length(lam0)> 0) ncov = (3*(ndrugs+1)+1)^2
    covv = nlme.results[[loc]]$covv1
    for (i in seq(2,ncov)) {
      strr = paste("covv <- c(covv, nlme.results[[",loc,"]]$covv",i,")", sep="")
      eval(parse(text=strr))
      }
    covv = array(covv,c(sqrt(ncov),sqrt(ncov)))
    covv = covv[1:(2*(ndrugs+1)),1:(2*(ndrugs+1))]
    crit = qt(.975,degf)
    }
  # drugs analyzed separately by nlme?
  if (length(loc) > 1) {
    cnt = 0
    covvL = list()
    g0 = numeric()
    p0 = numeric()
    u0 = numeric()
    lam0 = numeric()
    df0 = numeric()
    for (i in loc) {
      cnt = cnt+1
      lam00 = numeric()
      g0 = c(g0, eval(parse(text= paste("nlme.results[[",i,"]]$coeff.fixed.g",sep=""))))
      p0 = c(p0, eval(parse(text= paste("nlme.results[[",i,"]]$coeff.fixed.p",sep=""))))
      u0 = c(u0, eval(parse(text= paste("nlme.results[[",i,"]]$coeff.fixed.u",sep=""))))
      lam00 = eval(parse(text= paste("nlme.results[[",i,"]]$coeff.fixed.lambda",sep="")))
      lam0 = c(lam0, eval(parse(text= paste("nlme.results[[",i,"]]$coeff.fixed.lambda",sep=""))))
      df0 = c(df0, eval(parse(text= paste("nlme.results[[",i,"]]$df.X.g",sep=""))))
      if (length(lam00) == 0) ncov = 9
      if (length(lam00) == 1) ncov = 16
      covv0 = nlme.results[[i]]$covv1
      for (j in seq(2,ncov)){
        strr = paste("covv0 <- c(covv0, nlme.results[[",i,"]]$covv",j,")", sep="")
        eval(parse(text=strr)) 
        }
      covvL[[cnt]] = covv0
      }
    degf = min(df0) # using minimum of df's for df
    crit = qt(.975,degf)
    
    covv = matrix(data = 0.0, nrow = (2*(ndrugs+1)) , ncol = (2*(ndrugs+1)), byrow = FALSE, dimnames = NULL)
    for (j in seq(0,(length(covvL)-1))){
      m1 =  array(covvL[[j+1]],c(sqrt(length(covvL[[j+1]])),sqrt(length(covvL[[j+1]]))))    
      covv[((j*2+1):(j*2+2)), ((j*2+1):(j*2+2))] =  m1[1:2,1:2]
      }
    # change order in cov matrix to g1,g2,g3...p1,p2,p3...
    N = length(covvL)*2
    s1 = c(seq(1,N,by=2), seq(2,N,by=2))
    covv = covv[s1,s1]
    } 
  return(list(degf=degf, covv=covv, crit=crit, g0=g0, p0=p0, lam0=lam0,u0=u0))  
  }

