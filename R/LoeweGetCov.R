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
    
    if (length(lam0)==0) ncov = (2*(ndrugs+1)+1)^2
    if (length(lam0)> 0) ncov = (3*(ndrugs+1)+1)^2
    
    Names = names(nlme.results[[loc]])
    covv = numeric(0)
    for (i in seq(1,length(Names))) if (substr(Names[i],1,4) == "covv") covv = c(covv,nlme.results[[loc]][i])
    covv = as.numeric(covv)
    n = sqrt(length(covv))
    covv = array(covv,c(n,n))

    # delete row/column for u
    if (length(lam0)==0) {
      covv = covv[,-n]
      covv = covv[-n,]
      }
    if (length(lam0)>0) {
      covv = covv[,-(n-(ndrugs+1))]
      covv = covv[-(n-(ndrugs+1)),]
      }

    crit = qt(.975,degf)
    }
  
  
  # ---------------------------------------------------------------------------------
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
      # make sure all drugs have a lambda
      if (length(lam00) == 0) lam00 = 0
      lam0 = c(lam0, lam00)
      df0 = c(df0, eval(parse(text= paste("nlme.results[[",i,"]]$df.X.g",sep=""))))


      Names = names(nlme.results[[i]])
      covv = numeric(0)
      for (j in seq(1,length(Names))) if (substr(Names[j],1,4) == "covv") covv = c(covv,nlme.results[[i]][j])
      covv0 = as.numeric(covv)
      n = sqrt(length(covv0))
      covv0 = array(covv0,c(n,n))

      # drop u, 3rd row in either case
      covv0 = covv0[,-3]
      covv0 = covv0[-3,]
      n = sqrt(length(covv0))
      
      # copy covv0 into full size cov matrix
      covv = matrix(0,nrow=3,ncol=3)
      covv[1:n,1:n] = covv0
      
      covvL[[cnt]] = covv
      }
    degf = min(df0) # using minimum of df's for df
    crit = qt(.975,degf)
    
    
    if (all(lam0==0)) {
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

    

    if (all(lam0==0)==FALSE) {
      covv = matrix(data = 0.0, nrow = (3*(ndrugs+1)) , ncol = (3*(ndrugs+1)), byrow = FALSE, dimnames = NULL)
     
      for (j in seq(0,(length(covvL)-1))){
        m1 =  array(covvL[[j+1]],c(sqrt(length(covvL[[j+1]])),sqrt(length(covvL[[j+1]]))))    
        covv[((j*3+1):(j*3+3)), ((j*3+1):(j*3+3))] =  m1[1:3,1:3]
        }
      # change order in cov matrix to g1,g2,g3...p1,p2,p3...,lam1,lam2,lam3... 
      N = length(covvL)*3
      s1 = c(seq(1,N,by=3), seq(2,N,by=3), seq(3,N,by=3))
      covv = covv[s1,s1]
      }    
    
    } 
  
  return(list(degf=degf, covv=covv, crit=crit, g0=g0, p0=p0, lam0=lam0,u0=u0))  
  }

