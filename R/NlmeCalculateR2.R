`NlmeCalculateR2` <-
function(mbest, dat1, verbose)
  {
  ## calculates R^2 and rc values, where rc= average model concordance correlation
  pred0 = mbest$fitted
  rcdrug = numeric(0)
  r2drug = numeric(0)
  u = unique(dat1$drug)
  for (jj in u)
    {
    dr = levels(dat1$drug)[jj]
    resu <- cbind(dat1$adj_resp[dat1$drug==dr],pred0[dat1$drug==dr,2])
    ym = mean(resu[, 1])
    ymhat = mean(resu[, 2])
    resu = cbind(resu, (resu[,1]-resu[,2])^2)
    resu = cbind(resu, (resu[,1]-ym)^2)
    resu = cbind(resu, (resu[,2]-ymhat)^2)
    
    # do rc
    n1 = sum(resu[,3])
    n2 = sum(resu[,4])
    n3 = sum(resu[,5])
    n0 = length(resu[,1])
    rc = 1 - n1/(n2 + n3 + n0*(ym-ymhat)^2)
    strr = paste("drug = ", dr, "  rc = ",rc, sep="",collapse="")
    if (verbose==TRUE) print(strr)
    rcdrug = c(rcdrug,rc)
    
    # do R2
    a <- cor(resu[, c(1:2)])
    r2 = a[1,2]^2
    r2drug = c(r2drug,r2)
    }
  return (list(rcdrug=rcdrug, r2drug=r2drug))
  }

