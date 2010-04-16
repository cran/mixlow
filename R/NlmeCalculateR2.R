`NlmeCalculateR2` <-
function(mbest, dat1, verbose)
  {
  ## calculates R^2 and rc values, where rc= average model concordance correlation
  pred0 = mbest$fitted
  rcdrug = numeric(0)
  r2drug = numeric(0)
  u = as.vector(unique(dat1$drug))
  
  for (dr in u) {
    yhat = as.vector(predict(mbest)[dat1$drug==dr])
    y = dat1$adj_resp[dat1$drug==dr]
    y_mean = mean(y)
    yhat_mean = mean(yhat)
    
    resu = (y-yhat)^2
    resu = cbind(resu, (y-y_mean)^2)
    resu = cbind(resu, (yhat-yhat_mean)^2)

    # do rc
    n1 = sum(resu[,1])
    n2 = sum(resu[,2])
    n3 = sum(resu[,3])
    n0 = length(resu[,1])
    rc = 1 - n1/(n2 + n3 + n0*(y_mean-yhat_mean)^2)
    strr = paste("\ndrug = ", dr, "  rc = ",rc, "\n",sep="",collapse="")
    if (verbose==TRUE) writeLines(strr)
    rcdrug = c(rcdrug,rc)
    
    # do R2
    a <- cor(y, yhat)
    r2 = a^2
    r2drug = c(r2drug,r2)
    }
  
  return (list(rcdrug=rcdrug, r2drug=r2drug))
  }

