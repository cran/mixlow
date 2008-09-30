`plotLoeweData` <-
function(loeweData, ...) {
  ## graphs results from the Loewe analysis
  if (!inherits(loeweData, "loeweData")) 
    stop("use only with \"loeweData\" objects")
  
  arglist = list(...)    
  
  dat0 = loeweData$dat0
  mix = loeweData$mix
  ciL = loeweData$ciL
  psi = loeweData$p0
  gamm = loeweData$g0
  drugs = loeweData$drugs
  uu = loeweData$uu
  Units = loeweData$Units
  

  # make plots ------------------------------------------
  tit = paste("Interaction Index\nMixture= ", mix,sep="")
  #if (formatFlags$titles == 0) tit = ""
  yLab = "Estimated Loewe Index"
  xLab = "Fraction Affected"
  if (is.na(match("xlab", names(arglist))) == FALSE)
      xLab = arglist$xlab
  if (is.na(match("ylab", names(arglist))) == FALSE)
      yLab = arglist$ylab

  if (is.na(match("main", names(arglist))) == FALSE)
    tit = arglist$main

  cols = rep(c("black", "blue", "red", "purple", "green", "cyan"),4,each=1)
  typs = rep(1:6,4)

  plot(ciL[,2] ~ ciL[,1], type="l", col=cols[1],
    ylim=c(min(ciL[,4]),max(ciL[,5])), main= tit,
    ylab= yLab, xlab= xLab,lty=1)
  lines(ciL[,4] ~ ciL[,1], lty=2,col= cols[3],type="l")
  lines(ciL[,5] ~ ciL[,1], lty=2,col= cols[3], type="l")
  abline(1,0,lwd=2,col=cols[2])
  grid(col= "black", lwd = 1)
  
  leg = c("Loewe Index","Confidence Intervals","Additivity Reference Line")
  if (is.na(match("legend", names(arglist))) == FALSE) {
    leg = arglist$legend
    if (length(leg) == 3) 
      legend(x="topright", leg, lty=c(1,2,1),
    lwd=c(1,1,2),col=c(cols[1],cols[3],cols[2]))
    }
  if (is.na(match("legend", names(arglist))) == TRUE)
    legend(x="topright", leg, lty=c(1,2,1),
    lwd=c(1,1,2),col=c(cols[1],cols[3],cols[2]))
  
  
  cols = rep(c("black", "blue", "red", "purple", "green", "cyan"),4,each=1)
  typs = rep(1:6,4)
  xx = exp(seq(log(min(dat0$conc[dat0$conc > 0])),1.1*log(max(dat0$conc[dat0$conc > 0])),length.out = 500))
  
  y1 = exp(uu)/(1+ (exp(log(xx)-psi[1]))^exp(gamm[1]))
  
  tit = paste("Estimated Concentration-Response:\nMixture= ", mix, sep="")
  yLab = "Response"
  xLab = paste("Concentration, w/adj for zero, ", Units, sep="")

  plot( y1 ~ xx,
    main= tit,
    ylab= yLab, xlab= xLab, log="x", type="l", lty=typs[1], col=cols[1])
  for (k in seq(2,length(psi))) 
    {
    yk = exp(uu)/(1+ (exp(log(xx)-psi[k]))^exp(gamm[k]))
    lines(yk ~ xx, lty=typs[k], col=cols[k], type= "l", lwd=2)
    }
  grid(col= "black", lwd = 1)
  nam = c(drugs,mix)
  legend(x="bottomleft", legend=nam, lty=typs[1:length(psi)], col= cols[1:length(psi)])
  

  return ()
  }

