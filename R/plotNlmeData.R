`plotNlmeData` <-
function(nlmeData, ask = prod(par("mfcol")) < length(nlmeData$nlmeGraph) && dev.interactive(), ...) {
  ## graphs results from the NLME analysis
  
  arglist = list(...)
  if (!inherits(nlmeData, "nlmeData")) 
    stop("use only with \"nlmeData\" objects")

  if (ask) {
    opar <- par(ask=TRUE)
    on.exit(par(opar))
    }
  
  nlmeGraph = nlmeData$nlmeGraph    
  
  # loop through each nlme result set -------------------
  for (iii in seq(1,length(nlmeGraph))) {  
    # collect predicted values, drug order, data, residuals, etc.
    pred0 = nlmeGraph[[iii]]$pred0
    wid = length(pred0[1,])
    ord = nlmeGraph[[iii]]$ord
    dat1 = nlmeGraph[[iii]]$dat1
    residu = nlmeGraph[[iii]]$residu
    best = nlmeGraph[[iii]]$best
    
    # make drug string
    ord2 = paste(ord,"+",sep="",collapse="")
    ord2 = substr(ord2,1,nchar(ord2)-1)
    
    # loop through each drug ----------------
    for (jj in 1:length(ord))
      {
      # get list of unique trays for current drug
      tmp2 = dat1[dat1$drug==ord[jj], ]
      utray = as.vector(unique(tmp2$tray))

      # make graph for each tray and each drug
      yLab = "Response"
      xLab = "Concentration, w/adj for zero"
      Cols = c("blue","red")
      if (is.na(match("xlab", names(arglist))) == FALSE)
          xLab = arglist$xlab
      if (is.na(match("ylab", names(arglist))) == FALSE)
          yLab = arglist$ylab

      for (ii in utray)
        {
        tit = paste("NLME: ",ord[jj], ", Tray= ",ii, "\nModel= ", best, sep="")
        
        if (is.na(match("main", names(arglist))) == FALSE)
          tit = arglist$main
        
        
        ran = c(tmp2$adj_resp[tmp2$tray==ii], pred0[tmp2$tray==ii,2], pred0[tmp2$tray==ii,1])
        
        plot(pred0[dat1$drug==ord[jj] & dat1$tray==ii, 1] ~ dat1$adj_conc[dat1$drug==ord[jj] & dat1$tray==ii],
          log="x",type="n", ylim=range(ran), xlim=range(dat1$adj_conc[dat1$drug==ord[jj] & dat1$tray==ii]),
          main= tit, xlab= xLab, ylab= yLab)
        
        points((dat1$adj_resp[dat1$drug==ord[jj] & dat1$tray==ii])~dat1$adj_conc[dat1$drug==ord[jj] & dat1$tray==ii], pch=1, col=Cols[1])
        lines((pred0[dat1$drug==ord[jj] & dat1$tray==ii, wid])     ~dat1$adj_conc[dat1$drug==ord[jj] & dat1$tray==ii], lty=1, col=Cols[1])
        lines((pred0[dat1$drug==ord[jj] & dat1$tray==ii, wid-1])   ~dat1$adj_conc[dat1$drug==ord[jj] & dat1$tray==ii], lty=2, col=Cols[2])
        grid(col= "black", lwd = 1)
        
        leg = c("tray estimate", "population estimate" )
        if (is.na(match("legend", names(arglist))) == FALSE) {
          leg = arglist$legend
          if (length(leg) == 2) 
            legend(x="bottomleft", leg, lty=c(1,2), col=c(Cols[1],Cols[2]))
          }
        if (is.na(match("legend", names(arglist))) == TRUE)
          legend(x="bottomleft", leg, lty=c(1,2), col=c(Cols[1],Cols[2]))   
        
        }
      }
    
    # make qqnorm plots for each drug  
    for (jj in 1:length(ord)) {
      tit = paste("NLME: ",ord[jj],", Cell=", as.character(dat1$cell[1]), "\n Model= ",  best, sep="")
      qqnorm(residu[dat1$drug==ord[jj]], main= tit)
      qqline(residu[dat1$drug==ord[jj]] )
      grid(col= "black", lwd = 1)
      }
      
    }
  return()
  }

