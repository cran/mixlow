`plot.mixlowData` <-
function(x, ...) {
  ## plots the adjusted data by tray
  
  if (!inherits(x, "mixlowData")) 
    stop("Use only with \"mixlowData\" objects.  Run prepareData() first.")
  
  arglist = list(...)  
  trays = arglist$trays
  ask = arglist$ask
  showBlanks = arglist$showBlanks
  
  if (is.null(arglist$trays)) trays = getTrays(x)
  if (is.null(arglist$ask)) ask = prod(par("mfcol")) < length(trays) && dev.interactive()
  if (is.null(arglist$showBlanks)) showBlanks= FALSE
  
  opar <- par(ask=ask)
  on.exit(par(opar))
  
  mixlowData = x
  
  # collect data for plotting
  graphingData = mixlowPlotGetData(mixlowData,trays)
    
  # loop through each tray that is wanted -------------------------
  for (ii in seq(1,length(trays))) {
    tr = trays[ii]
    graphdata = graphingData[[ii]]
    Cols = c("blue","red")
    yLab = "Response"
    xLab = paste("Concentration, w/adj for zero, ", graphdata$Units,sep="")

    plot(graphdata$adjResp ~ graphdata$adjConc, log="x", type="n", main= graphdata$title, ylab= yLab, xlab= xLab,
          ylim= graphdata$ylim)

    points(graphdata$adjResp ~ graphdata$adjConc, pch=19, col=Cols[1])
    grid(col= "black", lwd = 1)
    
    # graph line through means
    lines(graphdata$meanLineY ~ graphdata$meanLineX, col=Cols[1], lty = 2)
    
    # graph adjustments to response
    if (showBlanks == TRUE){
      points(graphdata$blankY ~ graphdata$blankX, col=Cols[2])
      lines(graphdata$blankLineY ~ graphdata$blankLineX, col=Cols[2])
      }
    
    leg1 = c("means", "adjustments")
    leg2 = c("means")
    if (showBlanks==TRUE) legend(x="topright", legend=leg1, lty= c(2,1),
      col= c(Cols[1], Cols[2]), lwd= c(1.5,2,1) )
    if (showBlanks==FALSE) legend(x="topright", legend=leg2, lty= c(2),
      col= c(Cols[1]), lwd= c(1.5,2) )
    }

  
  }
  
  
  
  

      

