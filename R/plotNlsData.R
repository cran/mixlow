`plotNlsData` <-
function(nlsData, mixlowData, trays = getTrays(mixlowData), 
    ask = prod(par("mfcol")) < length(trays) && dev.interactive(), showBlanks= FALSE) {
  ## graphs results from the NLS analysis

  if (!inherits(nlsData, "nlsData")) 
    stop("use only with \"nlsData\" objects")

  if (ask) {
    opar <- par(ask=TRUE)
    on.exit(par(opar))
    }
  
  # collect data for graphing
  graphData = NlsPlotGetData (mixlowData, nlsData, trays)
  
  
  # loop through each tray that is wanted -------------------------
  trayCnt = 0
  useDrugs = NULL
  for (ii in seq(1,length(graphData$trays))) {  
    trayData = graphData$trays[[ii]]
    tr = trayData$tray
    drug = trayData$drug
    Units = trayData$Units
    # collect drugs that are in desired trays
    useDrugs = c(useDrugs, drug)
      
    trayCnt = trayCnt + 1
    
    # make empty plot for graphing points
    tit = trayData$title
    yLab = "Response"
    xLab = paste("Concentration, w/adj for zero, ", Units,sep="")
    Cols = c("blue","red")
    
    # plot data points
    plot(trayData$adjResp ~ trayData$adjConc,
          log="x",type="n", main= trayData$title,
          ylab= yLab, xlab= xLab,
          ylim= trayData$ylim)

    points(trayData$adjResp ~ trayData$adjConc, pch=1, col=Cols[1])
    grid(col= "black", lwd = 1)
    
    # graph line through means
    lines(trayData$linesY ~ trayData$linesX, col=Cols[1], lty = 2)
    
    # graph adjustments to response

    if (showBlanks){
      points(trayData$blankPointsY ~ trayData$blankPointsX, col=Cols[2])
      lines(trayData$blankLinesY ~ trayData$blankLinesX, col=Cols[2])
      }
    
    # graph predictions
    lines(trayData$predLinesY ~ trayData$predLinesX, col=Cols[1], lwd = 1.5 )

    leg1 = c("predicted", "means", "adjustments")
    leg2 = c("predicted", "means")
    if (showBlanks==TRUE) legend(x="topright", legend=leg1, lty= c(1,2,1),
      col= c(Cols[1], Cols[1], Cols[2]), lwd= c(1.5,2,1) )
    if (showBlanks==FALSE) legend(x="topright", legend=leg2, lty= c(1,2),
      col= c(Cols[1], Cols[1]), lwd= c(1.5,2) )
    }

  # graphs curves for each drug ------------------------------------------------
  udrug = as.vector(unique(useDrugs))
  
  # loop through each unique drug
  for (ii in seq(1,length(graphData$drugs))) { 

    drugData = graphData$drugs[[ii]] 
    

    
    # get tray list
    Units = drugData$Units
    drugTrays = drugData$trays
    # make emtpy plot
    yLab = "Response"
    xLab = paste("Concentration, w/adj for zero, ", Units,sep="")
    cols = rep(c("black", "blue", "red", "purple", "green", "cyan"),4,each=1)
    #if (formatFlags$color == 0) cols = rep("black",24)
    typs = rep(1:6,4)
    plot(drugData$y ~ drugData$x,
          log="x", type="n", main= drugData$title,
          ylab= yLab, xlab= xLab, ylim= drugData$ylim)
    
    cnt = 0
    # loop through each tray
    for (tr in drugTrays) {
      cnt = cnt + 1
      lines(drugData$lineList[[cnt]]$linesY ~ drugData$lineList[[cnt]]$linesX, type="l", col=cols[cnt], lty=typs[cnt], lwd= 2)
      }
    
    leg = drugData$legend
    legend(x="topright", legend=leg, lty=typs[1:length(drugTrays)], col= cols[1:length(drugTrays)])
    grid(col= "black", lwd = 1)
    
    }
  return(graphData)
  }

    
    
    
