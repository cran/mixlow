`mixlowPlotGetData` <-
function(data, trays) 
  ## collect adjusted response data for plotting
  {
  data0 = data$concentrationResponse
  yAdjust= data$plottingData
  drugRatios = data$drugRatios
  
  trayCnt = 0
  graphingData = list(0) 
  
  # loop through each tray -------------------------
  for (tr in trays)
    {
    trayCnt = trayCnt + 1
    drug = unique(as.vector(drugRatios$drug[drugRatios$tray==tr]))
    cell = unique(as.vector(drugRatios$cell[drugRatios$tray==tr]))
    Units = unique(as.vector(drugRatios$Units[drugRatios$tray==tr]))
    
    # make empty plot for graphing points
    dat2 = data0[data0$tray==tr & data0$label == "rx",]
    tit = paste("Raw Data: Tray= ",tr, ",\nDrug= ", drug, ",\nCell Line= ",cell, sep="")
     
    # graph line through means
    tmp = as.data.frame(cbind(x= dat2$adj_conc[dat2$tray==tr], y = dat2$adj_resp[dat2$tray==tr]))
    tmp = aggregate(tmp$y, list(X=tmp$x) ,mean)
    names(tmp) <- c("x","y")
    
    
    # graph adjustments to response
    for (j in seq(1,length(yAdjust))){
      if (yAdjust[[j]]$tray == tr) {
        Adj = yAdjust[[j]]
        break
        }
      }

    xA = Adj$x
    yA = Adj$y
    xxA = Adj$xx
    predA = Adj$pred
    
    graphingData[[trayCnt]] = list(drug=drug, cell=cell, Units= Units, title= tit, 
      adjResp= dat2$adj_resp[dat2$tray %in% trays], adjConc= dat2$adj_conc[dat2$tray %in%trays], 
      ylim= c(0, max(data0$adj_resp[dat2$tray==tr])), 
      meanLineY= as.vector(tmp$y), meanLineX= as.vector(tmp$x),
      blankY= yA, blankX= xA, blankLineY= predA, blankLineX= xxA)
    
    }
    
  return (graphingData= graphingData)
  }

