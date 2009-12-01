`NlsPlotGetData` <-
function(mixlowData, nlsData, trays) 
  {
  ## collect NLS data for plotting (by tray and by drug)

  data0 = mixlowData$concentrationResponse
  drugRatios = mixlowData$drugRatios
  yAdjust = mixlowData$plottingData
  nls.graphing = nlsData$nlsGraphing
  nls.estimates = nlsData$nlsEstimates
  
  graphDataTray = list(0) 
  graphDataDrug = list(0)
  trayCnt = 0
  usedDrugs = numeric(0)
  
   
  # loop through each tray -------------------------
  for (tr in trays)
    {
    trayCnt = trayCnt + 1
    drug = unique(as.vector(drugRatios$drug[drugRatios$tray==tr]))
    cell = unique(as.vector(drugRatios$cell[drugRatios$tray==tr]))
    Units = unique(as.vector(drugRatios$Units[drugRatios$tray==tr]))
    usedDrugs = c(usedDrugs, drug)
    
    # make empty plot for graphing points
    dat2 = data0[data0$tray==tr & data0$label == "rx",]
    title = paste("NLS: Tray= ",tr, ", \nDrug= ", drug, ", \nCell Line= ",cell, sep="")
    Cols = c("blue","red")
    
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
    
    pred = nls.graphing[[trayCnt]]$pred
    xx = nls.graphing[[trayCnt]]$xx
    
    
    graphDataTray[[trayCnt]] = list(tray= tr, drug=drug, cell=cell, Units= Units, title= title,
    adjResp= dat2$adj_resp[dat2$tray==tr], adjConc= dat2$adj_conc[dat2$tray==tr], 
    ylim= c(0,max(dat2$adj_resp[dat2$tray==tr])), 
    linesY= as.vector(tmp$y), linesX= as.vector(tmp$x), 
    blankPointsY= yA, blankPointsX= xA,
    blankLinesY= predA, blankLinesX= xxA,
    predLinesY= pred, predLinesX= xx)

    }

  # graph for each drug ------------------------------------------------
  udrug = unique(usedDrugs)
  tmp4 = 0
  tmp5 = 0
  drugCnt = 0
  # loop through each unique drug
  for (ud in udrug)
    {
    drugCnt = drugCnt + 1
    # get tray list
    drugTrays = as.vector(unique(drugRatios$tray[drugRatios$drug==ud]))
    Units = as.vector(unique(drugRatios$Units[drugRatios$drug==ud]))
    # get all predicted values to obtain range
    yall = numeric(0)
    for (tr in drugTrays)
      {
      u = exp(as.numeric(as.vector(nls.estimates$u[nls.estimates$tray==tr])))
      for (j in seq(1,length(nls.graphing)))
        {
        if (nls.graphing[[j]]$tray == tr) {
          tmp4 = nls.graphing[[j]]
          break
          }
        }
      yall = c(yall,tmp4$pred/u)
      }
    
    y = nls.graphing[[trayCnt]]$y
    x2 = nls.graphing[[trayCnt]]$x2
    
    # make emtpy plot
    tit = paste("NLS Estimated Scaled Response: \nDrug= ", ud, ", \nCell Line= ",cell, sep="")
    
    cnt = 0
    # loop through each tray
    lineList = list(0)
    for (tr in drugTrays)
      {
      u = exp(as.numeric(as.vector(nls.estimates$u[nls.estimates$tray==tr])))
      for (j in seq(1,length(nls.graphing))){
        if (nls.graphing[[j]]$tray == tr) {
          tmp5 = nls.graphing[[j]]
          break
          }
        }
      cnt = cnt + 1
      lineList[[cnt]] = list(linesY= tmp5$pred/u, linesX= tmp5$xx)
      }
    leg = paste("Tray: ", as.vector(drugTrays),sep="")
    
    graphDataDrug[[drugCnt]] = list(trays= drugTrays, Units= Units, title= tit,drug= ud,
    y= y, x= x2, ylim= c(0,max(yall)), xlim= c(min(xx),max(xx)),
    lineList=lineList, legend= leg)
    
    }
  graphData = list(trays=graphDataTray, drugs= graphDataDrug) 
  return(graphData)
  }

