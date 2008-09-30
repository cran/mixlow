`prepareData` <-
function(trayData, drugs= NULL, trays=NULL, cellLines=NULL) 
  {
  ## adjusts all responses to account for "blanks", adds two new columns to data frame
  ## returns adjusted data along with data frame for plotting "blanks" values
  ## new columns are: tray, row, col, label, conc, resp, adj_conc(no zero, for plotting), adj_resp(adj for blanks)

  if ( !inherits(trayData, "trayData")) 
    stop("use only with \"trayData\" objects")

  if (is.null(drugs)) drugs = getDrugs(trayData)
  if (is.null(trays)) trays = getTrays(trayData) 
  if (is.null(cellLines)) cellLines = getCellLines(trayData)
  
  crData = trayData$concentrationResponse
  crHead = names(crData)
  dRatios = trayData$drugRatios
  dHead = names(dRatios)
  data0 = merge(crData, dRatios, by= "tray")   
  
  # select subset
  crData = data0[((data0$drug %in% drugs) & (data0$tray %in% trays) & (data0$cell %in% cellLines)),crHead]
  dRatios = unique(data0[((data0$drug %in% drugs) & (data0$tray %in% trays) & (data0$cell %in% cellLines)),dHead])
  row.names(dRatios) = NULL
  
  # adjust drugs, cellLines, and trays vectors
  drugs = unique(as.vector(dRatios$drug))
  trays = unique(as.vector(dRatios$tray))
  
  # check that all ratios sum to one
  for (i in seq(1,dim(dRatios)[1]))  {
    ratios = as.vector(dRatios[i,7:dim(dRatios)[2]])
    if (!sum(ratios)==1)
      stop("Ratios (by row in drugRatios) should sum to one.")
    }
    
  # check that only one cell line is selected
  if (length(cellLines) > 1) 
    stop("Only one cell line should be selected.")
  
  
  crData$adj_conc = crData$conc
  crData$adj_resp = crData$resp

  # adjust conc values for adj_conc column (to allow printing of zero conc values)
  for (tr in trays)
    {
    # get controls
    tmp = ""
    tmp = crData[crData$tray==tr & crData$conc > 0, ]
    minn = min(tmp$conc)
    # fix dat concentrations
    crData$adj_conc[crData$tray==tr & crData$conc == 0] = minn/1000.
    }

  # adjust for blanks, assumes all rows in column receive the same conc
  # below, concentrations of zero are adjusted to 1/1000 of the lowest used conc to allow graphing
  # below, responses are adjusted to account for bbt or bbc data (blanks)
  plottingData = vector("list", length(trays) )
  trayCnt = 0
  for (tr in trays)
    {
    trayCnt = trayCnt + 1
    # check that both bbt and bbc are not used in the same tray
    tmp1 = as.vector(crData$label[crData$tray==tr & crData$label != "rx"])

    if (("bbt" %in% tmp) & ("bbc" %in% tmp))
      {stop("Cannot use both bbc and bbt in data for a given tray")}
    if (unique(tmp1) == "bbt")
      {
      # adjust conc for bbt blanks: -----------------------------------------
      tmp = as.list(crData[crData$tray==tr & crData$label == "bbt", c("col", "adj_conc", "resp")])
      blank = mean(tmp$resp)
      crData$adj_resp[crData$tray==tr & crData$label == "rx"] = crData$resp[crData$tray==tr & crData$label == "rx"] - blank
      
      col_lookup = unique(crData[crData$tray==tr & crData$label == "rx", c("col","adj_conc")])
      tmp = merge(tmp, col_lookup, by= "col")

      # save  for plotting
      tmp2 = crData[crData$tray==tr & crData$label == "rx", ]
      x2 = tmp2$adj_conc
      y = tmp$resp
      x = rep(0,length(y))
      xx = exp(seq(log(min(x2)),log(max(x2)),length.out = 500))
      pred = rep(blank,length(xx))
      plottingData[[trayCnt]] = list(tray=tr,x=x,y=y, xx=xx,pred=pred)
      }
    if (unique(tmp1) == "bbc")
      {
      # adjust conc for bbc blanks: -----------------------------------------
      tmp = crData[crData$tray==tr & crData$label == "bbc", c("col", "adj_conc", "resp")]

      col_lookup = unique(crData[crData$tray==tr & crData$label == "rx", c("col","adj_conc")])
      tmp = merge(tmp,col_lookup, by= "col")
      
      # fit response~concentration in "blanks" wells using polynomial
      y = tmp$resp
      x = tmp$adj_conc.y
      model = nls(y~ a + b*x + c*x^2 + d*x^3 + e*x^4, data.frame(y,x), start= c(a=min(y),b=0,c=0,d=0,e=0))
      if (model$convInfo$isConv == FALSE)
        {stop("nls model did not converge for bbc")}
      sm = summary(model)

      param = sm$parameters[,1]

      tmp2 = crData[crData$tray==tr & crData$label == "rx", ]
      x2 = tmp2$adj_conc
      xx = exp(seq(log(min(x2)),log(max(x2)),length.out = 500))
      pred = predict(model,list(x=xx))
      plottingData[[trayCnt]] = list(tray=tr,x=x,y=y, xx=xx,pred=pred)

      # adjust response
      crData$adj_resp[crData$tray==tr & crData$label == "rx"] =  crData$adj_resp[crData$tray==tr & crData$label == "rx"] -
        predict(model, list(x= crData$adj_conc[crData$tray==tr & crData$label == "rx"])) 
      }
    }
  
  returnList = list(concentrationResponse=crData, drugRatios=dRatios, plottingData=plottingData)
  class(returnList) <- "mixlowData"
  return (returnList)
  }
  
  

