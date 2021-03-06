`prepareData` <-
function(trayData, drugs= NULL, trays=NULL, cellLines=NULL, degree=3) 
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
    if (abs(sum(ratios)-1) > .000001){
      writeLines(paste("\nsum= ", sum(ratios), sep=""))  
      stop("Ratios (by row in drugRatios) should sum to one.")
      }
    }
    
  # check that only one cell line is selected
  if (length(cellLines) > 1) 
    stop("Only one cell line should be selected.")
  
  
  crData$adj_conc = crData$conc
  crData$adj_resp = crData$resp

  # adjust conc values for adj_conc column (to allow printing of zero conc values)
  # concentrations of zero are adjusted to 1/1000 of the lowest used conc to allow graphing
  for (tr in trays)
    {
    # get controls
    tmp = ""
    tmp = crData[crData$tray==tr & crData$conc > 0, ]
    minn = min(tmp$conc)
    # fix dat concentrations
    crData$adj_conc[crData$tray==tr & crData$conc == 0] = minn/1000.
    }

  # adjust for blanks
  
  # below, responses are adjusted to account for bbt or bbc data (blanks)
  plottingData = vector("list", length(trays) )
  trayCnt = 0
  for (tr in trays)
    {
    trayCnt = trayCnt + 1
    # check that both bbt and bbc are not used in the same tray
    tmp1 = as.vector(crData$label[crData$tray==tr & crData$label != "rx"])

    if (("bbt" %in% tmp1) & ("bbc" %in% tmp1))
      {stop("Cannot use both bbc and bbt in data for a given tray")}
    if (unique(tmp1) == "bbt")
      {
      # adjust conc for bbt blanks: -----------------------------------------
      tmp = as.list(crData[crData$tray==tr & crData$label == "bbt", c("adj_conc", "resp")])
      blank = mean(as.numeric(as.vector(tmp$resp)))
      crData$adj_resp[crData$tray==tr & crData$label == "rx"] = crData$resp[crData$tray==tr & crData$label == "rx"] - blank
      
      # save  for plotting
      tmp2 = crData[crData$tray==tr & crData$label == "rx", ]
      x2 = tmp2$adj_conc
      y = tmp2$resp
      x = rep(0,length(y))
      xx = exp(seq(log(min(x2)),log(max(x2)),length.out = 500))
      pred = rep(blank,length(xx))
      plottingData[[trayCnt]] = list(tray=tr,x=x,y=y, xx=xx,pred=pred)
      }
    if (unique(tmp1) == "bbc")
      {
      # adjust conc for bbc blanks: -----------------------------------------
      tmp = crData[crData$tray==tr & crData$label == "bbc", c("adj_conc", "resp")]
      tmp = tapply(tmp$resp, tmp$adj_conc, mean)       
      
      # fit response~concentration in "blanks" wells using polynomial
      y = as.numeric(tmp)
      x = as.numeric(names(tmp))
      
      if (length(x) == 1)
        {stop("Only one summarized concentration-response point available---did you include a proper concentration for each bbc well?")}
      
      if (degree == 4) model = nls(y~ a + b*x + c*x^2 + d*x^3 + e*x^4, data.frame(y,x), start= c(a=min(y),b=0,c=0,d=0,e=0))
      if (degree == 3) model = nls(y~ a + b*x + c*x^2 + d*x^3 , data.frame(y,x), start= c(a=min(y),b=0,c=0,d=0))
      if (degree == 2) model = nls(y~ a + b*x + c*x^2, data.frame(y,x), start= c(a=min(y),b=0,c=0))
      if (degree == 1) model = nls(y~ a + b*x, data.frame(y,x), start= c(a=min(y),b=0))
      if (model$convInfo$isConv == FALSE)
        {stop("nls model did not converge for bbc")}
      sm = summary(model)

      param = sm$parameters[,1]

      # save  for plotting
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
  
  

