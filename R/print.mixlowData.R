`print.mixlowData` <-
function(x, ...) {
  ## print results from read data
  
  if (!inherits(x, "mixlowData")) 
    stop("use only with \"mixlowData\" objects")

  arglist = list(...) 
  
  concentrationResponse = x$concentrationResponse
  drugRatios = x$drugRatios
  plottingData = x$plottingData
  
  writeLines("\n\n ====================== Drug Ratios ======================\n") 
  print (drugRatios)
  
  writeLines("\n\n ====================== Concentration-Response Data ======================\n") 
  print (concentrationResponse)
 
  
  # loop through each plottingData set -------------------
  writeLines("\n\n ====================== Plotting Data for Blank Wells ======================\n") 
  for (iii in seq(1,length(plottingData))) {  
    writeLines(paste("\n\n -------- plottingData item ", iii, " --------\n", sep=""))
    writeLines(paste("tray:  ", plottingData[[iii]]$tray, sep=""))
    
    tmp = data.frame(blank.conc=plottingData[[iii]]$x, blank.response=plottingData[[iii]]$y)
    writeLines("\nactual blanks data:")
    print(tmp)
    
    tmp = data.frame(blank.conc=plottingData[[iii]]$xx, blank.response=plottingData[[iii]]$pred)
    writeLines("\npredicted blanks data:")
    print(tmp)    
    }

  
  }

