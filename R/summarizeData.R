`summarizeData` <-
function(object) {
  ## print summary of data frame

  if ((!inherits(object, "trayData")) && (!inherits(object, "mixlowData"))) 
    stop("use only with \"trayData\" or \"mixlowData\" objects")
  
  rawData = object$concentrationResponse
  drugRatios = object$drugRatios

  
  writeLines ("\n\n  ***  Data Summary  ***\n")
  writeLines ( paste("number of wells = ",dim(rawData)[1],sep=""))
  writeLines (paste("number of drugs = ", length(getDrugs(object)),sep=""))
  writeLines (paste("number of trays = ", length(getTrays(object)),sep=""))
  writeLines (paste("number of cell lines = ", length(getCellLines(object)),"\n",sep=""))
  writeLines ("\n            Table of drug ratios")
  print(object$drugRatios)
  writeLines(" ")
  
  }
  
