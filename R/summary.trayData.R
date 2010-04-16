`summary.trayData` <-
function(object, ...) {
  ## summarize trayData object returned by the \code{readDataFile} function
  
  if (!inherits(object, "trayData")) 
    stop("use only with \"trayData\" objects")

  rawData = object$concentrationResponse
  drugRatios = object$drugRatios

  
  writeLines ("\n\n  ***  Data Summary  ***\n")
  writeLines ( paste("number of wells = ",dim(rawData)[1],sep=""))
  writeLines (paste("number of drugs = ", length(getDrugs(object)),sep=""))
  writeLines (paste("number of trays = ", length(getTrays(object)),sep=""))
  writeLines (paste("number of cell lines = ", length(getCellLines(object)),"\n",sep=""))
  writeLines ("\n  ====================== Drug Ratios ======================  ")
  print(object$drugRatios)
  writeLines(" ")
  
  }

