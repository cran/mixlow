`getTrays` <-
function(data, drugs=NULL, cellLines=NULL) {
  ## scans dataframe and returns vector of trays, subset for given drugs and cell lines if desired
  
  if ((!inherits(data, "trayData")) && (!inherits(data, "mixlowData"))) 
    stop("use only with \"trayData\" or \"mixlowData\" objects")
  
  drugRatios = data$drugRatios
  
  if ( is.null(drugs) & (is.null(cellLines)==FALSE) ) {
    trays = as.vector(drugRatios$tray[(drugRatios$cell %in% cellLines)])
    return (trays)
    }
  
  if ( (is.null(drugs)==FALSE) & is.null(cellLines) ) {
    trays = as.vector(drugRatios$tray[(drugRatios$drug %in% drugs)])
    return (trays)
    }

  if ( (is.null(drugs)==FALSE) & (is.null(cellLines)==FALSE) ) {
    trays = as.vector(drugRatios$tray[(drugRatios$drug %in% drugs) & (drugRatios$cell %in% cellLines)])
    return (trays)
    }
  
  if ( is.null(drugs) & is.null(cellLines) ) {
    trays = as.vector(drugRatios$tray)
    return (trays)
    }  
  }
