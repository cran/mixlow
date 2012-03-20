`getCellLines` <-
function(data, drugs=NULL, trays=NULL) {
  ## scans dataframe and returns vector of cell lines, subset for given drugs and trays if desired

  if ((!inherits(data, "trayData")) && (!inherits(data, "mixlowData"))) 
    stop("use only with \"trayData\" or \"mixlowData\" objects")
      
  drugRatios = data$drugRatios

  if ( is.null(drugs) & (is.null(trays)==FALSE) ) {
    cellLines = as.vector(drugRatios$cell[(drugRatios$tray %in% trays)])
    return (sort(unique(cellLines)))
    }
  
  if ( (is.null(drugs)==FALSE) & is.null(trays) ) {
    cellLines = as.vector(drugRatios$cell[(drugRatios$drug %in% drugs)])
    return (sort(unique(cellLines)))
    }

  if ( (is.null(drugs)==FALSE) & (is.null(trays)==FALSE) ) {
    cellLines = as.vector(drugRatios$cell[(drugRatios$drug %in% drugs) & (drugRatios$tray %in% trays)])
    return (sort(unique(cellLines)))
    }
  
  if ( is.null(drugs) & is.null(trays) ) {
    cellLines = as.vector(drugRatios$cell)
    return (sort(unique(cellLines)))
    }  

  }
