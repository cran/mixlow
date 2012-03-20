`getDrugs` <-
function(data, trays=NULL, cellLines=NULL) {
  ## scans dataframe and returns vector of drugs, subset for given trays and cell lines if desired
  
  if ((!inherits(data, "trayData")) && (!inherits(data, "mixlowData"))) 
    stop("use only with \"trayData\" or \"mixlowData\" objects")
  
  drugRatios = data$drugRatios
  
  if ( is.null(trays) & (is.null(cellLines)==FALSE) ) {
    drugs = as.vector(drugRatios$drug[(drugRatios$cell %in% trays)])
    return (sort(unique(drugs)))
    }
  
  if ( (is.null(trays)==FALSE) & is.null(cellLines) ) {
    drugs = as.vector(drugRatios$drug[(drugRatios$tray %in% trays)])
    return (sort(unique(drugs)))
    }

  if ( (is.null(trays)==FALSE) & (is.null(cellLines)==FALSE) ) {
    drugs = as.vector(drugRatios$drug[(drugRatios$tray %in% trays) & (drugRatios$cell %in% cellLines)])
    return (sort(unique(drugs)))
    }
  
  if ( is.null(trays) & is.null(cellLines) ) {
    drugs = as.vector(drugRatios$drug)
    return (sort(unique(drugs)))
    }  

  }
