`getNlsParameterDefaults` <- 
function(trays) {
  ## creates a data frame to hold default parameters for the NLS function 
  ## std= standard deviation of parameter for use in generating random values
  ## lambda= fraction of max for nonzero asymptote.  
  ## Null means let data determine best values 
  
  
  numberOfTrays <- length(trays)
  
  parameterDefaults <- data.frame(tray=rep("trayNum",numberOfTrays),
                       param.g=rep(NaN, numberOfTrays),
                       param.g.std = rep(NaN, numberOfTrays),
                       param.p=rep(NaN, numberOfTrays),
                       param.p.std = rep(NaN, numberOfTrays),
                       param.lambda=rep(NaN, numberOfTrays))
  parameterDefaults$tray = trays
  row.names(parameterDefaults) = parameterDefaults$tray
  return (parameterDefaults)
  }



