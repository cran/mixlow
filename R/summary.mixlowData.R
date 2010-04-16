`summary.mixlowData` <-
function(object, ...) {
  ## summarize results from read data
  
  if (!inherits(object, "mixlowData")) 
    stop("use only with \"mixlowData\" objects")

  arglist = list(...) 
    
  print(object, verbose=FALSE)
  }

