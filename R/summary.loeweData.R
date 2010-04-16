`summary.loeweData` <-
function(object, ...) {
  ## summarize results from Loewe analysis
  
  if (!inherits(object, "loeweData")) 
    stop("use only with \"loeweData\" objects")

  arglist = list(...) 
    
  print(object, verbose=FALSE)
  }

