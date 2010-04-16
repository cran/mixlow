`summary.nlsData` <-
function(object, ...) {
  ## summarize results from the NLS analysis
  
  if (!inherits(object, "nlsData")) 
    stop("use only with \"nlsData\" objects")

  arglist = list(...) 
    
  print(object, verbose=FALSE)
  }

