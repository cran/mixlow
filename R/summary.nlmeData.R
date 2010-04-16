`summary.nlmeData` <-
function(object, ...) {
  ## summarize results from the NLME analysis
  
  if (!inherits(object, "nlmeData")) 
    stop("use only with \"nlmeData\" objects")
  
  arglist = list(...) 
    
  print(object, verbose=FALSE)
  }

