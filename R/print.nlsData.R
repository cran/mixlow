`print.nlsData` <-
function(x, ...) {
  ## print results from the NLS analysis
  
  if (!inherits(x, "nlsData")) 
    stop("use only with \"nlsData\" objects")

  arglist = list(...) 
  verbose = arglist$verbose
  
  if (is.null(arglist$verbose)) verbose = TRUE
  
  nlsEstimates = x$nlsEstimates
  nlsGraphing = x$nlsGraphing
  

  
  writeLines("\n ====================== NLS Parameter Estimates ======================\n") 
  print(nlsEstimates)
  
  
  if (verbose == TRUE){
    # loop through each nlsGraphing set -------------------
    writeLines("\n\n ====================== Graphing Data for NLS ======================\n") 
    for (iii in seq(1,length(nlsGraphing))) {  
      writeLines(paste("\n\n -------- graphingData item ", iii, " --------\n", sep=""))
      writeLines(paste("tray:  ", nlsGraphing[[iii]]$tray, sep=""))
      writeLines(paste("drug:  ", nlsGraphing[[iii]]$drug, sep=""))
      
      tmp = data.frame(adj.conc=nlsGraphing[[iii]]$x2, adj.response=nlsGraphing[[iii]]$y)
      writeLines("\nadjusted concentration-response data:")
      print(tmp)
      
      tmp = data.frame(adj.conc=nlsGraphing[[iii]]$xx, predicted.response=nlsGraphing[[iii]]$pred)
      writeLines("\npredicted concentration-response data:")
      print(tmp)    
      }
    }

  
  }

