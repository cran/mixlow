`print.loeweData` <-
function(x, ...) {
  ## print results from Loewe analysis
  
  if (!inherits(x, "loeweData")) 
    stop("use only with \"loeweData\" objects")

  arglist = list(...) 
  verbose = arglist$verbose
  
  if (is.null(arglist$verbose)) verbose = TRUE

  loeweData = x

  if (verbose == TRUE){
    writeLines("\ndrugs:  ")
    print(loeweData$drugs)

    writeLines(paste("\nmixture:  ", loeweData$mix, sep=""))

    writeLines("\ngamma:  ")
    print(loeweData$g0)

    writeLines("\npsi:  ")
    print(loeweData$p0)  
    
    writeLines("\nlambda:  ")
    print(loeweData$lam0)
    
    writeLines(paste("\nmu:  ", loeweData$uu, sep=""))
    
    writeLines(paste("\nunits:  ", loeweData$Units, sep=""))
    }
    
  writeLines("\nLog Loewe index and confidence intervals:  ")
  if (verbose == TRUE) print(loeweData$ciL)  
  if (verbose == FALSE) {
    n = length(loeweData$ciL$Fraction.Affected)
    print(loeweData$ciL[seq(1,n,n/20),])
    }
  
  if (verbose == TRUE){
    writeLines("\ncovariance:  ")
    print(loeweData$covv)

    writeLines("\ndata:  ")
    print(loeweData$dat0)

    writeLines("\nscore interval for antagonism/synergism summary:  ")
    print(loeweData$score.interval)

    writeLines(paste("\nsummary synergism score:  ", loeweData$syner.total, sep=""))

    writeLines(paste("\nsummary antagonism score:  ", loeweData$antag.total, sep=""))
    
    writeLines("\nsynergism score:  ")
    print(loeweData$syner)  
    
    writeLines("\nantagonism score:  ")
    print(loeweData$antag)    
  }

  
  }

