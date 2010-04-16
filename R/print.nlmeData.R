`print.nlmeData` <-
function(x, ...) {
  ## print results from the NLME analysis
  
  if (!inherits(x, "nlmeData")) 
    stop("use only with \"nlmeData\" objects")
  
  arglist = list(...) 
  verbose = arglist$verbose
  
  if (is.null(arglist$verbose)) verbose = TRUE

  
  nlmeResults = x$nlmeResults
  nlmeGraph = x$nlmeGraph    
  nlmeModels = x$nlmeModels
  
  # loop through each nlme model set -------------------
  writeLines("\n ====================== NLME Models ======================\n") 
  for (iii in seq(1,length(nlmeModels))) {  
    writeLines(paste("\n\n -------- NLME Model item ", iii, " --------\n", sep=""))
    if (verbose == TRUE){
      print(summary(nlmeModels[[iii]]))
      }
    else {
     print(nlmeModels[[iii]]) 
     }
  }
  
  if (verbose == TRUE){
    # loop through each nlme results set -------------------
    writeLines("\n ====================== NLME Results ======================\n") 
    for (iii in seq(1,length(nlmeResults))) {  
      writeLines(paste("\n\n -------- NLME Results item ", iii, " --------\n", sep=""))
      writeLines(paste("model name:  ", nlmeResults[[iii]]$nam, sep=""))
      writeLines(paste("method:  ", nlmeResults[[iii]]$method, sep=""))
      writeLines("\ndrugs:  ")
      print(unlist(strsplit(nlmeResults[[iii]]$drug, split="\\|")))
      writeLines(paste("\nresult set number:  ", nlmeResults[[iii]]$setNum, sep=""))
      writeLines(paste("cell line:  ", nlmeResults[[iii]]$cell, sep=""))
      writeLines(paste("likelihood:  ", nlmeResults[[iii]]$lik, sep=""))
      writeLines(paste("BIC:  ", nlmeResults[[iii]]$bic, sep=""))
      writeLines(paste("sigma:  ", nlmeResults[[iii]]$sig, sep=""))
      
      writeLines("\nmodel structure:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,11) == "modelstruct") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }
        
      writeLines("\nrc values:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,2) == "rc") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\nr-squared values:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,2) == "r2") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\ncoefficient gamma:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,13) == "coeff.fixed.g") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\ncoefficient psi:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,13) == "coeff.fixed.p") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\ncoefficient mu:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,13) == "coeff.fixed.u") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\ncoefficient lambda:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,18) == "coeff.fixed.lambda") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\ncoefficient random:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,12) == "coeff.random") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\nstandard error gamma:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,4) == "se.g") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\nstandard error psi:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,4) == "se.p") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\nstandard error mu:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,4) == "se.u") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\nstandard error lambda:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,9) == "se.lambda") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\ndegrees of freedom:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,3) == "df.") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }

      writeLines("\ncovariance:")
      for (j in names(nlmeResults[[iii]])) {
        if (substr(j,1,4) == "covv") writeLines (paste("  ", j, ":  ", eval(parse(text=(paste("nlmeResults[[",iii,"]]$",j,sep="")))), sep=""))
        }    
      }  



    # loop through each nlme graph set -------------------
    writeLines("\n ====================== NLME Graph Data ======================\n") 
    for (iii in seq(1,length(nlmeGraph))) {
      
      writeLines(paste("\n\n -------- NLME Graph Data item ", iii, " --------\n", sep=""))
      
      writeLines("\ndrugs:")
      print (nlmeGraph[[iii]]$ord)

      writeLines("\nname of best model:")
      print (nlmeGraph[[iii]]$best)
      
      writeLines("\ndata:")
      print (nlmeGraph[[iii]]$dat1)
        
      writeLines("\npredictions:")
      print (nlmeGraph[[iii]]$pred0)
      
      writeLines("\nresiduals:")
      print (nlmeGraph[[iii]]$residu)
      }
    }      


  }

