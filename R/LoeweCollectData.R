`LoeweCollectData` <-
function(nlmeResults, drugRatios) {
  # collect nlme data  -----------------------------------------------------
  # need to determine which mixture goes with which drugs
  # find all mixtures
  
  if (length(nlmeResults) == 1) {
    mixtureList = list(analysis="multiple")
    drugs = nlmeResults[[1]]$drug
    drugs = strsplit(drugs,"\\|")[[1]]
    # mixture is last in list
    mixtureList$mixture = drugs[length(drugs)]
    mixtureList$drugs = drugs[1:length(drugs)-1]
    mixtureList$locations = 1
    }

  if (length(nlmeResults) > 1) {
    mixtureList = list(analysis = "single")
    
    # get drug list
    drugs = numeric(0)
    for (i in seq(1, length(nlmeResults)-1)) 
      drugs = c(drugs, nlmeResults[[i]]$drug)
    mixtureList$drugs = drugs 
    # mixture is last 
    mixtureList$mixture = nlmeResults[[i+1]]$drug
    mixtureList$locations = seq(1,length(nlmeResults))
    }

  # check that ratios for single drugs is unique for mixtures and 1 for drug itself
  frac = numeric()
  for (dr in mixtureList$drugs) {
    ratios = unique(as.vector(drugRatios[,dr][drugRatios$drug==mixtureList$mixture]))
    if (length(ratios) > 1)
      stop("Ratios across trays for any drug in a mixture cannot have more than one value")
    ratios = as.vector(drugRatios[,dr][drugRatios$drug==dr])
    if (!all(ratios == 1))
      stop("Ratios for any single drug must be equal to one")
    }  


  # get fractions
  fractions = numeric(0)
  for (dr in drugs) {
    frac = as.numeric(unique(as.vector(drugRatios[drugRatios$drug==mixtureList$mixture,][names(drugRatios) == dr])))
    fractions = c(fractions, frac)
    }
  mixtureList$fractions = fractions
  
  return(mixtureList)

  }

