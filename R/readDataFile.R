`readDataFile` <-
function(filename, excludeWells=NULL) {
  ## reads a formatted data file 
  
  ## check that excludeWells contains valid data
  if (is.null(excludeWells)==FALSE) {
    if (length(excludeWells$row) != length(excludeWells$col)) 
      stop("Vector of excluded rows must be of same length as vector of excluded columns")
    }  
        
  txt = readLines(filename)
  
  ## make lists for drug names and tray labels
  drugList = NULL
  trayList = NULL
  
  if (length(grep("global_notes",txt[1])) ==0) 
    stop("First line of file is not formatted correctly")
  
  for(i in seq(1,length(txt))) {
    if (length(grep("drug_name_short\t", txt[i])) > 0) {
      tmp = strsplit(txt[i],"\t")
      dname = as.character(tmp[[1]][2])
      drugList = as.character(c(drugList, dname))
      }
    if (length(grep("tray_label\t", txt[i])) > 0) {
      tmp = strsplit(txt[i],"\t")
      tname = as.character(tmp[[1]][2])
      trayList = c(trayList, tname)
      }
    }
  drugList = sort(unique(drugList))
  trayList = sort(trayList)

  lenTrays = length(trayList)
    
  ## collect data for drug ratios dataframe and for dose-response (data0) dataframe
  drugRatios = 1
  rows = as.numeric(strsplit(txt[7],"\t")[[1]][2])
  cols = as.numeric(strsplit(txt[8],"\t")[[1]][2])
  Units = strsplit(txt[6],"\t")[[1]][2]
  
  n = lenTrays * rows * cols
  data0 = data.frame(tray=character(n), row=numeric(n), col=numeric(n), label=character(n), conc=numeric(n), resp=numeric(n))

  index = 1
  for (i in seq(1,length(txt))) 
    {
    if (length(grep("tray_label\t", txt[i])) > 0) {
      tray = as.character(strsplit(txt[i+ 0],"\t")[[1]][2])
      drug = as.character(strsplit(txt[i+5],"\t")[[1]][2])
      cell = as.character(strsplit(txt[i+3],"\t")[[1]][2])
      comp = strsplit(txt[i+6],"\t")[[1]]
      comp = comp[grep("=",comp)]
      tmp = list(tray=tray, drug=drug, cell=cell, Units=Units, rows=rows, cols=cols)
      
      ## add zero composition values for each drug
      for (j in drugList) {eval(parse(text= paste("tmp$",j,sep=""," = 0"))) }

      if (is.data.frame(drugRatios)==TRUE) {
        ##add new row in drug ratio data frame
        levels(drugRatios$tray) = c(levels(drugRatios$tray),as.vector(tray))
        levels(drugRatios$drug) = c(levels(drugRatios$drug),as.vector(drug))
        levels(drugRatios$cell) = c(levels(drugRatios$cell),as.vector(cell))
        levels(drugRatios$Units) = c(levels(drugRatios$Units),as.vector(Units))
        drugRatios = rbind(drugRatios,tmp)
        }
      if (is.data.frame(drugRatios)==FALSE) drugRatios = data.frame(tmp)

      for (j in comp) {
        ## adjust composition values from zero for current tray
        tmp2 = strsplit(j,"=")
        k2 = as.numeric(tmp2[[1]][2])
        strr = paste("drugRatios$",tmp2[[1]][1],"[drugRatios$tray=='", tray,"'] = ",tmp2[[1]][2],sep="")
        eval(parse(text=strr))
      }

      ## get "rows" rows of dose-response values
      conc = txt[(i+14):(i+14+rows-1)]
      label = txt[(i+14+rows):(i+14+rows+rows-1)]
      resp = txt[(i+14+rows+rows):(i+14+rows+rows+rows-1)]
      conc = strsplit(conc,"\t")
      label = strsplit(label,"\t")
      resp = strsplit(resp,"\t")
      
      ## save results in data0 dataframe
      for (r in seq(1,rows)){
        for (cnm in seq(2,cols+1)) {
          this.row = r
          this.col =  cnm - 1
          this.label = label[[r]][cnm]
          this.conc =  as.numeric(conc[[r]][cnm])
          this.resp =  resp[[r]][cnm]
          if (this.resp == ".") {next}
          this.resp = as.numeric(this.resp)

          ### check for excuded wells
          if (is.null(excludeWells)==FALSE) {
            excludeFlag = 0
            for (i in seq(1,length(excludeWells$row))){
              r1 = excludeWells$row[i]
              c1 = excludeWells$col[i]
              if ((this.row==r1) && (this.col==c1)) {
                excludeFlag = 1 
                }
              }
            if (excludeFlag == 1) next
            }
          
          # add well info to data frame
          tmp = list(tray= tray, row=this.row, col=this.col, label=this.label, conc=this.conc, resp=this.resp)
          levels(data0$tray) = c(levels(data0$tray),tray)
          levels(data0$label) = c(levels(data0$label),this.label)
          data0[index,1:6] = tmp
          index = index + 1
        }
      }
    }
  }
  
  data0 = data0[data0$label != 0,]
  returnList = list(concentrationResponse= data0, drugRatios= drugRatios)
  class(returnList) <- "trayData"
  return (returnList)

  }

