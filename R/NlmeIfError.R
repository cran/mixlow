`NlmeIfError` <-
function(x)
  {
  ## handler for NLME errors
  print(x$message)
  txt = x$message
  gr = grep("Singular precision",txt)
  if (length(gr)>0)  stop("warn to error: singlular precision...")
  }

