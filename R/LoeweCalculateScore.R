`LoeweCalculateScore` <-
function(ciL)
  {
  # calculates a interaction score as area under the index curve of significant antagoism - significant synergism
  
  
  maxx = max(ciL$Fraction.Affected) -.03
  score.interval = c(0.05,maxx)
  antag = numeric()
  syner = numeric()
  antag2 = 0
  syner2 = 0
  L = numeric()
  for (i in seq(length(ciL[,1]))) {
    if ((ciL[i,1] < score.interval[1]) | (ciL[i,1] > score.interval[2])) next
    low1 = ciL[i,4]
    hi1 = ciL[i,5]
    low2 = ciL[i+1,4]
    hi2 = ciL[i+1,5]
  
    x1 = ciL[i,1]
    x2 = ciL[i+1,1]
    L = c(L, ciL[i,2])
  
    {if ((hi1 < 1) & (hi2 < 1)){
      h = c(1-hi1, 1-hi2)
      h = sort(h)
      syner = c(syner, (h[1]*(x2-x1) + 0.5*(x2-x1)*(h[2]-h[1])))    # area under the curve inerval estimate
      syner2 = syner2 + (h[1]*(x2-x1) + 0.5*(x2-x1)*(h[2]-h[1]))
      }
    else {
      syner = c(syner,0)
      }}
  
    {if ((low1 > 1) & (low2 > 1)){
      h = c(low1-1, low2-1)
      h = sort(h)
      antag = c(antag, (h[1]*(x2-x1) + 0.5*(x2-x1)*(h[2]-h[1])))
      antag2 = antag2 + (h[1]*(x2-x1) + 0.5*(x2-x1)*(h[2]-h[1]))
      }
    else {
      antag = c(antag, 0)
      }}
    }
  return (list(score.interval= score.interval, antag=antag, antag2=antag2, syner=syner, syner2=syner2))
  }

