\name{mixlow-package}
\alias{mixlow-package}
\alias{mixlow}
\docType{package}
\title{
Calculate synergism/antagonism indices based on concentration-
response curve parameters 
}
\description{
Uses a mixed-effects model to calculate parameters of concentration-
response curves.  These parameters are used in calculating Loewe 
interaction indices to quantify drug synergism/antagonism. 
}
\details{
\tabular{ll}{
Package: \tab mixlow\cr
Type: \tab Package\cr
Version: \tab 0.02\cr
Date: \tab 2010-4-8\cr
License: \tab 	GPL-2\cr
}
Functions are generally called in this order: \code{readDataFile}, 
\code{prepareData}, \code{doNls}, \code{doNlme}, and \code{doLoewe}. 
Plot functions are available to view results graphically.
}
\author{
John Boik \email{john.boik@newearthbiomed.org}
}
\references{
  Boik J.C., Newman R.A., Boik R.J. (2008).
  Quantifying Synergism/Antagonism Using Nonlinear Mixed-Effects Modeling:
  A Simulation Study. \emph{Statistics in Medicine}, \bold{27}(7), 1040--1061

  Boik J.C., Narasimhan N. (2010).
  An R Package for Assessing Drug Synergism/Antagonism.
  \emph{Journal of Statistical Software}, \bold{34}(6), 1--18.
}
\examples{
# trayData data object is obtained using the readDataFile function
data(trayData)
mixlowData <- prepareData(trayData)
trays <- getTrays(trayData)
parameterDefaults <- getNlsParameterDefaults(trays[1:9])
parameterDefaults["vin_tr1","param.lambda"] = .2
parameterDefaults["vin_tr2","param.lambda"] = .2
parameterDefaults["vin_tr3","param.lambda"] = .2
nlsData = doNls(mixlowData, parameterDefaults) 
nlmeData = doNlme(mixlowData, nlsData)
loeweData = doLoewe(mixlowData, nlmeData)

}
\keyword{ package }

