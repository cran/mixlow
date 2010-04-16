\name{prepareData}
\alias{prepareData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prepares raw data for use by other functions}
\description{
Prepares data for use by \code{doNls} and other functions.  
}
\usage{
prepareData(trayData, drugs, trays, cellLines, degree=3)
}

\arguments{
    \item{trayData}{A list returned from \code{readDataFile}}
    \item{drugs}{A vector of drug names used to subset data}
    \item{trays}{A vector of tray names used to subset data}
    \item{cellLines}{A vector of cell line names used to subset data}
    \item{degree}{The degree of the polynomial used in the model for bbc blank wells. Must be in
    the set [1,2,3,4]}
}
\details{
Data contained in a trayData object (produced by use of the 
\code{readDataFile} function) must be adjusted before it can be used 
by other functions of this package.  The \code{prepareData} function adjusts
observed responses based on responses in the optical control (\dQuote{blank}) 
wells.  The method used to adjust observed responses is stated in the data 
file (methods are either \dQuote{bbt} or \dQuote{bbc}).  For method \dQuote{bbc}, 
a 1st to 4th degree polynomial is fit to the concentration-dependent optical control 
data, depending on the value of the argument \code{degree}.  Results of using different
polynomial degrees can be seen by using the function \code{plot.mixlowData} with argument 
\code{showBlanks} set to \code{TRUE}.  Predictions for each concentration are subtracted from 
the observed experimental responses.  For type \dQuote{bbt}, the average 
optical control response is subtracted from all observed responses.

To aid in forming the arguments \code{drugs}, \code{trays}, and 
\code{cellLines}, vectors of each can be obtained using the convenience functions 
\code{getDrugs}, \code{getTrays}, and \code{getCellLines}. 

}
\value{
A list of class \code{mixlowData} with the following components:
  \item{concentrationResponse}{A data frame containing adjusted 
    concentration-response data}
  \item{drugRatios}{A data frame containing drug ratios}
  \item{plottingData}{A list used by \code{plot.mixlowData} in plotting}
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{ 
    \code{\link{readDataFile}}, \code{\link{plot.mixlowData}}
}
\examples{
# trayData data object is obtained using the readDataFile function
data(trayData)
trays <- getTrays(trayData)
cellLines <- getCellLines(trayData)
mixlowData <- prepareData(trayData=trayData, trays=trays[1:9], cellLines=cellLines[1])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }

