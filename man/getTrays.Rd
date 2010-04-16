\name{getTrays}
\alias{getTrays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve a vector of tray names from data}
\description{
Convenience function called on output from \code{prepareData} or 
\code{readDataFile} to obtain a vector of tray names. 
}
\usage{
getTrays(data, drugs=NULL, cellLines=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{A list obtained from \code{prepareData} or \code{readDataFile}
    functions.}
    \item{drugs}{A vector of drug names.  If \code{NULL}, trays will not
        be subset by drugs.}
    \item{cellLines}{A vector of cell line names.  If \code{NULL}, trays
        will not be subset by cell lines.} 
}
\details{
A data file will contain experimental results for multiple trays, with 
one drug tested per tray.  Slices of the vector returned from 
\code{getTrays} can be used as an argument to the \code{prepareData} 
function in order to select a subset of the data.
}
\value{
A vector of strings indicating which trays are contained in the input data 
object. 
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{
    \code{\link{getCellLines}}, \code{\link{getDrugs}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data("mixlowData")
trays <- getTrays(mixlowData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming } % Fix me

