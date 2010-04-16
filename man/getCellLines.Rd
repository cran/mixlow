\name{getCellLines}
\alias{getCellLines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve a vector of cell line names from data}
\description{
Convenience function called on output from \code{prepareData} or
\code{readDataFile} to obtain a vector of cell line names. 
}
\usage{
getCellLines(data, drugs=NULL, trays=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{A list obtained from \code{prepareData} or \code{readDataFile} 
    functions}
    \item{drugs}{A vector of drug names.  If \code{NULL}, cell lines will 
        not be subset by drugs.}
    \item{trays}{A vector of tray names.  If \code{NULL}, cell lines will
        not be subset by trays.} 
}
\details{
A data file may contain experimental results for multiple cell lines.  
If so, the \code{prepareData} function requires that one of these be 
specified for further analysis.  Slices of the vector returned from 
\code{getCellLines} can be used as an argument to the \code{prepareData} 
function to subset the data.
}
\value{
A vector of strings indicating which cell lines are contained in the input
data object. 
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{
    \code{\link{getTrays}}, \code{\link{getDrugs}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data("mixlowData")
cellLines <- getCellLines(mixlowData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming } % Fix me

