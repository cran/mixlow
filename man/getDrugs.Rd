\name{getDrugs}
\alias{getDrugs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve a vector of drug names from data}
\description{
Convenience function called on output from \code{prepareData} or 
\code{readDataFile} to obtain a vector of drug names. 
}
\usage{
getDrugs(data, trays=NULL, cellLines=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{A list obtained from \code{prepareData} or \code{readDataFile} 
    functions.}
    \item{trays}{A vector of tray names.  If \code{NULL}, drugs will not
        be subset by trays.} 
    \item{cellLines}{A vector of cell line names.  If \code{NULL}, drugs
        will not be subset by cell lines.}
}
\details{
A data file may contain experimental results for multiple drugs.  
If so, slices of the vector returned from \code{getDrugs} can be 
used as an argument to the \code{prepareData} function in order to select
a subset of the data.
}
\value{
A vector of strings indicating which drugs are contained in the input data 
object. 
}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{
    \code{\link{getTrays}}, \code{\link{getCellLines}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data("mixlowData")
drugs <- getDrugs(mixlowData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print } % Fix me

