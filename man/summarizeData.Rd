\name{summarizeData}
\alias{summarizeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print summary of data object}
\description{
Prints a summary of a mixlowData or trayData data object.  
}
\usage{
summarizeData(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A list obtained from the \code{prepareData} or \code{readDataFile}
  functions.}  
}
\details{
Prints a summary of information for data objects returned by the 
\code{prepareData} or or \code{readDataFile} functions.
}
\value{
None 
}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{ 
    \code{\link{readDataFile}}, \code{\link{prepareData}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data("mixlowData")
summarizeData(mixlowData)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods } %% Fix me

