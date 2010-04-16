\name{summary.mixlowData}
\alias{summary.mixlowData}

\title{Summarizes a mixlowData object}
\description{
Summarizes a mixlowData object obtained from the \code{prepareData} function.
}
\usage{
\method{summary}{mixlowData}(object, \dots)
}

\arguments{
    \item{object}{A mixlowData object obtained from the \code{prepareData} function.}
    \item{\dots}{Optional arguments passed to \code{summary.default}}
}
\details{
Summarizes a mixlowData object obtained from the \code{prepareData} function.

}
\value{
None
}


\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{ 
    \code{\link{readDataFile}}, \code{\link{prepareData}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data(mixlowData)
summary(mixlowData)
}

\keyword{ methods } 

