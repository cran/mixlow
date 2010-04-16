\name{summary.trayData}
\alias{summary.trayData}

\title{Summarizes a trayData object}
\description{
Summarizes a trayData object obtained from the \code{readDataFile} function.
}
\usage{
\method{summary}{trayData}(object, \dots)
}

\arguments{
    \item{object}{A trayData object obtained from the \code{readDataFile} function.}
    \item{\dots}{Optional arguments passed to \code{summary.default}}
}
\details{
Summarizes a trayData object obtained from the \code{readDataFile} function.

}
\value{
None
}


\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{ 
    \code{\link{readDataFile}}
}
\examples{
# trayData data object is obtained using the readDataFile function
data(trayData)
summary(trayData)
}

\keyword{ methods } 












