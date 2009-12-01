\name{print.mixlowData}
\alias{print.mixlowData}

\title{Print a mixlowData object}
\description{
Prints a mixlowData object obtained from the \code{prepareData} function.
}
\usage{
\method{print}{mixlowData}(x, \dots)
}

\arguments{
    \item{x}{A mixlowData object obtained from the \code{prepareData} function.}
    \item{\dots}{Optional arguments passed to \code{print.default}; see
    the documentation on that method function.}
}
\details{
Prints a mixlowData object obtained from the \code{prepareData} function.

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
data(mixlowData)
print(mixlowData)
}

\keyword{ methods } 

