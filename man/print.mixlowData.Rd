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
    \item{\dots}{Optional arguments passed to \code{print.default}, including 
    \code{verbose}, a logical argument with default \code{TRUE}.  Verbose controls the
    amount of output printed. }
}
\details{
Prints a mixlowData object obtained from the \code{prepareData} function.

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
print(mixlowData)
}

\keyword{ methods } 

