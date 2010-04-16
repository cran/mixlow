\name{summary.nlsData}
\alias{summary.nlsData}

\title{Summarizes a nlsData object obtained from nonlinear least-squares analysis }
\description{
  Summarizes a nlsData object obtained from the \code{doNls} function
}
\usage{
\method{summary}{nlsData}(object, \dots)
}

\arguments{
    \item{object}{A nlsData object obtained from the \code{doNls} function }
    \item{\dots}{Optional arguments passed to \code{summary.default}}
}
\details{
Summarizes a nlsData object obtained from the \code{doNls} function
}

\value{
None
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{
    \code{\link{doNls}}
}
\examples{
# nlsData data object is obtained using the doNls function
data(nlsData)
summary(nlsData)
}

\keyword{ methods } 
