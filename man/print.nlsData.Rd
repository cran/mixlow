\name{print.nlsData}
\alias{print.nlsData}

\title{Prints a nlsData object obtained from nonlinear least-squares analysis }
\description{
  Prints a nlsData object obtained from the \code{doNls} function
}
\usage{
\method{print}{nlsData}(x, \dots)
}

\arguments{
    \item{x}{A nlsData object obtained from the \code{doNls} function }
    \item{\dots}{Optional arguments passed to \code{print.default}; see
    the documentation on that method function.}
}
\details{
Prints a nlsData object obtained from the \code{doNls} function
}

\value{
None
}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{
    \code{\link{doNls}}
}
\examples{
# nlsData data object is obtained using the doNls function
data(nlsData)
print(nlsData)
}

\keyword{ methods } 
