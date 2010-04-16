\name{summary.nlmeData}
\alias{summary.nlmeData}

\title{Summarizes a nlmeData object obtained from nonlinear mixed-effects analysis}
\description{
Summarizes a nlmeData object obtained from \code{doNlme} function
}
\usage{
\method{summary}{nlmeData}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{A nlmeData object obtained from the \code{doNlme} function }
    \item{\dots}{Optional arguments passed to \code{summary.default}}
}
\details{
Summarizes a nlmeData object obtained from \code{doNlme} function
}
\value{
None
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{ 
    \code{\link{doNlme}}
}
\examples{
# nlmeData data object is obtained using the doNlme function
data(nlmeData)
summary(nlmeData)
}

\keyword{ methods } 
