\name{print.nlmeData}
\alias{print.nlmeData}

\title{Print a nlmeData object obtained from nonlinear mixed-effects analysis}
\description{
Prints a nlmeData object obtained from \code{doNlme} function
}
\usage{
\method{print}{nlmeData}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{A nlmeData object obtained from the \code{doNlme} function }
    \item{\dots}{Optional arguments passed to \code{print.default}, including 
    \code{verbose}, a logical argument with default \code{TRUE}.  Verbose controls the
    amount of output printed. }
}
\details{
Prints a nlmeData object obtained from \code{doNlme} function
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
print(nlmeData)
}

\keyword{ methods } 
