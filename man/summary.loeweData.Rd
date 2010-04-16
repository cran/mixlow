\name{summary.loeweData}
\alias{summary.loeweData}

\title{ Summarizes a loeweData object obtained from Loewe analysis}
\description{
Summarizes a loeweData object returned from the \code{doLoewe} function.
}
\usage{
\method{summary}{loeweData}(object, \dots)
}

\arguments{
    \item{object}{A loeweData object returned from the \code{doLoewe} function.}
    \item{\dots}{Optional arguments passed to \code{summary.default}}
}
\details{
Summarizes a loeweData object returned from the \code{doLoewe} function.
}
\value{
None
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\examples{
# loeweData data object is obtained using the doLoewe function
data(loeweData)
summary(loeweData)
}
\keyword{ methods }


