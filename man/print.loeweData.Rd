\name{print.loeweData}
\alias{print.loeweData}

\title{ Prints a loeweData object obtained from Loewe analysis}
\description{
Prints a loeweData object returned from the \code{doLoewe} function.
}
\usage{
\method{print}{loeweData}(x, \dots)
}

\arguments{
    \item{x}{A loeweData object returned from the \code{doLoewe} function.}
    \item{\dots}{Optional arguments passed to \code{print.default}, including 
    \code{verbose}, a logical argument with default \code{TRUE}.  Verbose controls the
    amount of output printed. }
}
\details{
Prints a loeweData object returned from the \code{doLoewe} function.
}
\value{
None
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\examples{
# loeweData data object is obtained using the doLoewe function
data(loeweData)
print(loeweData)
}
\keyword{ methods }


