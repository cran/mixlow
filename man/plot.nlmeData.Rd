\name{plot.nlmeData}
\alias{plot.nlmeData}

\title{Plot nlmeData object obtained from nonlinear mixed-effects analysis}
\description{
Plots nlmeData object obtained from \code{doNlme}  function
}
\usage{
\method{plot}{nlmeData}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{A nlmeData object obtained from the \code{doNlme} function }
    \item{\dots}{Additional optional arguments sent to the \code{plot} function, including: \code{ask},
    \code{main}, \code{xlab}, \code{ylab}, and \code{legend}. Argument \code{ask} is logical. 
    If TRUE (and the R session is interactive) the user is asked for input before a new figure is drawn.}
}
\details{
Four sets of figures are generated per \code{nlme} analysis.  The first
is concentration-response curves for all trays treated with a given drug,
the second is qqnorm plots of residuals for each drug, the third is residual
versus fitted values for each drug, and the fourth is standardized residuals
versus fitted values for each drug. 
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
plot(nlmeData, ask=TRUE)
}

\keyword{ methods } %% Fix me
