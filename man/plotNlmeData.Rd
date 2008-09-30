\name{plotNlmeData}
\alias{plotNlmeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot results obtained from nonlinear mixed-effects analysis}
\description{
Plots results obtained from \code{doNlme}  function
}
\usage{
plotNlmeData(nlmeData, ask, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{nlmeData}{A list obtained from the \code{do.NLME} function }
    \item{ask}{Logical. If TRUE (and the R session is interactive) the user
        is asked for input, before a new figure is drawn.}
    \item{...}{Additional parameters, including \code{main} and \code{legend} 
        to be passed to the plot function}
}
\details{
Two sets of figures are generated per \code{nlme} analysis.  The first
is concentration-response curves for all trays treated with a given drug,
and the second is qqnorm plots of residuals for each drug. 
}
\value{
None
}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{ 
    \code{\link{doNlme}}
}
\examples{
# nlmeData data object is obtained using the doNlme function
data("nlmeData")
plotNlmeData(nlmeData=nlmeData, ask=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods } %% Fix me
