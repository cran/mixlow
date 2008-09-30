\name{plotLoeweData}
\alias{plotLoeweData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Loewe analysis results}
\description{
Plots results returned from the \code{doLoewe} function.
}
\usage{
plotLoeweData(loeweData, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{loeweData}{A list returned from the \code{doLoewe} function.}
    \item{...}{Additional parameters, including \code{main} and \code{legend}
        to be passed to the plot function}
}
\details{
Two figures are generated from the Loewe analysis. The first shows the
estimated index and its confidence intervals, along with the reference 
index at 1.0. Values below this reference line are indicative of synergism 
and values above are indicative of antagonism. When the upper and lower 
confidence interval limits both fall below 1.0, then statistically 
significant synergism is occurring. When the upper and lower confidence 
interval limits both fall above 1.0, then statistically significant 
antagonism is occurring. The second figure shows the estimated concentration-
response curves for each drug and the mixture. To save both figures to disk,
sandwich the \code{plotLoeweData} command between \code{pdf} and code{dev.off}
commands, for example.
}
\value{
None
}

\author{John Boik \email{jcboik@stanford.edu}}

\examples{
# loeweData data object is obtained using the doLoewe function
data("loeweData")
plotLoeweData(loeweData, ask=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods } % Fix this

