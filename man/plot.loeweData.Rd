\name{plot.loeweData}
\alias{plot.loeweData}

\title{ Plot a loeweData object obtained from Loewe analysis}
\description{
Plots the loeweData object returned from the \code{doLoewe} function.
}
\usage{
\method{plot}{loeweData}(x, \dots)
}

\arguments{
    \item{x}{A loeweData object returned from the \code{doLoewe} function.}
    \item{\dots}{Additional optional arguments sent to the \code{plot} function, including: \code{ask},
    \code{main}, \code{xlab}, \code{ylab}, and \code{legend}. Argument \code{ask} is logical. 
    If TRUE (and the R session is interactive) the user is asked for input before a new figure is drawn.}
}
\details{
Two figures are generated from the loeweData object. The first shows the
estimated index and its confidence intervals, along with the reference 
index at 1.0. Values below this reference line are indicative of synergism 
and values above are indicative of antagonism. When the upper and lower 
confidence interval limits both fall below 1.0, then statistically 
significant synergism is occurring. When the upper and lower confidence 
interval limits both fall above 1.0, then statistically significant 
antagonism is occurring. The second figure shows the estimated concentration-
response curves for each drug and the mixture. To save both figures to disk,
sandwich the \code{plot} command between \code{pdf} and \code{dev.off}
commands.
}
\value{
None
}

\author{John Boik \email{jcboik@stanford.edu}}

\examples{
# loeweData data object is obtained using the doLoewe function
data(loeweData)
plot(loeweData, ask=TRUE)
}
\keyword{ methods }


