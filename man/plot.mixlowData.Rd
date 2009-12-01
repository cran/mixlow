\name{plot.mixlowData}
\alias{plot.mixlowData}

\title{Plot a mixlowData object}
\description{
Plots mixlowData object obtained from the \code{prepareData} function.
}
\usage{
\method{plot}{mixlowData}(x, \dots)
}

\arguments{
    \item{x}{A mixlowData object obtained from the \code{prepareData} function.}
    \item{\dots}{Additional optional arguments sent to the \code{plot} function, including: \code{ask}, 
    \code{trays}, with default getTrays(x), and \code{showBlanks}, with default \code{FALSE}. Argument 
    \code{ask} is logical. If \code{TRUE} (and the R session is interactive) the user is asked for 
    input before a new figure is drawn. Argument \code{trays} is a vector of strings indicating the 
    tray names for which data should be plotted. Argument \code{showBlanks} is logical. If \code{TRUE}, 
    data from the optical control wells (\dQuote{blanks}) in a tray are included in the graph.  The 
    remining portions of the graph show data after adjustment for the optical control well values. Default 
    is \code{FALSE}}
}
\details{
This function plots a mixlowData object.  Adjusted concentration-response 
data are plotted, and mean responses over replicate wells 
receiving identical drug concentrations are connected by a line.  
The adjustments subtracted from the raw data (see \code{PrepareData} function)
can be included in the plot by setting \code{showBlanks} equal to \code{TRUE}.
For purposes of plotting on a log-scale plot, drug concentrations of zero are changed to 
1/1000 of the lowest nonzero drug concentration.  This adjustment only
affects plotting, not analysis. To save both figures to disk,
sandwich the \code{plot} command between \code{pdf} and \code{dev.off}
commands.

}
\value{
None
}


\author{John Boik \email{jcboik@stanford.edu}}

\seealso{ 
    \code{\link{readDataFile}}, \code{\link{prepareData}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data(mixlowData)
trays <- getTrays(mixlowData)
plot(mixlowData, trays= trays[1:9], ask=TRUE, showBlanks= FALSE)
}

\keyword{ methods } % Fix this

