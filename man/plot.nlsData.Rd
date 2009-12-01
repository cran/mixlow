\name{plot.nlsData}
\alias{plot.nlsData}

\title{Plot nlsData object obtained from nonlinear least-squares analysis }
\description{
  Plots nlsData object obtained from the \code{doNls} function
}
\usage{
\method{plot}{nlsData}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{A nlsData object obtained from the \code{doNls} function }
    \item{\dots}{Additional arguments sent to the \code{plot} function, including: \code{ask}, \code{mixlowData}, 
    \code{trays}, with default getTrays(mixlowData), and \code{showBlanks}, with default \code{FALSE}. 
    All arguments are optional except \code{mixlowData}. Argument \code{ask} is logical. If \code{TRUE} (and 
    the R session is interactive) the user is asked for input before a new figure is drawn. Argument 
    \code{mixlowData} is a mixlowData object obtained from the \code{prepareData} function. Argument \code{trays} 
    is a vector of strings indicating the tray names for which data should be plotted. Argument \code{showBlanks} 
    is logical.  If \code{TRUE}, data from the optical control wells (\dQuote{blanks}) in a tray are included in the graph.  
    The remining portions of the graph show data after adjustment for the optical control well values. Default 
    is \code{FALSE}}    
}
\details{
Plots concentration-response curves for each tray and each drug.  See 
\code{plot.mixlowData} for more information about the plots of 
concentration-response data and optical control data. 
}

\value{
None
}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{
    \code{\link{doNls}, \link{plot.mixlowData}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data(mixlowData)
trays <- getTrays(mixlowData)
# nlsData data object is obtained using the doNls function
data(nlsData)
plot(nlsData, mixlowData=mixlowData, 
  trays= trays[1:9], ask=FALSE, showBlanks= FALSE)
}

\keyword{ methods } 
