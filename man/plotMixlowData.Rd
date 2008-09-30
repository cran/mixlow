\name{plotMixlowData}
\alias{plotMixlowData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot prepared data}
\description{
Plots results obtained from the \code{prepareData} function.
}
\usage{
plotMixlowData(mixlowData, trays=getTrays(mixlowData), ask, 
    showBlanks= TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{mixlowData}{A list obtained from the \code{prepareData} function.}
    \item{trays}{A vector of strings indicating the trays for which data
        should be plotted.}
    \item{ask}{Logical. If TRUE (and the R session is interactive) the user
      is asked for input, before a new figure is drawn.}
    \item{showBlanks}{Logical. If TRUE, data from the optical control wells
        (\dQuote{blanks}) in a tray are included in the graph.  Default is 
        \code{TRUE}}  

}
\details{
This function is used to plot prepared data.  Adjusted concentration-response 
data are plotted, and mean responses over replicate wells 
receiving identical drug concentrations are connected by a line.  
The adjustments subtracted from the raw data (see \code{PrepareData} function)
can be included in the plot by setting \code{showBlanks} equal to \code{TRUE}.
For purposes of plotting on a log-scale plot, drug concentrations of zero are changed to 
1/1000 of the lowest nonzero drug concentration.  This adjustment only
affects plotting, not analysis.

}
\value{
A list of length equal to the number of trays with components:
    \item{drug}{Drug name}
    \item{cell}{Cell line name}
    \item{Units}{Units of measurement}
    \item{title}{Title for plot}
    \item{adjResp}{Adjusted response data}
    \item{adjConc}{Adjusted concentration data}
    \item{ylim}{y limits of plot}
    \item{meanLineY}{Response values for mean line between group replicates}
    \item{meanLineX}{Concentration values for mean line between group 
      replicates}
    \item{blankY}{Responses from optical control wells}
    \item{blankX}{Concentrations for optical control wells}
    \item{blankLineY}{Responses for optical control data line}
    \item{blankLineX}{Concentrations for optical control data line}

}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{ 
    \code{\link{readDataFile}}, \code{\link{prepareData}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data(mixlowData)
trays <- getTrays(mixlowData)
saved_mixlow_graphing_data <- plotMixlowData(mixlowData=mixlowData, 
    trays= trays[1:9], ask=FALSE, showBlanks= TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods } % Fix this

