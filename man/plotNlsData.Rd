\name{plotNlsData}
\alias{plotNlsData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot results obtained from nonlinear least-squares analysis }
\description{
  Plots results from the \code{doNls} function
}
\usage{
plotNlsData(nlsData, mixlowData, trays=getTrays(mixlowData), ask, 
    showBlanks=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{nlsData}{A list obtained from the \code{doNls} function }
    \item{mixlowData}{A list obtained from the \code{prepareData} function }
    \item{trays}{A vector of strings indicating the trays for which data
        should be plotted.}
    \item{ask}{Logical. If \code{TRUE} (and the R session is interactive) 
        the user is asked for input, before a new figure is drawn.}
    \item{showBlanks}{Logical. If \code{TRUE}, data from the optical control 
        wells (\dQuote{blanks}) in a tray are included in the graph.  
        Default is \code{TRUE}} 
}
\details{
Plots concentration-response curves for each tray and each drug.  See 
\code{plotMixlowData} for more information about the plots of 
concentration-response data and optical control data.
}
\value{
A list of length equal to the number of trays with components:
    
    \item{graphDataTray}{A list of graph data per tray with components:
        \describe{
        \item{tray}{Tray name}
        \item{drug}{Drug name}
        \item{cell}{Cell line name}
        \item{Units}{Units of measurement}
        \item{title}{Title for plot}
        \item{adjResp}{Adjusted response data}
        \item{adjConc}{Adjusted concentration data}
        \item{ylim}{y limits of plot}
        \item{linesY}{Response values for mean line between group replicates}
        \item{linesX}{Concentration values for mean line between group 
          replicates}
        \item{blankPointsY}{Responses from optical control wells}
        \item{blankPointsX}{Concentrations for optical control wells}
        \item{blankLineY}{Responses for optical control data line}
        \item{blankLineX}{Concentrations for optical control data line}
        \item{predLinesY}{Predicted responses}
        \item{predLinesX}{Concentrations for predicted responses}}}

    \item{graphDataDrug}{A list of graph data per drug with components:
        \describe{
        \item{trays}{Tray names}
        \item{Units}{Units of measurement}
        \item{title}{Title for plot}
        \item{drug}{Drug name}
        \item{y}{Adjusted responses}
        \item{x}{Adjusted concentrations}
        \item{ylim}{y limits of plot}
        \item{lineList}{List of x and y coordinates for each line of plot}
        \item{legend}{Legend for plot}}}

}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{
    \code{\link{doNls}, \link{plotMixlowData}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data(mixlowData)
trays <- getTrays(mixlowData)
# nlsData data object is obtained using the doNls function
data(nlsData)
saved_nls_GraphingData <- plotNlsData(nlsData=nlsData, mixlowData=mixlowData, 
  trays= trays[1:9], ask=FALSE, showBlanks= TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods } % Fix me
