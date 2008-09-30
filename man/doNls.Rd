\name{doNls}
\alias{doNls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conduct nonlinear least-squares analysis of sigmoidal 
concentration-response data}
\description{
Organizes input and calls the \code{nls} function to estimate the
parameters of sigmoidal concentration-response curves.  Users can supply
starting values for parameters or can use default values that are based 
on the data.
}
\usage{
doNls(mixlowData, parameterDefaults, lambdaThreshold=0.05, numRandomTries=10, 
verbose=FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{mixlowData}{A list obtained from \code{prepareData}.}
    \item{parameterDefaults}{A list obtained from 
        \code{getNlsParameterDefaults}. The default values can be edited 
        as desired before sending to the \code{doNls} function.}
    \item{lambdaThreshold}{A scalar designating any lower boundary
        for lambda.  If a model estimates a lambda value less than the 
        threshold, that model will be discarded and a model that does not 
        use a lambda parameter will be used instead.}
    \item{numRandomTries}{A scalar designating the number of attempts
        that should be made to reach convergence for each model.  For each
        attempt, starting values are randomly selected from a normal
        distribution with mean equal to the given parameter starting
        value and a standard deviation equal to the given standard deviation.}
    \item{verbose}{An optional logical value. If \code{TRUE}, information 
        from intermediate steps in the \code{doNls} analysis is printed.  
        Default is \code{FALSE}.}
}
\details{
See \code{doNlme} for a description of the parameters estimated.
See \code{getNlsParameterDefaults} for a description of default 
parameter values.  Parameter values of NaN indicate the following default 
values will be used: param.g=0, param.g.std=1, param.p = the closest tested 
log concentration corresponding to 1/2 the mean response of control wells, 
and param.p.std = the absolute value of param.p/2.  The default for 
param.lambda is NaN, which indicates that a lambda term will not be used.  
If a lambda term is used, its value should be between zero and 0.5. 
}
\value{
A list of type \code{nlsData} of length equal to the number of trays with 
the following components:
    \item{nlsEstimates}{Data frame containing parameter estimates}
    \item{nlsGraphing}{List used by \code{plotNlsData} for making graphs:
        \describe{
        \item{tray}{String tray name}
        \item{drug}{String drug name}
        \item{xx}{Concentration values for graphing predictions}
        \item{pred}{Predicted values}
        \item{y}{Adjusted observed responses}
        \item{x2}{Adjusted observed concentrations}}}
    
Adjusted responses take into account responses in optical control wells
(see \code{prepareData}).  In the adjusted concentrations, values of zero
are replaced by concentrations of 1/1000th of the lowest nonzero 
concentration (for log-plot graphing purposes only).
}

\references{Boik J.C., Newman R.A., Boik R.J. (2008) Quantifying 
synergism/antagonism using nonlinear mixed-effects modeling: a 
simulation study.  Statistics in Medicine 27(7), 1040-61}
\author{John Boik \email{jcboik@stanford.edu}}

\seealso{\code{\link{prepareData}}, \code{\link{plotNlsData}}}
\examples{

# mixlowData data object is obtained using the prepareData function
data(mixlowData)
trays <- getTrays(mixlowData)
parameterDefaults <- getNlsParameterDefaults(trays=trays[1:9])
nlsData <- doNls(mixlowData=mixlowData, parameterDefaults) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models } % Fix me

