\name{doLoewe}
\alias{doLoewe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate Loewe indices and confidence intervals}
\description{
Uses parameter and covariance estimates produced by
\code{doNlme} to estimate Loewe indices and confidence
intervals at different fraction affected values. 
}
\usage{
doLoewe(mixlowData, nlmeData, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{mixlowData}{A list obtained from the output of 
        \code{prepareData}}
    \item{nlmeData}{A list obtained from the output of \code{doNlme}}
    \item{verbose}{An optional logical value. If \code{TRUE}, information 
        from intermediate steps in the Loewe analysis is printed.  Default 
        is \code{FALSE}.}
}
\details{
Loewe indices are estimated based on parameter and covariance values generated
by use of the \code{doNlme} function. In addition to returning interaction 
indices, the \code{doLoewe} function returns a vector indicating the degree 
of statistically significant synergism at each fraction affected value.  The degree of 
statistically significant synergism at each fraction affected value is 
zero if the upper confidence limit is above 1.0, and is the difference 
between the upper confidence limit and 1.0 otherwise.  A vector of 
statistically significant antagonism is calculated in an analogous 
manner.  These vectors summarize the degree of synergism/antagonism
occurring over a range of fraction affected values. 

The width of the index confidence intervals is dependent in part on the 
critical value taken from the t-distribution, which in turn is dependent 
on the degrees of freedom.  For multiple drug analysis with \code{doNlme}, 
the degrees of freedom is taken directly from the \code{nlme} results.  
For single drug analysis, the degrees of freedom is taken as the minimum 
degrees of freedom over all drugs. 

}
\value{
A list with the following named components: 
    \item{drugs}{The drugs analyzed}
    \item{mix}{The mixture analyzed}
    \item{ciL}{A data frame containing Loewe index estimates and confidence 
        intervals based on the log of the index, for different fraction 
        affected values}
    \item{covv}{Covariance matrix of the parameter estimates}
    \item{g0}{gamma, shape parameters}
    \item{p0}{psi, IC50 parameters}
    \item{uu}{u parameter}
    \item{dat0}{Data frame containing concentration-response values}
    \item{score.interval}{Interval of fraction affected values over which
        synergism/antagonism scores were computed}
    \item{syner}{Synergism scores for fraction affected values in interval}
    \item{syner.total}{Sum of synergism scores over interval}
    \item{antag}{Antagonism scores for fraction affected values in interval}
    \item{antag.total}{Sum of antagonism scores over interval}

See \code{doNlme} for a description of gamma, psi, and u parameters.

}
\references{The model is described in Boik J.C., Newman R.A., Boik R.J.
(2008) Quantifying synergism/antagonism using nonlinear mixed-effects 
modeling: a simulation study.  Statistics in Medicine 27(7), 1040-61}
\author{John Boik \email{jcboik@stanford.edu}}

\seealso{
    \code{\link{doNlme}}, \code{\link{plotLoeweData}}
}
\examples{
# mixlowData data object is obtained using the prepareData function
data(mixlowData)
# nlmeData data object is obtained using the doNlme function
data(nlmeData)
loeweData <- doLoewe(mixlowData, nlmeData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models } %% Fix me

