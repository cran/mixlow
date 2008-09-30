\name{doNlme}
\alias{doNlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conduct nonlinear mixed-effects analysis of sigmoidal concentration-
response data
}
\description{
Organizes input and calls the \code{nlme} function to obtain parameter
estimates for sigmoidal concentration-response curves.  Starting values
for the fixed effects are obtained from the output of the
\code{doNls} function. 
}
\usage{
doNlme(mixlowData, nlsData, drugs=getDrugs(mixlowData), analysis="single", 
    varFunction= 1, method="ML", verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{mixlowData}{A list obtained from the output of
        \code{prepareData}.  It contains adjusted concentration-response
        data. }
    \item{nlsData}{A list obtained from the output of
        \code{doNls}.  It contains parameter estimates for the
        concentration-effect curves, which are used as fixed-effects
        starting values by the \code{nlme} function. }
    \item{drugs}{A vector of drug names}
    \item{analysis}{An optional character string either \dQuote{single}
        or \dQuote{multiple} where \dQuote{single} indicates that each 
        drug is to be analyzed separately and \dQuote{multiple} indicates 
        that all drugs should be analyzed together.  The default is 
        \dQuote{single}.} 
    \item{varFunction}{An optional numerical vector or list of lists 
        from the set (1,2,3,4) that specifies the variance function(s) 
        to be used for each analysis.  If \code{analysis} is \dQuote{multiple},
        \code{varFunction} can be a vector.  If \code{analysis} is 
        \dQuote{single} it can also be a vector, in which case all drugs
        will use the same varFunction vector.  Or, if \code{analysis} is 
        \dQuote{single}, a list can be supplied of length equal to the 
        number of drugs. Each entry in the list should be a list with named
        component varFunction.  This allows a different set of varFunctions
        to be used with each drug. The default is 1.}
    \item{method}{An optional character string.  If \dQuote{REML} the 
        model is fit by maximizing the restricted log-likelihood.  If 
        \dQuote{ML} the log-likelihood is maximized.  Default is  
        \dQuote{ML}.} 
    \item{verbose}{An optional logical value. If \code{TRUE} information 
        from intermediate steps in the \code{doNlme} analysis is printed.  
        The default is \code{FALSE}.}
    }

\details{
This function uses output from \code{doNls} as starting values for 
fixed effects.  It organizes input and sets up calls to the \code{nlme}
function in order to estimate parameters of sigmoidal concentration-
response curves.

The variance functions for analysis of single drugs are: 
    \describe{
    \item{1}{\code{sigma}}
    \item{2}{\code{sigma*E[response]}}
    \item{3}{\code{sigma*E[response]^beta}}
    \item{4}{\code{sigma*(beta1 + E[response])}}}

The variance functions for analysis of multiple drugs are:
    \describe{
    \item{1}{\code{sigma}}
    \item{2}{\code{sigma*alpha}, where alpha is drug-dependent}
    \item{3}{\code{sigma*E[response]}}
    \item{4}{\code{sigma*E[response]^beta}, where beta is drug-dependent}}
    
\code{E[response]} is used above to designate the expected response in a given 
well of a tray.  For heteroscedastic errors, varFunction=2 (for single) 
and varFunction=3 (for multiple) analysis may be appropriate in many cases.  
The simpler error functions could be tried if the \code{nlme} function 
does not converge.

See Boik et al., 2008 for details of the model.  The parameters
estimated by the model are:
    \describe{
    \item{g}{gamma, designates the steepness of the concentration-effect 
        curve at the IC50}
    \item{p}{psi, designates the IC50}
    \item{u}{\code{E[exp(u) + b_t]}, the expected value of all control wells 
        across trays, where b is a tray-dependent random variable}
    \item{lambda}{Optional. The non-zero asymptote of a \dQuote{modified}
        sigmoidal concentration-response curve.  The value of lambda 
        represents the fraction of \code{E[exp(u) + b_t]} that is associated 
        with the non-zero asymptote.}  }

Parameters u, g, and p are in log scale.

Use of analysis=\dQuote{multiple} for mixtures that contain more than 
two drugs can sometimes be problematic in that the estimation procedure 
may not converge or may take a long time to converge. If a mixture 
contains more than a few drugs, one alternative is to estimate 
concentration-response curve parameters for each drug/mixture separately
(i.e., use analysis=\dQuote{single}).  The disadvantage with this 
approach is that for analysis of any given drug, control-well data from 
the trays of other drugs are not used.  Control-well responses are then 
estimated based on only a few (replicate) trays.
}
\value{
Returns a list of class \code{nlmeData} for each analysis with the 
following components:
    
    \item{nlmeResults}{A list of results from the \code{doNlme} analysis:    
        \describe{
        \item{nam}{Name of the best model, chosen according to BIC score}
        \item{method}{Method used to estimate the model}
        \item{drug}{A string containing drug names}
        \item{setNum}{The sequential number of the analysis}
        \item{cell}{String indicating the cell line}
        \item{lik}{Log likelihood}
        \item{bic}{BIC score}
        \item{sig}{sigma}
        \item{modelstruct}{A list containing \code{log(stdev(b_t)/sigma)} and 
          other variance structures}
        \item{rc}{Goodness-of-fit statistic}
        \item{r2}{\code{r^2} goodness-of-fit statistic}
        \item{se}{Standard errors}
        \item{coeff}{Fixed effects}
        \item{df}{Degrees of freedom}
        \item{covv}{Parameter covariance values}}}
    \item{nlmeGraph}{A list of data used for graphing \code{nlme} results:
        \describe{
        \item{pred0}{Predicted values}
        \item{dat1}{Data values}
        \item{ord}{A string containing drug names}
        \item{residu}{Residuals from the best model}
        \item{best}{A string denoting the name of the best model}}}


}
\references{The model is described in Boik J.C., Newman R.A., Boik R.J.
(2008) Quantifying synergism/antagonism using nonlinear mixed-effects 
modeling: a simulation study.  Statistics in Medicine 27(7), 1040-61}
\author{John Boik \email{jcboik@stanford.edu}}
\seealso{\code{\link{doNls}}, \code{\link{plotNlmeData}}}
\examples{
# mixlowData data object is obtained using the prepareData function
data(mixlowData)
# nlsData data object is obtained using the doNls function
data(nlsData)
drugs = getDrugs(mixlowData)
nlmeData = doNlme(mixlowData, nlsData, drugs=drugs, varFunction= c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models } % Fix me
