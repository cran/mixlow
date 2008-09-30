\name{NlmePrintVarFunctions}
\alias{NlmePrintVarFunctions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print variance function options for nonlinear mixed-effects analysis}
\description{
Prints a table of variance function options that can be used in arguments
to the \code{doNlme} function.
}
\usage{
NlmePrintVarFunctions()
}
%- maybe also 'usage' for other objects documented here.
\arguments{None
}
\details{
One argument to the \code{doNlme} function is specification of the 
variance functions to be used in \code{nlme} analysis.  Eight different 
variance functions are allowed.  The convenience function 
\code{NlmePrintVarFunctions} prints information on each of the eight
possible choices.  There are four choices for analysis type \dQuote{single}
and four for analysis type \dQuote{multiple}.   
}
\value{
None
}

\author{John Boik \email{jcboik@stanford.edu}}

\seealso{
    \code{\link{doNlme}}
}
\examples{
NlmePrintVarFunctions()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }

