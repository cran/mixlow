\name{nlsData}
\alias{nlsData}
\docType{data}
\title{nlsData Dataset}
\description{
  Dataset obtained as a result of using function \code{doNls}, based on data 
  in the text file \dQuote{A549\_vin\_topo\_data.txt}.  
}
\usage{data(nlsData)}
\format{
See help file for function \code{doNls} for composition of nlsData
dataset.
}
\details{
  The text file \dQuote{A549\_vin\_topo\_data.txt} provides example data for 
  a cytotoxicity experiment in which two drugs, vincristine and topotecan, and 
  their 1:1 mixture is tested in A549 human lung cancer cells. Each drug is 
  tested in three replicate trays.  See reference for function \code{doLoewe}.
  }


\seealso{\code{\link{doNls}}}
            

\examples{
data(nlsData)
}
\keyword{datasets}
