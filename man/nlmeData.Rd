\name{nlmeData}
\alias{nlmeData}
\docType{data}
\title{nlmeData Dataset}
\description{
  Dataset obtained as a result of using function \code{doNlme}, based on data 
  in the text file \dQuote{A549_vin_topo_data.txt}.  
}
\usage{data(nlmeData)}
\format{
See help file for function \code{doNlme} for composition of nlmeData
dataset.
}
\details{
  The text file \dQuote{A549_vin_topo_data.txt} provides example data for 
  a cytotoxicity experiment in which two drugs, vincristine and topotecan, and 
  their 1:1 mixture is tested in A549 human lung cancer cells. Each drug is 
  tested in three replicate trays.  See reference for function \code{doLoewe}.
  }

\seealso{\code{\link{doNlme}}}

\examples{
data(nlmeData)
}
\keyword{datasets}
