\name{mixlowData}
\alias{mixlowData}
\docType{data}
\title{mixlowData Dataset}
\description{
  Dataset obtained as a result of using function \code{prepareData}, based on data 
  in the text file \dQuote{A549_vin_topo_data.txt}.  
}
\usage{data(mixlowData)}
\format{
See help file for function \code{prepareData} for composition of mixlowData
dataset.

}
\details{
  The text file \dQuote{A549_vin_topo_data.txt} provides example data for 
  a cytotoxicity experiment in which two drugs, vincristine and topotecan, and 
  their 1:1 mixture is tested in A549 human lung cancer cells. Each drug is 
  tested in three replicate trays.  See reference for function \code{doLoewe}.
  }


\seealso{\code{\link{prepareData}}}
           

\examples{
data(mixlowData)
}
\keyword{datasets}
