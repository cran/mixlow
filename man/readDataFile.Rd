\name{readDataFile}
\alias{readDataFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a formatted file of concentration-response data}
\description{
Reads a formatted data file that contains concentration-response information.  
}
\usage{
readDataFile(filename, excludeWells)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{filename}{A string designating the data file to be read.}
    \item{excludeWells}{An optional list that specifies cells 
        that should be skipped in the analysis.  Cells are identified by row and column numbers.  Paired components of 
        the list are:
        \describe{
        \item{row}{skip a cell in row}
        \item{col}{skip a cell in column}} 
        } 
}
\details{
The data file read by this function must be formatted according to a
specific structure.  See the example file mixlowExampleData.txt for
use as a template.  The data file is tab-delimited pure ASCII (without
enclosing quotations) and the first entry on each line is a line label.
The first eight lines contain general information about the experiments.
Of these lines, only the entries \code{conc_units}, \code{rows}, and 
\code{cols} are required.  The remainder of the data file is split into 
blocks, one block for each tray.  For each block, the required entries are
\code{tray_label}, \code{cell_line}, \code{drug_name_short}, 
\code{composition}, as well as the entries for \code{conc} (concentration), 
\code{label}, and \code{resp} (response) for each row and column
of a tray.  Even if the optional entries are left blank, the line
labels must still be included (including the \dQuote{space} labels). 

The \dQuote{ray} design should be used for each experiment.  In this design,
various dilutions of a drug or mixture are tested in the wells.  If a
mixture is being tested, a fixed ratio between component drugs occurs at 
every concentration tested.  A tray should contain data for only a single
drug or mixture.  Replicate trays should be analyzed for each drug or 
mixture (for a given cell line).  At least three replicate trays are 
recommended.  However, if only one tray is available for a drug,
\code{doNlme} will duplicate that tray (to make two replicates)
so that analysis can be done in \dQuote{experimental} mode.  Tray labels
can be alphanumeric, but avoid any characters other than digits and letters
(the underscore can be used).  The same is true for the abbreviated drug
names that are required.  

Each tray should contain wells used as treatment controls (i.e., those
that receive cells but no drug).  In addition, trays should contain
wells used as optical controls (called \dQuote{blanks} here).  These can
be of two types: 1) wells that contain only media, and 2) wells at each drug
concentration that contain media and drug.  The first type is referred
to as \dQuote{bbt}, or blanks-by-tray, and the second type is referred to as
\dQuote{bbc}, or blanks-by-concentration.  The use of \dQuote{bbc} is 
recommended when concentration-dependent responses can be induced by the 
drug alone (e.g., autofluorescence).  For type \dQuote{bbt}, raw responses 
are adjusted by subtracting the mean response of all \dQuote{bbt} wells 
in a tray.  For type \dQuote{bbc}, a 4th degree polynomial is fit to the
\dQuote{bbc} responses and predicted values for each concentration are 
subtracted from treatment well responses. 
 
In many cases, a reasonable arrangement for a 96-well tray in a
cytotoxicity experiment is to use the last two rows of every
column as \dQuote{bbc} blanks, use the first six rows of the first
two columns as treatment controls, and use the first six rows of the
third to twelfth columns as treatment wells, with each column receiving
a different drug concentration.  

When conducting large simulations, it is recommended that many small
files be read from disk rather than one large file.  The time required to read
large files can be prohibitive.  
}
\value{
A list of class trayData with components:

    \item{concentrationResponse}{Data frame containing concentration-
        response information}
    \item{drugRatios}{Data frame containing summary information for each
        tray} 
}

\author{John Boik \email{john.boik@newearthbiomed.org}}

\seealso{
    \code{\link{prepareData}}
}
\examples{
dataFile = system.file(package="mixlow", "exdata/mixlowExampleData.txt")
excludeWells= list(row= c(1,2), col= c(1,1))
trayData <- readDataFile(filename=dataFile, excludeWells=excludeWells)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities } %% Fix me

