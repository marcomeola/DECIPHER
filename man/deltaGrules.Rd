\name{deltaGrules}
\alias{deltaGrules}
\docType{data}
\title{
An 8-dimensional array containing the free energy of hybridization of probe/target quadruplets.
}
\description{
The 8D array works with four adjacent base pairs of the probe and target sequence at a time.  Each dimension has five elements defining the residue at that position ("A", "C", "G", "T", or "-").  The array contains the standard Gibb's free energy change of probe binding (dG, [kcal/mol]) for every quadruple base pairing.
}
\usage{data(deltaGrules)}
\format{
  The format is:
 num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.141 0 0 0 0 ...
 - attr(*, "dimnames")=List of 8
  ..$ : chr [1:5] "A" "C" "G" "T" ...
  ..$ : chr [1:5] "A" "C" "G" "T" ...
  ..$ : chr [1:5] "A" "C" "G" "T" ...
  ..$ : chr [1:5] "A" "C" "G" "T" ...
  ..$ : chr [1:5] "A" "C" "G" "T" ...
  ..$ : chr [1:5] "A" "C" "G" "T" ...
  ..$ : chr [1:5] "A" "C" "G" "T" ...
  ..$ : chr [1:5] "A" "C" "G" "T" ...
}
\details{
The first four dimensions correspond to the 4 probe positions from 5' to 3'.  The fifth to eighth dimensions correspond to the 4 positions from 5' to 3' of the target sequence.
}
\source{
Data obtained using NimbleGen microarrays and a Linear Free Energy Model developed by Yilmaz \emph{et al}.
}
\references{
Coming soon!
}
\examples{
data(deltaGrules)
# dG of probe = AGCT / target = A-CT pairing
deltaGrules["A", "G", "C", "T", "A", "-", "C", "T"]
}
\keyword{datasets}
