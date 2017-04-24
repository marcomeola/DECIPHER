\name{deltaSrules}
\alias{deltaSrules}
\docType{data}
\title{
Change in Entropy of Hybridization of Primer/Target Quadruplets in Solution
}
\description{
An 8D array with four adjacent base pairs of the primer and target sequences at a time.  Each dimension has five elements defining the residue at that position ("A", "C", "G", "T", or "-").  The array contains the standard entropy change of probe binding (dS, [kcal/mol]) for every quadruple base pairing.
}
\usage{data(deltaSrules)}
\format{
  The format is:
 num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.0226 0 0 0 0 ...
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
The first four dimensions correspond to the 4 primer positions from 5' to 3'.  The fifth to eighth dimensions correspond to the 4 positions from 5' to 3' of the target sequence.
}
\source{
Data from a variety of publications by SantaLucia \emph{et al}.
}
\references{
SantaLucia, J., Jr., & Hicks, D. (2004) The Thermodynamics of DNA Structural Motifs. Annual Review of Biophysics and Biomolecular Structure, 33(1), 415-440. doi:10.1146/annurev.biophys.32.110601.141800.
}
\examples{
data(deltaSrules)
# dS of primer = AGCT / target = A-CT pairing
deltaSrules["A", "G", "C", "T", "A", "-", "C", "T"]
}
\keyword{datasets}
