\name{RESTRICTION_ENZYMES}
\alias{RESTRICTION_ENZYMES}
\docType{data}
\title{
Common Restriction Enzyme's Cut Sites
}
\description{
A character vector of common restriction sites named by the restriction enzyme that cuts at each site.  Sequence specificity is listed in 5' to 3' orientation based on the \code{IUPAC_CODE_MAP}.  The cut site is either signified by a ``/'' for palindromic sites, or two numbers giving the position of the top and bottom cut positions relative to the site's 3'-end.
}
\usage{data(RESTRICTION_ENZYMES)}
\format{
  The format is:
 Named chr [1:224] "GACGT/C" "G/GTACC" "GT/MKAC" ...
 - attr(*, "names")= chr [1:224] "AatII" "Acc65I" "AccI" "AciI" ...
}
\source{
Restriction enzymes sold by \href{http://www.neb.com}{New England BioLabs}.
}
\examples{
data(RESTRICTION_ENZYMES)
}
\keyword{datasets}
