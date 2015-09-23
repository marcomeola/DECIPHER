\name{MIQS}
\alias{MIQS}
\docType{data}
\title{
MIQS Amino Acid Substitution Matrix
}
\description{
The MIQS amino acid substitution matrix defined by Yamada & Tomii (2014).
}
\usage{data("MIQS")}
\format{
  The format is:
 num [1:25, 1:25] 3.2 -1.3 -0.4 -0.4 1.5 -0.2 -0.4 0.4 -1.2 -1.3 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:25] "A" "R" "N" "D" ...
  ..$ : chr [1:25] "A" "R" "N" "D" ...
}
\details{
Substitution matrix values represent the log-odds of observing an aligned pair of amino acids versus the likelihood of finding the pair by chance.  Values in the MIQS matrix are in units of third-bits (\eqn{log(odds\ ratio)*3/log(2)}).
}
\source{
Yamada, K., & Tomii, K. (2014). Revisiting amino acid substitution matrices for identifying distantly related proteins. \emph{Bioinformatics}, \bold{30(3)}, 317-325. doi:10.1093/bioinformatics/btt694.
}
\examples{
data(MIQS)
MIQS["A", "R"] # score for A/R pairing

data(BLOSUM62)
plot(BLOSUM62[1:20, 1:20], MIQS[1:20, 1:20])
abline(a=0, b=1)
}
\keyword{datasets}
