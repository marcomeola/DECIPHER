\name{HEC_MI}
\alias{HEC_MI1}
\alias{HEC_MI2}
\docType{data}
\title{
Mutual Information for Protein Secondary Structure Prediction
}
\description{
Arrays containing values of mutual information for single residues (\code{HEC_MI1}) and pairs of residues (\code{HEC_MI2}) located within 10 residues of the position being predicted (position "0").  The arrays have dimensions corresponding to the 20 (standard) amino acids, positions (-10 to 10), and states (helix ("H"), sheet ("E"), or coil ("C")).
}
\usage{
data("HEC_MI1")
data("HEC_MI2")
}
\format{
  The format of HEC_MI1 is:
 num [1:20, 1:21, 1:3] 0.04264 -0.00117 0.02641 0.08264 -0.04876 ...
 - attr(*, "dimnames")=List of 3
  ..$ : chr [1:20] "A" "R" "N" "D" ...
  ..$ : chr [1:21] "-10" "-9" "-8" "-7" ...
  ..$ : chr [1:3] "H" "E" "C"

  The format of HEC_MI2 is:
 num [1:20, 1:20, 1:21, 1:21, 1:3] 2.56 -Inf -Inf -Inf -Inf ...
 - attr(*, "dimnames")=List of 5
  ..$ : chr [1:20] "A" "R" "N" "D" ...
  ..$ : chr [1:20] "A" "R" "N" "D" ...
  ..$ : chr [1:21] "-10" "-9" "-8" "-7" ...
  ..$ : chr [1:21] "-10" "-9" "-8" "-7" ...
  ..$ : chr [1:3] "H" "E" "C"
}
\details{
The values in each matrix were derived based on a set of 15,201 proteins in the ASTRAL Compendium (Chandonia, 2004).  The 8-states assigned by the Dictionary of Protein Secondary Structure (DSSP) were reduced to 3-states via H = G, H, or I;  E = E;  and C = B, S, C, or T.
}
\references{
Chandonia, J. M. (2004). The ASTRAL Compendium in 2004. \emph{Nucleic Acids Research}, \bold{32(90001)}, 189D-192. doi:10.1093/nar/gkh034.
}
\examples{
data(HEC_MI1)
# the contribution of an arginine ("R")
# located 3 residues left of center
# to a helical ("H") state at the center
HEC_MI1["R", "-3", "H"]

data(HEC_MI2)
# the contribution of arginine and lysine ("K")
# located at positions -1 and +1, respectively
# to a coil ("C") state at the center position
HEC_MI2["R", "K", "-1", "1", "C"]

matplot(-10:10, t(HEC_MI1[,, "H"]),
       type="l", col=1:8, lty=rep(1:3, each=8),
       xlab="Amino Acid Position Relative to Center",
       ylab="Log-Odds of Helix at Center Position")
legend("bottomleft",
       lwd=1, col=1:8, lty=rep(1:3, each=8),
       legend=dimnames(HEC_MI1)[[1]], ncol=2)
}
\keyword{datasets}
