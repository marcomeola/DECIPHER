\name{CalculateEfficiencyPCR}
\alias{CalculateEfficiencyPCR}
\title{
Predicts Amplification Efficiency of Primer Sequences
}
\description{
Calculates the amplification efficiency of primers from their hybridization efficiency and elongation efficiency at the target site.
}
\usage{
CalculateEfficiencyPCR(primer,
                       target,
                       temp,
                       P,
                       ions,
                       batchSize = 1000,
                       taqEfficiency = TRUE,
                       maxDistance = 0.4,
                       maxGaps = 2)
}
\arguments{
  \item{primer}{
A \code{DNAStringSet} object or character vector with unaligned primer sequences in 5' to 3' orientation.
}
  \item{target}{
A \code{DNAStringSet} object or character vector with unaligned target or non-target sequences in 5' to 3' orientation.
}
  \item{temp}{
Numeric specifying the annealing temperature used in the PCR reaction.
}
  \item{P}{
Numeric giving the molar concentration of primers in the reaction.
}
  \item{ions}{
Numeric giving the molar sodium equivalent ionic concentration.  Values may range between 0.01M and 1M.
}
  \item{batchSize}{
Integer specifying the number of primers to simulate hybridization per batch.  See the Description section below.
}
  \item{taqEfficiency}{
Logical determining whether to make use of elongation efficiency and maxDistance to increase predictive accuracy for \emph{Taq} DNA Polymerase amplifying primers with mismatches near the 3' terminus.  Note that this should be set to FALSE if using a high-fidelity polymerase with 3' to 5' exonuclease activity.
}
  \item{maxDistance}{
Numeric specifying the maximal fraction of mismatched base pairings on a rolling basis beginning from the 3' end of the primer.  Only used if \code{taqEfficiency} is \code{TRUE}.
}
  \item{maxGaps}{
Integer specifying the maximum number of insertions or deletions (indels) in the primer/target alignment.  Only used if \code{taqEfficiency} is \code{TRUE}.
}
}
\details{
Amplification of pairwise primer/target pairs is simulated \emph{in silico}.  A complex model of hybridization is employed that takes into account the side reactions resulting in probe-folding, target-folding, and primer-dimer formation.  The resulting hybridization efficiency is multiplied by the elongation efficiency to predict the overall efficiency of amplification.

Free energy is obtained from system calls to UNAFold, which must be properly installed (see the Notes section below).  Primer/target pairs are sent to UNAFold in batches of \code{batchSize}, which prevents systems calls from being too many characters.
}
\value{
A vector of predicted efficiencies for amplifying each primer/target pair of sequences.
}
\references{
Coming Soon!
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\note{
The program UNAFold (\url{http://mfold.rna.albany.edu/?q=DINAMelt/software}) must be installed in a location accessible by the system.  For example, the following code should print the installed UNAFold version when executed from the R console:  \code{system("hybrid -V")}.
}
\seealso{
\code{\link{DesignPrimers}}, \code{\link{TileSeqs}}
}
\examples{
primers <- c("AAAAACGGGGAGCGGGGGG", "AAAAACTCAACCCGAGGAGCGCGT")
targets <- reverseComplement(DNAStringSet(primers))
# not run (must have UNAFold installed first):
#CalculateEfficiencyPCR(primers, targets, temp=75, P=4e-7, ions=.225)
}