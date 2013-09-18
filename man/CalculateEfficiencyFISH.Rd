\name{CalculateEfficiencyFISH}
\alias{CalculateEfficiencyFISH}
\title{
Predicts Thermodynamic Parameters of Probe/Target Sequence Pairs
}
\description{
Calculates the Gibb's free energy, formamide melt point, and hybridization efficiency of probe/target (DNA/RNA) pairs.
}
\usage{
CalculateEfficiencyFISH(probe,
                        target,
                        temp,
                        P,
                        ions,
                        FA,
                        batchSize = 1000)
}
\arguments{
  \item{probe}{
A \code{DNAStringSet} object or character vector with unaligned probe sequences in 5' to 3' orientation.
}
  \item{target}{
A \code{DNAStringSet} object, \code{RNAStringSet}, or character vector with unaligned target or non-target sequences in 5' to 3' orientation.  The DNA base Thymine will be treated the same as Uracil.
}
  \item{temp}{
Numeric specifying the hybridization temperature, typically \code{46} degrees Celsius.
}
  \item{P}{
Numeric giving the molar concentration of probes during hybridization.
}
  \item{ions}{
Numeric giving the molar sodium equivalent ionic concentration.  Values may range between 0.01M and 1M.  Note that salt correction is not available for thermodynamic rules of RNA/RNA interactions, which were determined at \code{1} molar concentration.
}
  \item{FA}{
Numeric concentration (as percent v/v) of the denaturant formamide in the hybridization buffer.
}
  \item{batchSize}{
Integer specifying the number of probes to simulate hybridization per batch.  See the Description section below.
}
}
\details{
Hybridization efficiency of pairwise probe/target (DNA/RNA) pairs is simulated \emph{in silico}.  Gibb's free energies are obtained from system calls to OligoArrayAux, which must be properly installed (see the Notes section below).  Probe/target pairs are sent to OligoArrayAux in batches of \code{batchSize}, which prevents systems calls from being too many characters.
}
\value{
A matrix of predicted hybridization efficiency (\code{HybEff}), formamide melt point (\code{FAm}), and free energy (\code{ddG1} and \code{dG1}) for each probe/target pair.
}
\references{
Coming Soon!
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\note{
The program OligoArrayAux (\url{http://mfold.rna.albany.edu/?q=DINAMelt/OligoArrayAux}) must be installed in a location accessible by the system.  For example, the following code should print the installed OligoArrayAux version when executed from the R console:

\code{system("hybrid-min -V")}
}
\seealso{
\code{\link{DesignProbes}}, \code{\link{TileSeqs}}
}
\examples{
probe <- c("GGGCTTTCACATCAGACTTAAGAAACC", "CCCCACGCTTTCGCGCC")
target <- reverseComplement(DNAStringSet(probe))
# not run (must have OligoArrayAux installed first):
#CalculateEfficiencyFISH(probe, target, temp=46, P=250e-9, ions=1, FA=35)
}