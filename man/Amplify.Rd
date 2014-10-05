\name{Amplify}
\alias{Amplify}
\title{
Simulate Amplification of DNA by PCR
}
\description{
Predicts the amplification efficiency of theoretical PCR products (amplicons) given one or more primer sequences.
}
\usage{
Amplify(primers,
        myDNAStringSet,
        annealingTemp,
        P,
        ions,
        maxProductSize,
        includePrimers=TRUE,
        minEfficiency = 0.001,
        ...)
}
\arguments{
  \item{primers}{
A \code{DNAStringSet} object or character vector with one or more unaligned primer sequences in 5' to 3' orientation.
}
  \item{myDNAStringSet}{
A \code{DNAStringSet} object or character vector with unaligned template DNA sequences in 5' to 3' orientation.
}
  \item{annealingTemp}{
Numeric specifying the annealing temperature used in the PCR reaction.
}
  \item{P}{
Numeric giving the molar concentration of primers in the reaction.
}
  \item{ions}{
Numeric giving the molar sodium equivalent ionic concentration.  Values may range between 0.01M and 1M.
}
  \item{maxProductSize}{
Integer specifying the maximum length of PCR products (amplicons) in nucleotides.
}
  \item{includePrimers}{
Logical indicating whether to include the primer sequences in the theoretical PCR products.  (See details section below.)
}
  \item{minEfficiency}{
Numeric giving the minimum amplification efficiency of PCR products to include in the output (default \code{0.1\%}).  (See details section below.)
}
  \item{\dots}{
Additional arguments to be passed directly to \code{\link{CalculateEfficiencyPCR}}, including \code{batchSize}, \code{taqEfficiency}, \code{maxDistance}, \code{maxGaps}, and \code{processors}.
}
}
\details{
Exponential amplification in PCR requires the annealing and elongation of two primers from target sites on opposing strands of the template DNA.  If the template DNA sequence (e.g., chromosome) is known then predictions of theoretical amplicons can be obtained from in silico simulations of amplification.  \code{Amplify} first searches for primer target sites on the template DNA, and then calculates an amplification efficiency from each target site using \code{\link{CalculateEfficiencyPCR}}.

If \code{taqEfficiency} is \code{TRUE} (the default), the amplification efficiency of each primer is defined as the product of hybridization efficiency and elongation efficiency.  Amplification efficiency must be at least \code{minEfficiency} for a primer to be amplified in silico.  Overall amplification efficiency of the PCR product is then calculated as the geometric mean of the two (i.e., forward and reverse) primers' efficiencies.  Finally, amplicons are generated if the two primers are within \code{maxProductSize} nucleotides downstream of each other.

Potential PCR products are returned, either with or without including the primer sequences in the amplicon.  The default (\code{includePrimers=TRUE}) is to incorporate the primer sequences as would normal occur during amplification.  The alternative is to return the complete template sequence including the target sites, which may not exactly match the primer sequences.
}
\value{
A \code{DNAStringSet} object containing potential PCR products, sorted by highest-to-lowest amplification efficiency.  The sequences are named by their predicted amplification efficiency followed by the index of each primer in \code{myDNAStringSet}.  (See examples section below.)
}
\references{
ES Wright et al. (2013) "Exploiting Extension Bias in PCR to Improve Primer Specificity in Ensembles of Nearly Identical DNA Templates." Environmental Microbiology, doi:10.1111/1462-2920.12259.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\note{
The program OligoArrayAux (\url{http://mfold.rna.albany.edu/?q=DINAMelt/OligoArrayAux}) must be installed in a location accessible by the system.  For example, the following code should print the installed OligoArrayAux version when executed from the R console:

\code{system("hybrid-min -V")}
}
\seealso{
\code{\link{CalculateEfficiencyPCR}}, \code{\link{DesignPrimers}}
}
\examples{
data(yeastSEQCHR1)

# not run (must have OligoArrayAux installed first):

# match a single primer that acts as both the forward and reverse
primer1 <- "TGGAAGCTGAAACG"
\dontrun{Amplify(primer1, yeastSEQCHR1, annealingTemp=55, P=4e-7, ions=0.225, maxProductSize=500)}

# perform a typical amplification with two primer sequences:
primer2 <- c("GGCTGTTGTTGGTGTT", "TGTCATCAGAACACCAA")
\dontrun{Amplify(primer2, yeastSEQCHR1, annealingTemp=55, P=4e-7, ions=0.225, maxProductSize=500)}

# perform a multiplex PCR amplification with multiple primers:
primers <- c(primer1, primer2)
\dontrun{Amplify(primers, yeastSEQCHR1, annealingTemp=55, P=4e-7, ions=0.225, maxProductSize=500)}
}