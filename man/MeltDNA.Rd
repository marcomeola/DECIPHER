\name{MeltDNA}
\alias{MeltDNA}
\title{
Simulate Melting of DNA
}
\description{
The denaturation of double-stranded DNA occurs over a range of temperatures.  Beginning from a helical state, DNA will transition to a random-coil state as temperature is increased.  \code{MeltDNA} predicts the positional helicity, melt curve, or its negative derivate at different temperatures.
}
\usage{
MeltDNA(myDNAStringSet,
        type = "derivative",
        temps = 50:100,
        ions = 0.2)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object or character vector with one or more sequences in 5' to 3' orientation.
}
  \item{type}{
Character string indicating the type of results desired.  This should be (an abbreviation of) one of \code{"derivative curves"}, \code{"melt curves"}, or \code{"positional probabilities"}.
}
  \item{temps}{
Numeric vector of temperatures (in degrees Celsius).
}
  \item{ions}{
Numeric giving the molar sodium equivalent ionic concentration.  Values must be at least 0.01M.
}
}
\details{
When designing a high resolution melt (HRM) assay, it is useful to be able to predict the results before performing the experiment.  Multi-state models of DNA melting can provide near-qualitative agreement with experimental DNA melt curves obtained with quantitative PCR (qPCR).  \code{MeltDNA} employs the algorithm of Tostesen et al. (2003) with an approximation for loop entropy that runs in nearly linear time and memory, which allows very long DNA sequences (up to 100,000 base pairs) to be analyzed.

Denaturation is a highly cooperative process whereby regions of double-stranded DNA tend to melt together.  For short sequences (< 100 base pairs) there is typically a single transition from a helical to random-coil state.  Longer sequences may exhibit more complex melting behavior with multiple peaks, as domains of the DNA melt at different temperatures.  The melting curve represents the average fractional helicity (Theta) at each temperature, and can be used for genotyping with high resolution melt analysis.
}
\value{
\code{MeltDNA} can return three \code{type}s of results:  positional helicity, melting curves, or the negative derivative of the melting curves.  If \code{type} is \code{"position"}, then a list is returned with one component for each sequence in \code{myDNAStringSet}.  Each list component contains a matrix with the probability of helicity (Theta) at each temperature (rows) and every position in the sequence (columns).

If \code{type} is \code{"melt"}, then a matrix with the average Theta across the entire sequence is returned.  This matrix has a row for each input temperature (\code{temps}), and a column for each sequence in \code{myDNAStringSet}.  For example, the value in element \code{[3, 4]} is the average helicity of the fourth input sequence at the third input temperature.  If \code{type} is \code{"derivative"} then the values in the matrix are the derivative of the melt curve at each temperature.
}
\references{
SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proceedings of the National Academy of Sciences, 95(4), 1460-1465.

Tostesen, E., et al. (2003). Speed-up of DNA melting algorithm with complete nearest neighbor properties. Biopolymers, 70(3), 364-376. doi:10.1002/bip.10495.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\note{
\code{MeltDNA} uses nearest neighbor parameters from SantaLucia (1998).
}
\seealso{
\code{\link{AmplifyDNA}}, \code{\link{CalculateEfficiencyPCR}}, \code{\link{DesignSignatures}}
}
\examples{
fas <- system.file("extdata", "IDH2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

# plot the melt curve for the two alleles
temps <- seq(85, 100, 0.2)
m <- MeltDNA(dna,
             type="melt", temps=temps, ions=0.1)
matplot(temps, m,
        type="l", xlab="Temperature (\u00B0C)", ylab="Average Theta")
legend("topright", names(dna), lty=seq_along(dna), col=seq_along(dna))

# plot the negative derivative curve for a subsequence of the two alleles
temps <- seq(80, 95, 0.25)
m <- MeltDNA(subseq(dna, 492, 542),
             type="derivative", temps=temps)
matplot(temps, m,
        type="l", xlab="Temperature (\u00B0C)", ylab="-d(Theta)/dTemp")
legend("topright", names(dna), lty=seq_along(dna), col=seq_along(dna))

# plot the positional helicity profile for the IDH2 allele
temps <- seq(90.1, 90.5, 0.1)
m <- MeltDNA(dna[1],
             type="position", temps=temps, ions=0.1)
matplot(seq_len(dim(m[[1]])[2]), t(m[[1]]),
        type="l", xlab="Nucleotide Position", ylab="Theta")
temps <- formatC(temps, digits=1, format="f")
legend("topright", legend=paste(temps, "\u00B0C", sep=""),
        col=seq_along(temps), lty=seq_along(temps), bg="white")
}