\name{DesignArray}
\alias{DesignArray}
\title{
Designs a set of DNA Microarray Probes for Detecting Sequences
}
\description{
Chooses the set of microarray probes maximizing sensitivity and specificity to each target consensus sequence.
}
\usage{
DesignArray(myDNAStringSet,
            maxProbeLength = 24,
            minProbeLength = 20,
            maxPermutations = 4,
            numRecordedMismatches = 500,
            numProbes = 10,
            start = 1,
            end = NULL,
            maxOverlap = 5,
            hybridizationFormamide = 10,
            minMeltingFormamide = 15,
            maxMeltingFormamide = 20,
            minScore = -1e+12,
            processors = NULL,
            verbose = TRUE)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object of aligned consensus sequences.
}
  \item{maxProbeLength}{
The maximum length of probes, not including the poly-T spacer.  Ideally less than 27 nucleotides.
}
  \item{minProbeLength}{
The minimum length of probes, not including the poly-T spacer.  Ideally more than 18 nucleotides.
}
  \item{maxPermutations}{
The maximum number of probe permutations required to represent a target site.  For example, if a target site has an 'N' then 4 probes are required because probes cannot be ambiguous.  Typically fewer permutations are preferably because this requires less space on the microarray and simplifies interpretation of the results.
}
  \item{numRecordedMismatches}{
The maximum number of recorded potential cross-hybridizations for any target site.
}
  \item{numProbes}{
The target number of probes on the microarray per input consensus sequence.
}
  \item{start}{
Integer specifying the starting position in the alignment where potential forward primer target sites begin.  Preferably a position that is included in most sequences in the alignment.
}
  \item{end}{
Integer specifying the ending position in the alignment where potential reverse primer target sites end.  Preferably a position that is included in most sequences in the alignment.
}
  \item{maxOverlap}{
Maximum overlap in nucleotides between target sites on the sequence.
}
  \item{hybridizationFormamide}{
The formamide concentration (\%, vol/vol) used in hybridization at 42 degrees Celsius.  Note that this concentration is used to approximate hybridization efficiency of cross-amplifications.
}
  \item{minMeltingFormamide}{
The minimum melting point formamide concentration (\%, vol/vol) of the designed probes.  The melting point is defined as the concentration where half of the template is bound to probe.
}
  \item{maxMeltingFormamide}{
The maximum melting point formamide concentration (\%, vol/vol) of the designed probes.  Must be greater than the \code{minMeltingFormamide}.
}
  \item{minScore}{
The minimum score of designed probes before exclusion.  A greater \code{minScore} will accelerate the code because more target sites will be excluded from consideration.  However,  if the \code{minScore} is too high it will prevent any target sites from being recorded.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The algorithm begins by determining the optimal length of probes required to meet the input constraints while maximizing sensitivity to the target consensus sequence at the specified hybridization formamide concentration.  This set of potential target sites is then scored based on the possibility of cross-hybridizing to the other non-target sequences.  The set of probes is returned with the minimum possibility of cross-hybridizing is chosen to represent each consensus sequence.
}
\value{
A \code{data.frame} with the optimal set of probes matching the specified constraints.  Each row lists the probe's target sequence (\code{name}), \code{start} position, \code{length} in nucleotides, start and end position in the sequence alignment, number of \code{permutations}, \code{score}, melt point in percent \code{formamide} at 42 degrees Celsius, hybridization efficiency (\code{hyb_eff}), target site, and probe(s).  Probes are designed such that the stringency is determined by the equilibrium hybridization conditions and not subsequent washing steps.
}
\references{
ES Wright et al. (2013) Identification of Bacterial and Archaeal Communities From Source to Tap. Water Research Foundation, Denver, CO.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{Array2Matrix}}, \code{\link{NNLS}}
}
\examples{
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
names(dna) <- 1:length(dna)
probes <- DesignArray(dna)
probes[1,]
}
