\name{AlignProfiles}
\alias{AlignProfiles}
\title{
Aligns Two Sets of Aligned Sequences
}
\description{
Aligns two sets of one or more aligned sequences by first generating representative profiles, then aligning the profiles with dynamic programming, and finally merging the two aligned sequence sets.
}
\usage{
AlignProfiles(pattern,
              subject,
              p.weight = 1,
              s.weight = 1,
              perfectMatch = 2,
              misMatch = -3,
              gapOpening = -6,
              gapExtension = -4,
              terminalGap = -1)
}
\arguments{
  \item{pattern}{
A \code{DNAStringSet} object of aligned sequences to use as the pattern.
}
  \item{subject}{
A \code{DNAStringSet} object of aligned sequences to use as the subject.
}
  \item{p.weight}{
A numeric vector of weights for each sequence in the pattern to use in generating a profile, or a single number implying equal weights.
}
  \item{s.weight}{
A numeric vector of weights for each sequence in the subject to use in generating a profile, or a single number implying equal weights.
}
  \item{perfectMatch}{
Numeric giving the reward for aligning two matching nucleotides in the alignment.
}
  \item{misMatch}{
Numeric giving the cost for aligning two mismatched nucleotides in the alignment.
}
  \item{gapOpening}{
Numeric giving the cost for opening a gap in the alignment.
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap in the alignment.
}
  \item{terminalGap}{
Numeric giving the cost for allowing leading and trailing gaps in the alignment.  Either two numbers, the first for leading gaps and the second for trailing gaps, or a single number for both.
}
}
\details{
Profiles are aligned using dynamic programming, a variation of the Needleman-Wunsch algorithm.
}
\value{
A \code{DNAStringSet} of aligned sequences.
}
\references{
Needleman S., Wunsch, C. (1970).  A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology, \bold{48(3)}, 443-453.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AlignSeqs}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna1 <- SearchDB(db, remove="common", limit=100) # the first 100 sequences
dna2 <- SearchDB(db, remove="common", limit="100,100") # the rest
alignedDNA <- AlignProfiles(dna1, dna2)
BrowseSequences(alignedDNA, highlight=1)
}
