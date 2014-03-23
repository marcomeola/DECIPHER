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
              perfectMatch = 6,
              misMatch = 0,
              gapOpening = -9,
              gapExtension = -3,
              terminalGap = -2,
              restrict = -1000,
              anchor = 0.7,
              substitutionMatrix = "BLOSUM62",
              processors = NULL)
}
\arguments{
  \item{pattern}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of aligned sequences to use as the pattern.
}
  \item{subject}{
A \code{XStringSet} object of aligned sequences to use as the subject.  Must match the type of the pattern.
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
  \item{restrict}{
Numeric specifying the lowest relative score to consider when aligning.  The default (\code{-1000}) will align most inputs that can reasonably be globally aligned without any loss in accuracy.  Input sequences with high similarity could be more restricted (e.g., \code{-500}), whereas a \code{pattern} and \code{subject} with little overlap should be less restricted (e.g., \code{-10000}).  (See details section below.)
}
  \item{anchor}{
Numeric giving the fraction of sequences with identical k-mers required to become an anchor point, or \code{NA} to not use anchors.  (See details section below.)
}
  \item{substitutionMatrix}{
Name of the amino acid substitution matrix to use in alignment.  Must be one of:  ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250''.  Only applicable for \code{AAStringSet} inputs.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
}
\details{
Profiles are aligned using dynamic programming, a variation of the Needleman-Wunsch algorithm for global alignment.  This method works by filling in a matrix of the possible ``alignment space'' by considering all matches, insertions, and deletions between two sequence profiles.  The highest scoring alignment is then used to add gaps to each of the input sequence sets.  The default input parameters have been optimized to align amino acid and RNA sequences using the structural benchmarks.  Careful consideration should be given to the choice of \code{substitutionMatrix} (for \code{AAStringSet} inputs) and sequence weightings so that they reflect the divergence between input sequences.

The dynamic programming method requires order \code{N*M} time and memory space where \code{N} and \code{M} are the width of the pattern and subject.  Therefore heuristics can be useful to improve performance on long input sequences.  The \code{restrict} parameter can be used to dynamically constrain the possible ``alignment space'' to only paths that will likely include the final alignment, which in the best case can improve the speed from quadratic time to linear time.  The degree of restriction is critically important, and if the sequences are not necessarily mostly overlapping then \code{restrict} should be relaxed (lower than the default).  For example, if aligning a short profile to a long profile then restrict should be set to \code{-Inf}.

The argument \code{anchor} can be used to split the global alignment into multiple sub-alignments.  This can greatly decrease the memory requirement for long sequences when appropriate anchor points can be found.  Anchors are 15-mer (for DNA/RNA) or 7-mer (for AA) subsequences that are shared between between at least \code{anchor} fraction of \code{pattern}(s) and \code{subject}(s).  Anchored ranges are extended along the length of each sequence in a manner designed to split the alignment into sub-alignments that can be separately solved.  For most input sequences \code{anchor} has little or no effect on accuracy, but anchoring can be disabled by setting \code{anchor=NA}.
}
\value{
An \code{XStringSet} of aligned sequences.
}
\references{
Needleman S., Wunsch, C. (1970).  A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology, \bold{48(3)}, 443-453.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AlignDB}}, \code{\link{AlignSeqs}}, \code{\link{AlignTranslation}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna1 <- SearchDB(db, remove="common", limit=100) # the first 100 sequences
dna2 <- SearchDB(db, remove="common", limit="100,100") # the rest
alignedDNA <- AlignProfiles(dna1, dna2)
BrowseSequences(alignedDNA, highlight=1)
}
