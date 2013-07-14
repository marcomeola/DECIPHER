\name{AlignSeqs}
\alias{AlignSeqs}
\title{
Aligns A Set of Unaligned Sequences
}
\description{
Performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.
}
\usage{
AlignSeqs(myDNAStringSet,
          perfectMatch = 2,
          misMatch = -3,
          gapOpening = -6,
          gapExtension = -4,
          terminalGap = -1,
          doNotAlign = 0.75,
          guideTree = NULL,
          orient = TRUE,
          verbose = TRUE)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object of unaligned sequences to be aligned.
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
  \item{doNotAlign}{
Numeric specifying the threshold above which to not attempt to align sequences.  Roughly equivalent to the fraction dissimilarity between unalignable sequences.  (See details section below.)
}
  \item{guideTree}{
Either \code{NULL} or a \code{data.frame} giving the ordered tree structure in which to align profiles.  If \code{NULL} then a guide tree will be constructed.  (See details section below.)
}
  \item{orient}{
Logical specifying whether some sequences may need to be reoriented (reverse complemented) before alignment.  If \code{TRUE}, an attempt to determine orientation will be performed with sequences reoriented as necessary.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The profile-to-profile method aligns a sequence set by merging profiles along a guide tree until all sequences are aligned.  If \code{guideTree=NULL}, an initial UPGMA guide tree is constructed based on a distance matrix of shared k-mers.  A second guide tree is built based on the initial alignment, and the alignment is refined using this tree.  If a \code{guideTree} is provided then sequences are only aligned once.  The \code{guideTree} should be provided in the output given by \code{IdClusters} with increasing levels of \code{cutoff}.

If \code{doNotAlign} is less than the height of the \code{guideTree}, then sequences more than \code{doNotAlign} distance apart are not aligned in order to save time.  For example, setting \code{doNotAlign} to \code{0.75} (the default) would not align groups of sequences that are more than 75\% dissimilar, whereas setting it to 1.0 would attempt to align all sequences.  The \code{doNotAlign} argument can be used to discourage attempts at aligning completely unrelated DNA sequences, such as random DNA which would be 75\% dissimilar on average.  In such cases, groups of sequences are not aligned to each other, and instead have gaps added to make them the same length.  In most reasonable alignment scenarios the sequences being aligned are generally not more than 50\% dissimilar.
}
\value{
A \code{DNAStringSet} of aligned sequences.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AlignProfiles}}, \code{\link{IdClusters}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db, remove="all")
alignedDNA <- AlignSeqs(dna)
BrowseSequences(alignedDNA, highlight=1)
}
