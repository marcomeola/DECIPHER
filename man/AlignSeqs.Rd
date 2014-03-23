\name{AlignSeqs}
\alias{AlignSeqs}
\title{
Aligns A Set of Unaligned Sequences
}
\description{
Performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.
}
\usage{
AlignSeqs(myXStringSet,
          guideTree = NULL,
          orient = FALSE,
          processors = NULL,
          verbose = TRUE,
          ...)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of unaligned sequences.
}
  \item{guideTree}{
Either \code{NULL} or a \code{data.frame} giving the ordered tree structure in which to align profiles.  If \code{NULL} then a guide tree will be constructed.
}
  \item{orient}{
Logical specifying whether some sequences may need to be reoriented before alignment.  If \code{TRUE}, an attempt to determine the best orientation (reverse and/or complement) will be performed with sequences reoriented as necessary to match the orientation of the longest sequence.  Not applicable for an \code{AAStringSet} input.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{AlignProfiles}}, including \code{perfectMatch}, \code{misMatch}, \code{gapOpening}, \code{gapExtension}, \code{terminalGap}, \code{restrict}, and \code{anchor}.
}
}
\details{
The profile-to-profile method aligns a sequence set by merging profiles along a guide tree until all sequences are aligned.  If \code{guideTree=NULL}, an initial UPGMA guide tree is constructed based on a distance matrix of shared k-mers.  A second guide tree is built based on the initial alignment, and the alignment is refined using this tree.  If a \code{guideTree} is provided then sequences are only aligned once.  The \code{guideTree} should be provided in the output given by \code{IdClusters} with ascending levels of \code{cutoff}.
}
\value{
An \code{XStringSet} of aligned sequences.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AlignDB}}, \code{\link{AlignProfiles}}, \code{\link{AlignTranslation}}, \code{\link{IdClusters}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db, limit=10, remove="all")
alignedDNA <- AlignSeqs(dna)
BrowseSequences(alignedDNA, highlight=1)
}
