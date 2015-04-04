\name{StaggerAlignment}
\alias{StaggerAlignment}
\title{
Produce a Staggered Alignment
}
\description{
Staggers overlapping characters in a multiple sequence alignment that are better explained by multiple insertions than multiple deletions.
}
\usage{
StaggerAlignment(myXStringSet,
                 tree = NULL,
                 threshold = 3,
                 fullLength = FALSE,
                 processors = NULL,
                 verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of aligned sequences.
}
  \item{tree}{
A bifurcating \code{dendrogram} representing the evolutionary relationships between sequences, such as that created by \code{\link{IdClusters}}.  The root should be the topmost node of the \code{tree}.
}
  \item{threshold}{
Numeric giving the ratio of insertions to deletions that must be met to stagger a region of the alignment.  Specifically, the number of insertions divided by deletions must be less than threshold to stagger.
}
  \item{fullLength}{
Logical specifying whether the sequences are full-length (\code{TRUE}), or terminal gaps should be treated as missing data (\code{FALSE}, the default).  Either a single logical, a vector with one logical per sequence, or a list with \code{right} and \code{left} components containing logicals for the right and left sides of the alignment.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Multiple sequence aligners typically maximize true homologies at the expense of increased false homologies.  \code{StaggerAlignment} creates a ``staggered alignment'' which separates regions of the alignment that are likely not homologous into separate regions.  This re-balances the trade-off between true positives and false positives by decreasing the number of false homologies at the loss of some true homologies.  The resulting alignment is less aesthetically pleasing because it is widened by the introduction of many gaps.  However, in an evolutionary sense a staggered alignment is more correct because each aligned position represents a hypothesis about evolutionary events.  Overlapping characters in the alignment represent original positions in the ancestral (root) sequence that may have evolved through mutation, whereas non-overlapping regions are the result of insertion or deletion events.

The single parameter \code{threshold} controls the degree of staggering.  Its value represents the ratio of insertions to deletions that must be crossed in order to stagger a region.  A \code{threshold} of \code{1} would mean that any region which could be better explained by separate insertions than deletions should be staggered.  A higher value \code{threshold} makes it more likely to stagger, and vise-versa.  A very high value would conservatively stagger most regions, resulting in few false homologies but also few true homologies.  The default value (\code{3}) is intended to remove more false homologies than it eliminates in true homologies.  It may be preferable to tailor the \code{threshold} depending on the purpose of the alignment, as some downstream procedures (such as tree building) may be more or less sensitive to false homologies.
}
\value{
An \code{XStringSet} of aligned sequences.
}
\references{
Coming soon!
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AdjustAlignment}}, \code{\link{AlignSeqs}}, \code{\link{IdClusters}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db, remove="all")
alignedDNA <- AlignSeqs(dna)
staggerDNA <- StaggerAlignment(alignedDNA)
BrowseSeqs(staggerDNA, highlight=1)
}