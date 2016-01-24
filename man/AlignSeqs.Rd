\name{AlignSeqs}
\alias{AlignSeqs}
\title{
Align A Set of Unaligned Sequences
}
\description{
Performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.
}
\usage{
AlignSeqs(myXStringSet,
         guideTree = NULL,
         iterations = 1,
         refinements = 1,
         gapOpening=c(-16, -12),
         gapExtension=c(-2, -1),
         structures = NULL,
         FUN = AdjustAlignment,
         levels = c(0.95, 0.7, 10, 5),
         processors = 1,
         verbose = TRUE,
         \dots)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of unaligned sequences.
}
  \item{guideTree}{
Either \code{NULL} or a \code{data.frame} giving the ordered tree structure in which to align profiles.  If \code{NULL} then a guide tree will be automatically constructed based on the order of shared k-mers.
}
  \item{iterations}{
Number of iteration steps to perform.  During each iteration step the guide tree is regenerated based on the alignment and the sequences are realigned.
}
  \item{refinements}{
Number of refinement steps to perform.  During each refinement step groups of sequences are realigned to rest of the sequences, and the best of these two alignments (before and after realignment) is kept.
}
  \item{gapOpening}{
Single numeric giving the cost for opening a gap in the alignment, or two numbers giving the minimum and maximum costs.  In the latter case the cost will be varied depending upon whether the groups of sequences being aligned are nearly identical or maximally distant.
}
  \item{gapExtension}{
Single numeric giving the cost for extending an open gap in the alignment, or two numbers giving the minimum and maximum costs.  In the latter case the cost will be varied depending upon whether the groups of sequences being aligned are nearly identical or maximally distant.
}
  \item{structures}{
Either a list of secondary structure probabilities matching the \code{structureMatrix}, such as that output by \code{PredictHEC}, or \code{NULL} to generate the structures automatically.  Only applicable if \code{myXStringSet} is an \code{AAStringSet}.
}
  \item{FUN}{
A function to be applied after each profile-to-profile alignment.  (See details section below.)
}
  \item{levels}{
Numeric with four elements specifying the levels above which to apply \code{FUN}.  (See details section below.)
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{AlignProfiles}}, including \code{perfectMatch}, \code{misMatch}, \code{gapPower}, \code{terminalGap}, \code{restrict}, \code{anchor}, \code{normPower}, \code{substitutionMatrix}, and \code{structureMatrix}.
}
}
\details{
The profile-to-profile method aligns a sequence set by merging profiles along a guide tree until all the input sequences are aligned.  This process has three main steps:  (1)  If \code{guideTree=NULL}, an initial single-linkage guide tree is constructed based on a distance matrix of shared k-mers.  If an initial \code{guideTree} is provided then  the \code{guideTree} should be provided in the output given by \code{IdClusters} with ascending levels of \code{cutoff}.  (2)  If \code{iterations} is greater than zero, then a UPGMA guide tree is built based on the initial alignment and the sequences are re-aligned along this tree.  This process repeated \code{iterations} times or until convergence.  (3)  If \code{refinements} is greater than zero, then groups of sequences are iteratively realigned to the full-alignment.  This process generates two alignments, the best of which is chosen based on its sum-of-pairs score.  This refinement process is repeated \code{refinements} times, or until no improvement can be made.

The \code{FUN} function is applied during each of the three steps based on \code{levels}.  The purpose of \code{levels} is to speed-up the alignment process by not running \code{FUN} on the alignment when it is unnecessary.  The default \code{levels} specify that \code{FUN} should be run on the sequences when the initial tree is above 0.95 average dissimilarity, when the iterative tree is above 0.7 average dissimilarity, and after every tenth improvement made during refinement.  The final element of levels prevents \code{FUN} from being applied at any point to less than 5 sequences.  The \code{FUN} function is always applied just before returning the alignment, independently of the first three values of \code{levels}.  The default \code{FUN} is \code{AdjustAlignment}, but \code{FUN} accepts any function that takes in an \code{XStringSet} as its first argument, and \code{weights}, \code{processors}, and \code{substitutionMatrix} as optional arguments.  For example, the default \code{FUN} could be altered to not perform any function by setting it equal to \code{FUN=function(x, ...) return(x)} where \code{x} is an \code{XStringSet}.
}
\value{
An \code{XStringSet} of aligned sequences.
}
\references{
ES Wright (2015) "DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment". BMC Bioinformatics, doi:10.1186/s12859-015-0749-z.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AdjustAlignment}}, \code{\link{AlignDB}}, \code{\link{AlignProfiles}}, \code{\link{AlignSynteny}}, \code{\link{AlignTranslation}}, \code{\link{IdClusters}}, \code{\link{StaggerAlignment}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db, remove="all")
alignedDNA <- AlignSeqs(dna)
BrowseSeqs(alignedDNA, highlight=1)
}
