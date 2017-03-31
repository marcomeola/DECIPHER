\name{AdjustAlignment}
\alias{AdjustAlignment}
\title{
Improve An Existing Alignment By Adjusting Gap Placements
}
\description{
Makes small adjustments by shifting groups of gaps left and right to find their optimal positioning in a multiple sequence alignment.
}
\usage{
AdjustAlignment(myXStringSet,
                perfectMatch = 5,
                misMatch = 0,
                gapLetter = -3,
                gapOpening = -0.1,
                gapExtension = 0,
                substitutionMatrix = NULL,
                shiftPenalty = -0.2,
                threshold = 0.1,
                weight = 1,
                processors = 1)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of aligned sequences.
}
  \item{perfectMatch}{
Numeric giving the reward for aligning two matching nucleotides in the alignment.  Only used for \code{DNAStringSet} or \code{RNAStringSet} inputs.
}
  \item{misMatch}{
Numeric giving the cost for aligning two mismatched nucleotides in the alignment.  Only used for \code{DNAStringSet} or \code{RNAStringSet} inputs.
}
  \item{gapLetter}{
Numeric giving the cost for aligning gaps to letters.  A lower value (more negative) encourages the overlapping of gaps across different sequences in the alignment.
}
  \item{gapOpening}{
Numeric giving the cost for opening or closing a gap in the alignment.
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap in the alignment.
}
  \item{substitutionMatrix}{
Either a substitution matrix representing the substitution scores for an alignment or the name of the amino acid substitution matrix to use in alignment.  The latter may be one of the following:  ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250'', or ``MIQS''.  The default (NULL) will use the \code{perfectMatch} and \code{misMatch} penalties for DNA/RNA or ``MIQS'' for AA.  (See examples section below.)
}
  \item{shiftPenalty}{
Numeric giving the cost for every additional position that a group of gaps is shifted.
}
  \item{threshold}{
Numeric specifying the improvement in score required to permanently apply an adjustment to the alignment.
}
  \item{weight}{
A numeric vector of weights for each sequence, or a single number implying equal weights.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
}
\details{
The process of multiple sequence alignment often results in the integration of small imperfections into the final alignment.  Some of these errors are obvious by-eye, which encourages manual refinement of automatically generated alignments.  However, the manual refinement process is inherently subjective and time consuming.  \code{AdjustAlignment} refines an existing alignment in a process similar to that which might be applied manually, but in a repeatable and must faster fashion.  This function shifts all of the gaps in an alignment to the left and right to find their optimal positioning.  The optimal position is defined as the position that maximizes the alignment ``score'', which is determined by the input parameters.  The resulting alignment will be similar to the input alignment but with many imperfections eliminated.  Note that the affine gap penalties here are different from the more flexible penalties used in \code{\link{AlignProfiles}}, and have been optimized independently.
}
\value{
An \code{XStringSet} of aligned sequences.
}
\references{
Wright, E. S. (2015). DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC Bioinformatics, 16, 322. http://doi.org/10.1186/s12859-015-0749-z
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignSeqs}}, \code{\link{AlignTranslation}}, \code{\link{StaggerAlignment}}
}
\examples{
# a trivial example
aa <- AAStringSet(c("ARN-PK", "ARRP-K"))
aa
AdjustAlignment(aa)

# a real example
fas <- system.file("extdata", "Streptomyces_ITS_aligned.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
adjustedDNA <- AdjustAlignment(dna)
BrowseSeqs(adjustedDNA, highlight=1)
adjustedDNA==dna # most sequences were adjusted
}