\name{AlignProfiles}
\alias{AlignProfiles}
\title{
Align Two Sets of Aligned Sequences
}
\description{
Aligns two sets of one or more aligned sequences by first generating representative profiles, then aligning the profiles with dynamic programming, and finally merging the two aligned sequence sets.
}
\usage{
AlignProfiles(pattern,
              subject,
              p.weight = 1,
              s.weight = 1,
              p.struct = NULL,
              s.struct = NULL,
              perfectMatch = 5,
              misMatch = 0,
              gapOpening = -13,
              gapExtension = -1,
              gapPower = -1,
              terminalGap = -5,
              restrict = -1000,
              anchor = 0.7,
              normPower = 1,
              substitutionMatrix = NULL,
              structureMatrix = NULL,
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
  \item{p.struct}{
Either \code{NULL} (the default), or a list of matrices with one list element per sequence in the pattern.  (See details section below.)
}
  \item{s.struct}{
Either \code{NULL} (the default), or a list of matrices with one list element per sequence in the subject.  (See details section below.)
}
  \item{perfectMatch}{
Numeric giving the reward for aligning two matching nucleotides in the alignment.  Only used for \code{DNAStringSet} or \code{RNAStringSet} inputs.
}
  \item{misMatch}{
Numeric giving the cost for aligning two mismatched nucleotides in the alignment.  Only used for \code{DNAStringSet} or \code{RNAStringSet} inputs.
}
  \item{gapOpening}{
Numeric giving the cost for opening a gap in the alignment.
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap in the alignment.
}
  \item{gapPower}{
Numeric specifying the exponent to use in the gap cost function.  (See details section below.)
}
  \item{terminalGap}{
Numeric giving the cost for allowing leading and trailing gaps ("-" or "." characters) in the alignment.  Either two numbers, the first for leading gaps and the second for trailing gaps, or a single number for both.
}
  \item{restrict}{
Numeric specifying the lowest relative score to consider when aligning.  The default (\code{-1000}) will align most inputs that can reasonably be globally aligned without any loss in accuracy.  Input sequences with high similarity could be more restricted (e.g., \code{-500}), whereas a \code{pattern} and \code{subject} with little overlap should be less restricted (e.g., \code{-10000}).  (See details section below.)
}
  \item{anchor}{
Numeric giving the fraction of sequences with identical k-mers required to become an anchor point, or \code{NA} to not use anchors.  Alternatively, a matrix specifying anchor regions.  (See details section below.)
}
  \item{normPower}{
Numeric giving the exponent that controls the degree of normalization applied to scores by column occupancy.  A \code{normPower} of \code{0} does not normalize the scores, which results in all columns of the profiles being weighted equally, and is the optimal value for aligning fragmentary sequences.  A \code{normPower} of \code{1} normalizes the score for aligning two positions by their column occupancy (1 - fraction of gaps).  A \code{normPower} greater than \code{1} more strongly discourages aligning with ``gappy'' regions of the alignment.
}
  \item{substitutionMatrix}{
Either a substitution matrix representing the substitution scores for an alignment or the name of the amino acid substitution matrix to use in alignment.  The latter may be one of the following:  ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250'', or ``MIQS''.  The default (NULL) will use the \code{perfectMatch} and \code{misMatch} penalties for DNA/RNA or ``MIQS'' for AA.  (See examples section below.)
}
  \item{structureMatrix}{
A structure matrix if \code{p.struct} and \code{s.struct} are supplied, or \code{NULL} otherwise.  (See examples section below.)
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
}
\details{
Profiles are aligned using dynamic programming, a variation of the Needleman-Wunsch algorithm for global alignment.  The dynamic programming method requires order \code{N*M} time and memory space where \code{N} and \code{M} are the width of the pattern and subject.  This method works by filling in a matrix of the possible ``alignment space'' by considering all matches, insertions, and deletions between two sequence profiles.  The highest scoring alignment is then used to add gaps to each of the input sequence sets.

Heuristics can be useful to improve performance on long input sequences.  The \code{restrict} parameter can be used to dynamically constrain the possible ``alignment space'' to only paths that will likely include the final alignment, which in the best case can improve the speed from quadratic time to linear time.  The degree of restriction is important, and if the sequences are not mostly overlapping then \code{restrict} should be relaxed (more negative than the default).  For example, if aligning a pattern to a much longer subject then restrict should be set to \code{-1e10} to prevent restriction.

The argument \code{anchor} can be used to split the global alignment into multiple sub-alignments.  This can greatly decrease the memory requirement for long sequences when appropriate anchor points can be found.  Anchors are 15-mer (for DNA/RNA) or 7-mer (for AA) subsequences that are shared between at least \code{anchor} fraction of \code{pattern}(s) and \code{subject}(s).  Anchored ranges are extended along the length of each sequence in a manner designed to split the alignment into sub-alignments that can be separately solved.  For most input sequences \code{anchor} has no effect on accuracy, but anchoring can be disabled by setting \code{anchor=NA}.

Alternatively, \code{anchor} can be a matrix with \code{4} rows and one column per anchor.  The first two rows correspond to the anchor start and end positions in the pattern sequence(s), and the second two rows are the equivalent anchor region in the subject sequence(s).  Anchors specified in this manner must be strictly increasing (non-overlapping) in both sequences, and have an anchor width of at least two positions.  Note that the anchors do not have to be equal length, in which case the anchor regions will also be aligned.  Manually splitting the alignment into more subtasks can sometimes make it more efficient, but typically automatic anchoring is effective.

The arguments \code{p.struct} and \code{s.struct} may be used to provide secondary structure probabilities.  Each list element must contain a matrix with dimensions \code{q*w}, where \code{q} is the number of possible secondary structure states, and \code{w} is the width of the unaligned pattern sequence.  Values in each matrix represent the probability of the given state at that position in the sequence.  A \code{structureMatrix} must be supplied along with the structures.  The function \code{PredictHEC} can be used to predict secondary structure probabilities in the format required by \code{AlignProfiles}.

The gap cost function is based on the observation that gap lengths are best approximated by a Zipfian distribution (Chang & Benner, 2004).  The cost of inserting a gap of length \code{L} is equal to:
\code{gapOpening + gapExtension*sum(seq_len(L - 1)^gapPower)}
when \code{L > 1}, and \code{gapOpen} when \code{L = 1}.  This function effectively penalizes shorter gaps significantly more than longer gaps when \code{gapPower < 0}, and is equivalent to the affine gap penalty when \code{gapPower} is \code{0}.
}
\value{
An \code{XStringSet} of aligned sequences.
}
\references{
Chang, M. S. S., & Benner, S. A. (2004). Empirical Analysis of Protein Insertions and Deletions Determining Parameters for the Correct Placement of Gaps in Protein Sequence Alignments. Journal of Molecular Biology, \bold{341(2)}, 617-631.

Needleman S., Wunsch, C. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology, \bold{48(3)}, 443-453.

ES Wright (2015) "DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment". BMC Bioinformatics, doi:10.1186/s12859-015-0749-z.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AlignDB}}, \code{\link{AlignSeqs}}, \code{\link{AlignSynteny}}, \code{\link{AlignTranslation}}, \code{\link{MIQS}}
}
\examples{
# align two sets of sequences
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna1 <- SearchDB(db, remove="common", limit=100) # the first 100 sequences
dna2 <- SearchDB(db, remove="common", limit="100,100") # the rest
alignedDNA <- AlignProfiles(dna1, dna2)
BrowseSeqs(alignedDNA, highlight=1)

# specify a DNA substitution matrix
subMatrix <- matrix(0,
                    nrow=4, ncol=4,
                    dimnames=list(DNA_BASES, DNA_BASES))
diag(subMatrix) <- 5 # perfectMatch
alignedDNA.defaultSubM <- AlignProfiles(dna1, dna2, substitutionMatrix=subMatrix)
all(alignedDNA.defaultSubM==alignedDNA)

# specify a different DNA substitution matrix
subMatrix2 <- matrix(c(12, 3, 5, 3, 3, 12, 3, 6, 5, 3, 11, 3, 3, 6, 3, 9),
                    nrow=4, ncol=4,
                    dimnames=list(DNA_BASES, DNA_BASES))
alignedDNA.alterSubM <- AlignProfiles(dna1, dna2, substitutionMatrix=subMatrix2)
all(alignedDNA.alterSubM==alignedDNA)

# anchors are found automatically by default, but it is also
# possible to specify anchor regions between the sequences
anchors <- matrix(c(774, 788, 752, 766), nrow=4)
anchors
subseq(dna1, anchors[1, 1], anchors[2, 1])
subseq(dna2, anchors[3, 1], anchors[4, 1])
alignedDNA2 <- AlignProfiles(dna1, dna2, anchor=anchors)
}
