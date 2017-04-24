\name{CorrectFrameshifts}
\alias{CorrectFrameshifts}
\title{
Corrects Frameshift Errors In Protein Coding Sequences
}
\description{
Corrects the reading frame to mitigate the impact of frameshift errors caused by insertions or deletions in unaligned nucleotide sequences.
}
\usage{
CorrectFrameshifts(myXStringSet,
                   myAAStringSet,
                   type = "indels",
                   acceptDistance = 0.01,
                   rejectDistance = 0.60,
                   maxComparisons = 10,
                   gapOpening = -13,
                   gapExtension = -1,
                   frameShift = -15,
                   geneticCode = GENETIC_CODE,
                   substitutionMatrix = "MIQS",
                   verbose = TRUE,
                   processors = 1)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} of unaligned protein coding sequences in 5' to 3' orientation.
}
  \item{myAAStringSet}{
An \code{AAStringSet} of reference sequences having accurate translations.
}
  \item{type}{
Character string indicating the type of result desired.  This should be (an abbreviation of) one of \code{"indels"}, \code{"sequences"}, or \code{"both"}.  (See details section below.)
}
  \item{acceptDistance}{
Numeric giving the maximum distance from a reference sequence that is acceptable to skip any remaining comparisons.
}
  \item{rejectDistance}{
Numeric giving the maximum distance from a reference sequence that is allowed when correcting frameshifts.  Sequences in \code{myXStringSet} that are greater than \code{rejectDistance} from the nearest reference sequence will only have their length trimmed from the 3'-end to a multiple of three nucleotides without any frameshift correction.
}
  \item{maxComparisons}{
The number of reference comparisons to make before stopping the search for a closer reference sequence.
}
  \item{gapOpening}{
Numeric giving the cost for opening a gap between the query and reference sequences.
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap between the query and reference sequences.
}
  \item{frameShift}{
Numeric giving the cost for shifting between frames of the query sequence.
}
  \item{geneticCode}{
Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code.
}
  \item{substitutionMatrix}{
Either a substitution matrix representing the substitution scores for matching two amino acids or the name of the amino acid substitution matrix.  The latter may be one of the following:  ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250'', or ``MIQS'' (the default).
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
}
\details{
Accurate translation of protein coding sequences can be greatly disrupted by one or two nucleotide phase shifts that occasionally occur during DNA sequencing.  These frameshift errors can potentially be corrected through comparison with other unshifted protein sequences.  This function uses a set of reference amino acid sequences (\code{AAStringSet}) to find and correct frameshift errors in a set of nucleotide sequences (\code{myXStringSet}).  First, three frame translation of the nucleotide sequences is performed, and the nearest reference sequence is selected.  Then the optimal reading frame at each position is determined based on a variation of the Guan & Uberbacher (1996) method.  Putative insertions and/or deletions (indels) are returned in the result, typically with close proximity to the true indel locations.  For a comparison of this method to others, see Wang et al. (2013).

If \code{type} is \code{"sequences"} or \code{"both"}, then frameshifts are corrected by adding \code{N}'s and/or removing nucleotides.  Note that this changes the nucleotide sequence, and the new sequence often has minor errors because the exact location of the indel(s) cannot be determined.  However, the original frameshifts that disrupted the entire downstream sequence are reduced to local perturbations.  All of the returned nucleotide sequences will have a reading frame starting from the first position.  This allows direct translation, and in practice works well if there is a similar reference \code{myAAStringSet} with the correct reading frame.  Hence it is more important that \code{myAAStringSet} contain a wide variety of sequences than it is that it contain a lot of sequences.

Multiple inputs control the time required for frameshift correction.  The number of sequences in the reference set (\code{myAAStringSet}) will affect the speed of the first search for similar sequences.  Assessing frameshifts in the second step requires order \code{N*M} time, where \code{N} and \code{M} are the lengths of the query (\code{myXStringSet}) and reference sequences.  Two parameters control the number of assessments that are made for each sequence:  (1) \code{maxComparisons} determines the maximum number of reference sequences to compare to each query sequence, and (2) \code{acceptDist} defines the maximum distance between a query and reference that is acceptable before continuing to the next query sequence.  A lower value for \code{maxComparisons} or a higher value for \code{acceptDist} will accelerate frameshift correction, potentially at the expense of some accuracy.
}
\value{
If \code{type} is \code{"indels"} then the returned object is a list with the same length as \code{myXStringSet}.  Each element is a list with four components:
  \item{insertions}{
Approximate positions of inserted nucleotides, which could be removed to correct the reading frame, or excess nucleotides at the 3'-end that make the length longer than a multiple of three.
}
  \item{deletions}{
Approximate positions of deleted nucleotides, which could be added back to correct the reading frame.
}
  \item{distance}{
The amino acid distance from the nearest reference sequence, between 0 and 1.
}
  \item{index}{
The integer index of the reference sequence that was used for frame correction, or \code{0} if no reference sequence was within \code{rejectDistance}.
}
Note that positions in \code{insertions} and \code{deletions} are sometimes repeated to indicate that the same position needs to be shifted successively more than once to correct the reading frame.

If \code{type} is \code{"sequences"} then the returned object is an \code{XStringSet} of the same type as the input (\code{myXStringSet}).  Nucleotides are added or deleted as necessary to correct for frameshifts.  The returned sequences all have a reading frame starting from position 1, so that they can be translated directly.

If \code{type} is \code{"both"} then the returned object is a list with two components:  one for the \code{"indels"} and the other for the \code{"sequences"}.
}
\references{
Guan, X., & Uberbacher, E. C. (1996). Alignments of DNA and protein sequences containing frameshift errors. Computer Applications in the Biosciences : CABIOS, \bold{12(1)}, 31-40.

Wang, Q., et al. (2013). Ecological Patterns of nifH Genes in Four Terrestrial Climatic Zones Explored with Targeted Metagenomics Using FrameBot, a New Informatics Tool. mBio, \bold{4(5)}, e00592-13-e00592-13.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignTranslation}}, \code{\link{OrientNucleotides}}
}
\examples{
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

# introduce artificial indels
n_ins <- 2 # insertions per sequence
shifted <- replaceAt(dna,
	lapply(width(dna),
		sample,
		n_ins),
	sample(DNA_BASES,
		n_ins,
		replace=TRUE))
n_dels <- 1 # deletions per sequence
shifted <- replaceAt(shifted,
	RangesList(lapply(width(shifted),
		function(x) {
			IRanges(sample(x,
					n_dels),
				width=1)
		})))

# to make frameshift correction more challenging,
# only supply 20 reference amino acid sequences
s <- sample(length(dna), 20)
x <- CorrectFrameshifts(shifted,
	translate(dna[s]),
	type="both")

# there was a wide range of distances
# to the nearest reference sequence
quantile(unlist(lapply(x[[1]], `[`, "distance")))

# none of the sequences were > rejectDistance
# from the nearest reference sequence
length(which(unlist(lapply(x[[1]], `[`, "index"))==0))

# the number of indels was generally correct
table(unlist(lapply(x[[1]], function(x) {
	length(x$insertions)})))/length(shifted)
table(unlist(lapply(x[[1]], function(x) {
	length(x$deletions)})))/length(shifted)

# align and display the translations
AA <- AlignTranslation(x$sequences,
	readingFrame=1,
	type="AAStringSet")
BrowseSeqs(AA)
}
