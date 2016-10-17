\name{DistanceMatrix}
\alias{DistanceMatrix}
\title{
Calculate the Distances Between Sequences
}
\description{
Calculates a distance matrix for an \code{XStringSet}.  Each element of the distance matrix corresponds to the dissimilarity between two sequences in the \code{XStringSet}.
}
\usage{
DistanceMatrix(myXStringSet,
               includeTerminalGaps = FALSE,
               penalizeGapLetterMatches = TRUE,
               penalizeGapGapMatches = FALSE,
               correction = "none",
               processors = 1,
               verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
An \code{XStringSet} object of aligned sequences (\code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet}).
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" characters on each end of the sequence) into the calculation of distance.
}
  \item{penalizeGapLetterMatches}{
Logical specifying whether or not to consider gap-to-letter matches as mismatches.  If \code{FALSE}, then gap-to-letter matches are not included in the total length used to calculate distance.
}
  \item{penalizeGapGapMatches}{
Logical specifying whether or not to consider gap-to-gap matches as mismatches.  If \code{FALSE} (the default), then gap-to-gap matches are not included in the total length used to calculate distance.
}
  \item{correction}{
The substitution model used for distance correction.  This should be (an abbreviation of) either \code{"none"} or \code{"Jukes-Cantor"}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The uncorrected distance matrix represents the hamming distance between each of the sequences in \code{myXStringSet}.  Ambiguity can be represented using the characters of the \code{IUPAC_CODE_MAP} for \code{DNAStringSet} and \code{RNAStringSet} inputs, or using the \code{AMINO_ACID_CODE} for an \code{AAStringSet} input.  For example, the distance between an 'N' and any other nucleotide base is zero.  The letters B (N or D), J (I or L), Z (Q or E), and X (any letter) are degenerate in the \code{AMINO_ACID_CODE}.

If \code{includeTerminalGaps = FALSE} then terminal gaps ("-" or "." characters) are not included in sequence length.  This can be faster since only the positions common to each pair of sequences are compared.    Sequences with no overlapping region in the alignment are given a value of \code{NA}, unless \code{includeTerminalGaps = TRUE}, in which case distance is 100\%.

Penalizing gap-to-gap and gap-to-letter mismatches specifies whether to penalize these special mismatch types and include them in the total length when calculating distance.  Both "-" and "." characters are interpreted as gaps.  The default behavior is to calculate distance as the fraction of positions that differ across the region of the alignment shared by both sequences (not including gap-to-gap matches).

The elements of the distance matrix can be referenced by \code{dimnames} corresponding to the \code{names} of the \code{XStringSet}.  Additionally, an attribute named "correction" specifying the method of correction used can be accessed using the function \code{attr}.
}
\value{
A symmetric matrix where each element is the distance between the sequences referenced by the respective row and column.  The \code{dimnames} of the matrix correspond to the \code{names} of the \code{XStringSet}.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{IdClusters}}
}
\examples{
# defaults compare intersection of internal ranges:
dna <- DNAStringSet(c("ANGCT-","-ACCT-"))
d <- DistanceMatrix(dna)
# d[1,2] is 1 base in 4 = 0.25

# compare the entire sequence ranges:
dna <- DNAStringSet(c("ANGCT-","-ACCT-"))
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE,
                    penalizeGapGapMatches=TRUE)
# d[1,2] is now 3 bases in 6 = 0.50

# compare union of internal ranges:
dna <- DNAStringSet(c("ANGCT-","-ACCT-"))
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE,
                    penalizeGapGapMatches=FALSE)
# d[1,2] is now 2 bases in 5 = 0.40

# gap ("-") and unknown (".") characters are interchangeable:
dna <- DNAStringSet(c("ANGCT.",".ACCT-"))
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE,
                    penalizeGapGapMatches=FALSE)
# d[1,2] is still 2 bases in 5 = 0.40
}