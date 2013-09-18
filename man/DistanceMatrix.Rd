\name{DistanceMatrix}
\alias{DistanceMatrix}
\title{
Calculate the Distance Between DNA Sequences
}
\description{
Calculates a distance matrix for a \code{DNAStringSet}.  Each element of the distance matrix corresponds to the dissimilarity between two sequences in the \code{DNAStringSet}.
}
\usage{
DistanceMatrix(myDNAStringSet,
               includeTerminalGaps = FALSE,
               penalizeGapLetterMatches = TRUE,
               penalizeGapGapMatches = FALSE,
               removeDuplicates = FALSE,
               correction = "none",
               processors = NULL,
               verbose = TRUE)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object of aligned sequences.
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" characters on each end of the sequence) into the calculation of distance.
}
  \item{penalizeGapLetterMatches}{
Logical specifying whether or not to consider gap-to-letter matches as mismatches.
}
  \item{penalizeGapGapMatches}{
Logical specifying whether or not to consider gap-to-gap matches as mismatches.
}
  \item{removeDuplicates}{
Logical specifying whether to remove any identical sequences from the input sequences before calculating distance.  If \code{FALSE} (the default) then the distance matrix is calculated with the entire \code{DNAStringSet} provided as input.
}
  \item{correction}{
The substitution model used for distance correction.  This should be (an unambiguous abbreviation of) one of \code{"none"} or \code{"Jukes-Cantor"}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The uncorrected distance matrix represents the percent distance between each of the sequences in the \code{DNAStringSet}.  Ambiguity can be represented using the characters of the \code{IUPAC_CODE_MAP}.  For example, the distance between an 'N' and any other base is zero.

If \code{includeTerminalGaps = FALSE} then terminal gaps are not included in sequence length.  This can be faster since only the positions common to each pair of sequences are compared.  Similarly, if \code{removeDuplicates = TRUE} then the distance matrix will only represent unique sequences in the \code{DNAStringSet}.  This is can be faster because less sequences need to be compared.  For example, if only two sequences in the set are exact duplicates then one is removed and the distance is calculated on the remaining set.  Note that the distance matrix can still contain values of 100\% after removing duplicates because only exact duplicates are removed without taking into account ambiguous matches represented by the \code{IUPAC_CODE_MAP} or the treatment of gaps.

The elements of the distance matrix can be referenced by \code{dimnames} corresponding to the \code{names} of the \code{DNAStringSet}.  Additionally, an attribute named "correction" specifying the method of correction used can be accessed using the function \code{attr}.
}
\value{
A symmetric matrix where each element is the distance between the sequences referenced by the respective row and column.  The \code{dimnames} of the matrix correspond to the \code{names} of the \code{DNAStringSet}.  Sequences with no overlapping positions in the alignment are given a value of \code{NA}.
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
# d[1,2] is still 1 base in 4 = 0.25

# compare union of internal ranges:
dna <- DNAStringSet(c("ANGCT-","-ACCT-"))
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE)
# d[1,2] is now 2 bases in 5 = 0.40

# compare the entire sequence ranges:
dna <- DNAStringSet(c("ANGCT-","-ACCT-"))
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE,
                    penalizeGapGapMatches=TRUE)
# d[1,2] is now 3 bases in 6 = 0.50
}