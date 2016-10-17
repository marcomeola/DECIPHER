\name{ConsensusSequence}
\alias{ConsensusSequence}
\title{
Create a Consensus Sequence
}
\description{
Forms a consensus sequence representing a set of sequences.
}
\usage{
ConsensusSequence(myXStringSet,
                  threshold = 0.05,
                  ambiguity = TRUE,
                  noConsensusChar = "+",
                  minInformation = 1 - threshold,
                  ignoreNonBases = FALSE,
                  includeTerminalGaps = FALSE)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of aligned sequences.
}
  \item{threshold}{
Numeric specifying that less than \code{threshold} fraction of sequence information can be lost at any position of the consensus sequence.
}
  \item{ambiguity}{
Logical specifying whether to consider ambiguity as split between their respective nucleotides.  Degeneracy codes are specified in the \code{IUPAC_CODE_MAP}.
}
  \item{noConsensusChar}{
Single character from the sequence's alphabet giving the base to use when there is no consensus in a position.
}
  \item{minInformation}{
Minimum fraction of information required to form consensus in each position.
}
  \item{ignoreNonBases}{
Logical specifying whether to count gap ("-"), mask ("+"), and unknown (".") characters towards the consensus.
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" or "." characters on each end of the sequence) into the formation of consensus.
}
}
\details{
\code{ConsensusSequence} removes the least frequent characters at each position, so long as they represent less than \code{threshold} fraction of the sequences in total.  If necessary, \code{ConsensusSequence} represents the remaining characters using a degeneracy code from the \code{IUPAC_CODE_MAP}.  Degeneracy codes are always used in cases where multiple characters are equally abundant.

Two key parameters control the degree of consensus:  \code{threshold} and \code{minInformation}.  The default \code{threshold} (\code{0.05}) means that at less than 5\% of sequences will not be represented by the consensus sequence at any given position.  The default \code{minInformation} (\code{1 - 0.05}) specifies that at least 95\% of sequences must contain the information in the consensus, otherwise the \code{noConsensusChar} is used.  This enables an alternative character (e.g., "+") to be substituted at positions that would otherwise yield an ambiguity code.

If \code{ambiguity = TRUE} (the default) then degeneracy codes in \code{myXStringSet} are split between their respective bases according to the \code{IUPAC_CODE_MAP} for DNA/RNA and \code{AMINO_ACID_CODE} for AA.  For example, an ``R'' in a \code{DNAStringSet} would count as half an ``A'' and half a ``G''.  If \code{ambiguity = FALSE} then degeneracy codes are not considered in forming the consensus.  For an \code{AAStringSet} input, the lack of degeneracy codes generally results in ``X'' at positions with mismatches, unless the \code{threshold} is set to a higher value than the default.

If \code{includeNonBases = TRUE} (the default) then gap ("-"), mask ("+"), and unknown (".") characters are counted towards the consensus, otherwise they are omitted from calculation of the consensus.  Note that gap ("-") and unknown (".") characters are treated interchangeably as gaps when forming the consensus sequence.  For this reason, the consensus of a position with all unknown (".") characters will be a gap ("-").  Also, note that if consensus is formed between different length sequences then it will represent only the longest sequences at the end.  For this reason the consensus sequence is generally based on a sequence alignment such that all of the sequences have equal lengths.
}
\value{
An \code{XStringSet} with a single consensus sequence matching the input type.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{Disambiguate}}, \code{\link{IdConsensus}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db, limit=10)
BrowseSeqs(dna) # consensus at bottom
BrowseSeqs(dna, threshold=0.5) # consensus at bottom

# controlling the degree of consensus
AAAT <- DNAStringSet(c("A", "A", "A", "T"))
ConsensusSequence(AAAT) # "W"
ConsensusSequence(AAAT, threshold=0.3) # "A"
ConsensusSequence(AAAT, threshold=0.3, minInformation=0.8) # "+"
ConsensusSequence(AAAT, threshold=0.3, minInformation=0.8, noConsensusChar="N") # "N"

# switch between degenerate-based and majority-based consensus
majority <- DNAStringSet(c("GTT", "GAA", "CTG"))
ConsensusSequence(majority) # degenerate-based
ConsensusSequence(majority, threshold=0.5) # majority-based
ConsensusSequence(majority, threshold=0.5, minInformation=0.75)

# behavior in the case of a tie
ConsensusSequence(DNAStringSet(c("A", "T"))) # "W"
ConsensusSequence(DNAStringSet(c("A", "T")), threshold=0.5) # "W"
ConsensusSequence(AAStringSet(c("A", "T"))) # "X"
ConsensusSequence(AAStringSet(c("A", "T")), threshold=0.5) # "X"
ConsensusSequence(AAStringSet(c("I", "L"))) # "J"
ConsensusSequence(AAStringSet(c("I", "L")), threshold=0.5) # "J"

# handling terminal gaps
dna <- DNAStringSet(c("ANGCT-","-ACCT-"))
ConsensusSequence(dna) # "ANSCT-"
ConsensusSequence(dna, includeTerminalGaps=TRUE) # "+NSCT-"

# the "." character is treated is a "-"
aa <- AAStringSet(c("ANQIH-", "ADELW."))
ConsensusSequence(aa) # "ABZJX-"

# internal non-bases are included by default
ConsensusSequence(DNAStringSet(c("A-+.A", "AAAAA")), noConsensusChar="N") # "ANNNA"
ConsensusSequence(DNAStringSet(c("A-+.A", "AAAAA")), ignoreNonBases=TRUE) # "AAAAA"

# degeneracy codes in the input are considered by default
ConsensusSequence(DNAStringSet(c("AWNDA", "AAAAA"))) # "AWNDA"
ConsensusSequence(DNAStringSet(c("AWNDA", "AAAAA")), ambiguity=FALSE) # "AAAAA"
}