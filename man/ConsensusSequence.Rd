\name{ConsensusSequence}
\alias{ConsensusSequence}
\title{
Create A Consensus Sequence
}
\description{
Forms a consensus sequence representing a set of sequences.
}
\usage{
ConsensusSequence(myXStringSet,
                  threshold = 0.05,
                  ambiguity = TRUE,
                  noConsensusChar = "+",
                  minInformation = 0.75,
                  ignoreNonBases = FALSE,
                  includeTerminalGaps = FALSE,
                  verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of aligned sequences.
}
  \item{threshold}{
Maximum fraction of sequence information that may be lost in forming the consensus.
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
Logical specifying whether to count gap ("-") or mask ("+") characters towards the consensus.
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" characters on each end of the sequence) into the formation of consensus.
}
  \item{verbose}{
Logical indicating whether to print the elapsed time upon completion.
}
}
\details{
Two key parameters control the degree of consensus.  The default \code{threshold} (0.05) indicates that at least 95\% of sequence information will be represented by the consensus sequence.  The default \code{minInformation} (0.75) specifies that at least 75\% of sequences must contain the information in the consensus, otherwise the \code{noConsensusChar} is used.

If \code{ambiguity = TRUE} (the default) then degeneracy codes are split between their respective bases according to the \code{IUPAC_CODE_MAP}, or \code{AMINO_ACID_CODE} for \code{AAStringSet}s.  For example, an "R" in a \code{DNAStringSet} would count as half an "A" and half a "G".  If \code{ambiguity = FALSE} then degeneracy codes are not considered in forming the consensus.  If \code{includeNonBases = TRUE} (the default) then gap ("-") and mask ("+") characters are counted towards the consensus, otherwise they are omitted from calculation of the consensus.  For an \code{AAStringSet} input, the lack of degeneracy codes generally results in ``X'' in positions with mismatches, unless the \code{threshold} is set higher than 0.05.
}
\value{
An \code{XStringSet} matching the input type with a single consensus sequence.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{IdConsensus}}, \code{\link{Seqs2DB}}
}
\examples{
dna <- DNAStringSet(c("ANGCT-","-ACCT-"))
ConsensusSequence(dna)
# returns "ANSCT-"
}