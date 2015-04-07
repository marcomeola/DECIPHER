\name{Disambiguate}
\alias{Disambiguate}
\title{
Expand Ambiguities into All Permutations of a DNAStringSet
}
\description{
Performs the inverse function of \code{ConsensusSequence} by expanding any ambiguities present in sequences.
}
\usage{
Disambiguate(myXStringSet)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} object of sequences.
}
}
\details{
Ambiguity codes in the \code{IUPAC_CODE_MAP} can be used to represent multiple nucleotides at a single position.  Using these letters, multiple oligonucleotide permutations can be represented with a single ambiguous sequence.  This function expands each sequence in the \code{DNAStringSet} input into all of its permutations.  Note that sequences with many ambiguities can result in a very large number of potential permutations.
}
\value{
A \code{DNAStringSetList} or \code{RNAStringSetList} with one element for each sequence in \code{myXStringSet}.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{ConsensusSequence}}
}
\examples{
dna <- DNAStringSet(c("ACST", "NNN"))
dna_list <- Disambiguate(dna)
dna_list[[1]]
dna_list[[2]]
unlist(dna_list)

rna <- RNAStringSet(c("ACGU", "AGAU")) # 2 permutations
rna <- ConsensusSequence(rna) # "ASRU"
Disambiguate(rna) # 4 permutations
}