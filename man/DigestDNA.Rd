\name{DigestDNA}
\alias{DigestDNA}
\title{
Simulate Restriction Digestion of DNA
}
\description{
Restriction enzymes can be used to cut double-stranded DNA into fragments at specific cut sites.  \code{DigestDNA} performs an \emph{in-silico} restriction digest of the input DNA sequence(s) given one or more restriction sites.
}
\usage{
DigestDNA(sites,
          myDNAStringSet,
          type = "fragments",
          strand = "both",
          processors = NULL)
}
\arguments{
  \item{sites}{
A character vector of DNA recognition sequences and their enzymes' corresponding cut site(s).
}
  \item{myDNAStringSet}{
A \code{DNAStringSet} object or character vector with one or more sequences in 5' to 3' orientation.
}
  \item{type}{
Character string indicating the type of results desired.  This should be (an abbreviation of) either \code{"fragments"} or \code{"positions"}.
}
  \item{strand}{
Character string indicating the strand(s) to cut.  This should be (an unambiguous abbreviation of) one of \code{"both"}, \code{"top"}, or \code{"bottom"}.  The top strand is defined as the input DNAStringSet sequence, and the bottom strand is its reverse complement.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
}
\details{
In the context of a restriction digest experiment with a known DNA sequence, it can be useful to predict the expected DNA fragments \emph{in-silico}.  Restriction enzymes make cuts in double-stranded DNA at specific positions near their recognition site.  The recognition site may be somewhat ambiguous, as represented by the \code{IUPAC_CODE_MAP}.  Cuts that occur at different positions on the top and bottom strands result in sticky-ends, whereas those that occur at the same position result in fragments with blunt-ends.  Multiple restriction \code{sites} can be supplied to simultaneously digest the DNA.  In this case, \code{sites} for the different restriction enzymes may be overlapping, which could result in multiple close-proximity cuts that would not occur experimentally.  Also, note that cut sites will not be matched to non-\code{DNA_BASES} in \code{myDNAStringSet}.
}
\value{
\code{DigestDNA} can return two \code{type}s of results:  cut \code{positions} or the resulting DNA \code{fragments} corresponding to the \code{top}, \code{bottom}, or \code{both} strands.  If \code{type} is \code{"positions"} then the output is a list with the cut location(s) in each sequence in \code{myDNAStringSet}.  The cut location is defined as the position after the cut relative to the 5'-end.  For example, a cut at \code{6} would occur between positions 5 and 6, where the respective strand's 5' nucleotide is defined as position 1.

If \code{type} is \code{"fragments"} (the default), then the result is a \code{DNAStringSetList}.  Each element of the list contains the \code{top} and/or \code{bottom} strand fragments after digestion of \code{myDNAStringSet}, or the original sequence if no cuts were made.  Sequences are named by whether they originated from the \code{top} or \code{bottom} strand, and list elements are named based on the input DNA sequences.  The \code{top} strand is defined by \code{myDNAStringSet} is it is input, whereas the \code{bottom} strand is its reverse complement.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{DesignSignatures}}
}
\examples{
# digest hypothetical DNA sequences with BamHI
data(RESTRICTION_ENZYMES)
site <- RESTRICTION_ENZYMES[c("BamHI")]
dna <- DNAStringSet(c("AAGGATCCAA", "GGGATCAT"))
dna # top strand
reverseComplement(dna) # bottom strand
names(dna) <- c("hyp1", "hyp2")
d <- DigestDNA(site, dna)
d # fragments in a DNAStringSetList
unlist(d) # all fragments as one DNAStringSet

# Restriction digest of Yeast Chr. 1 with EcoRI and EcoRV
data(yeastSEQCHR1)
sites <- RESTRICTION_ENZYMES[c("EcoRI", "EcoRV")]
seqs <- DigestDNA(sites, yeastSEQCHR1)
seqs[[1]]

pos <- DigestDNA(sites, yeastSEQCHR1, type="positions")
str(pos)
}