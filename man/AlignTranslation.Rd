\name{AlignTranslation}
\alias{AlignTranslation}
\title{
Aligns Sequences By Their Amino Acid Translation
}
\description{
Performs alignment of a set of DNA or RNA sequences by aligning their corresponding amino acid sequences.
}
\usage{
AlignTranslation(myXStringSet,
                 sense = "+",
                 direction = "5' to 3'",
                 readingFrame = NA,
                 asAAStringSet = FALSE,
                 ...)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} object of unaligned sequences.
}
  \item{sense}{
Single character specifying sense of the input sequences, either the positive (\code{"+"}) coding strand or negative (\code{"-"}) non-coding strand.
}
  \item{direction}{
Direction of the input sequences, either \code{"5' to 3'"} or \code{"3' to 5'"}.
}
  \item{readingFrame}{
Numeric vector giving a single reading frame for all of the sequences, or an individual reading frame for each sequence in \code{myXStringSet}.  The \code{readingFrame} can be either \code{1}, \code{2}, \code{3} to begin translating on the first, second, and third nucleotide position, or \code{NA} (the default) to guess the reading frame.  (See details section below.)
}
  \item{asAAStringSet}{
Logical specifying whether to return the aligned translation as an \code{AAStringSet} rather than an \code{XStringSet} of the input type.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{AlignSeqs}}, including \code{perfectMatch}, \code{misMatch}, \code{gapOpening}, \code{gapExtension}, \code{terminalGap}, \code{restrict}, \code{anchor}, \code{doNotAlign}, \code{guideTree}, \code{processors}, and \code{verbose}.
}
}
\details{
Alignment of proteins is often more accurate than alignment of their coding nucleic acid sequences.  This function aligns the input nucleic acid sequences via aligning their translated amino acid sequence.  First, the input sequences are translated according to the specified \code{sense}, \code{direction}, and \code{readingFrame}.  The resulting amino acid sequences are aligned with the function \code{AlignSeqs}, and this alignment is reverse translated into the original sequence type, \code{sense}, and \code{direction}.

If the \code{readingFrame} is \code{NA} (the default) then an attempt is made to guess the reading frame of each sequence based on the number of stop codons in the translated amino acid sequence.  The first reading frame will be chosen (\code{1}, \code{2}, or \code{3}) that has one or fewer stop codons.  If the number of stop codons is inconclusive for a sequence then reading frame will default to \code{1}.  The entire length of each sequence is translated in spite of any stop codons identified.  Note that this method is only constructive in circumstances where there is a substantially long coding sequence with at most a single stop codon expected, and therefore it is preferable to specify the reading frame of each sequence if it is known.
}
\value{
An \code{XStringSet} matching the input type.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AlignProfiles}}, \code{\link{AlignSeqs}}
}
\examples{
# First three sequences translate to MFITP*, last sequence is MF-TP*
rna <- RNAStringSet(c("AUGUUCAUCACCCCCUAA", "AUGUUCAUAACUCCUUGA",
	"AUGUUCAUUACACCGUAG", "AUGUUUACCCCAUAA"))
RNA <- AlignSeqs(rna)
BrowseSequences(RNA) # incorrect gap position

RNA <- AlignTranslation(rna)
BrowseSequences(RNA) # correct gap position
}
