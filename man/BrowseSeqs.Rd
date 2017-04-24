\name{BrowseSeqs}
\alias{BrowseSeqs}
\title{
View Sequences in a Web Browser
}
\description{
Opens an html file in a web browser to show the sequences in an \code{XStringSet}.
}
\usage{
BrowseSeqs(myXStringSet,
           htmlFile = paste(tempdir(), "/myXStringSet.html", sep = ""),
           openURL = interactive(),
           colorPatterns = TRUE,
           highlight = NA,
           patterns = c("-", alphabet(myXStringSet, baseOnly=TRUE)),
           colors = substring(rainbow(length(patterns),
                              v=0.8, start=0.9, end=0.7), 1, 7),
           colWidth = Inf,
           \dots)
}
\arguments{
  \item{myXStringSet}{
A \code{XStringSet} object of sequences.
}
  \item{htmlFile}{
Character string giving the location where the html file should be written.
}
  \item{openURL}{
Logical indicating whether the \code{htmlFile} should be opened in a web browser.
}
  \item{colorPatterns}{
Logical specifying whether to color matched \code{patterns}, or an integer vector providing pairs of start and stop boundaries for coloring.
}
  \item{highlight}{
Numeric specifying which sequence in the set to use for comparison or \code{NA} to color all sequences (default).  If \code{highlight} is \code{0} then positions differing from the consensus sequence are highlighted.
}
  \item{patterns}{
Either an \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object, a character vector containing regular expressions, a list of numeric matrices, or \code{NULL}.  (See details section below.)
}
  \item{colors}{
Character vector providing the color for each of the matched \code{patterns}.  Typically a character vector with elements of 7 characters:  ``#'' followed by the red, blue, green values in hexadecimal (after rescaling to 0 ... 255).  Ignored when \code{patterns} is a list of matrices.
}
  \item{colWidth}{
Integer giving the maximum number of nucleotides wide the display can be before starting a new page.  Must be a multiple of \code{20} (e.g., \code{100}), or \code{Inf} (the default) to display all the sequences in one set of rows.
}
  \item{\dots}{
Additional arguments to adjust the appearance of the consensus sequence at the base of the display.  Passed directly to \code{ConsensusSequence} for an \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet}, or to \code{consensusString} for a \code{BStringSet}.
}
}
\details{
\code{BrowseSeqs} converts an \code{XStringSet} into html format for viewing in a web browser.  The sequences are colored in accordance with the \code{patterns} that are provided, or left uncolored if \code{colorPatterns} is \code{FALSE} or \code{patterns} is \code{NULL}.  Character or \code{XStringSet} patterns are matched as regular expressions and colored according to \code{colors}.  If \code{patterns} is a list of matrices, then it must contain one element per sequence.  Each matrix is interpreted as providing the fraction red, blue, and green for each letter in the sequence.  Thus, \code{colors} is ignored when \code{patterns} is a list.  (See examples section below.)

Patterns are not matched across column breaks, so multi-character \code{patterns} should be carefully considered when \code{colWidth} is less than the maximum sequence length.  Patterns are matched sequentially in the order provided, so it is feasible to use nested \code{patterns} such as \code{c("ACCTG", "CC")}.  In this case the ``CC'' could be colored differently inside the previously colored ``ACCTG''.  Note that \code{patterns} overlapping the boundaries of a previously matched pattern will not be matched.  For example, ``ACCTG'' would not be matched if \code{patterns=c("CC", "ACCTG")}.

Some web browsers cannot quickly display a large amount colored text, so it is recommended to use \code{colorPatterns = FALSE} or to \code{highlight} a sequence when viewing a large \code{XStringSet}.  Highlighting will only show all of the characters in the highlighted sequence, and convert all matching positions in the other sequences into dots without \code{color}.  Also, note that some web browsers display small shifts between fixed-width characters that may become noticeable as color offsets between the ends of long sequences.
}
\value{
Creates an html file containing sequence data and (if \code{openURL} is \code{TRUE}) opens it in a web browser for viewing.  The layout has the sequence name on the left, position legend on the top, cumulative number of nucleotides on the right, and consensus sequence on the bottom.

Returns \code{htmlFile} if the html file was written successfully.
}
\references{
ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{BrowseDB}}, \code{\link{ConsensusSequence}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db) # non-coding RNA sequences
BrowseSeqs(dna)
BrowseSeqs(dna, colWidth=100, highlight=1)
dbn <- PredictDBN(dna, type="structures")
BrowseSeqs(dna, patterns=dbn) # color according to base-pairing

# color amino acid sequences
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
aa <- unique(translate(dna)) # the unique amino acid sequences
AA <- AlignSeqs(aa)
hec <- PredictHEC(AA, type="probabilities")
BrowseSeqs(AA, patterns=hec) # color according to secondary structure

# color bases in alternating groups with a different color scheme
BrowseSeqs(dna[1:5],
	colorPatterns=seq(1, width(dna)[1], 10),
	patterns=c("A", "C", "G", "T", "-"),
	colors=c("#1E90FF", "#32CD32", "#9400D3", "#000000", "#EE3300"))

# color all restriction sites
data(RESTRICTION_ENZYMES)
sites <- RESTRICTION_ENZYMES
sites <- gsub("[^A-Z]", "", sites) # remove non-letters
sites <- DNAStringSet(sites)
rc_sites <- DNAStringSet(sites)
w <- which(sites != rc_sites)
sites <- c(sites, rc_sites[w])
sites <- sites[order(nchar(sites))] # match shorter sites first

dna <- SearchDB(db, remove="all") # unaligned sequences
BrowseSeqs(dna, patterns=sites)
}