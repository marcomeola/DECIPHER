\name{BrowseSequences}
\alias{BrowseSequences}
\title{
View Sequences In A Web Browser
}
\description{
Opens an html file in a web browser to show the sequences in an \code{XStringSet}.
}
\usage{
BrowseSequences(myXStringSet,
                htmlFile = paste(tempdir(), "/myXStringSet.html", sep = ""),
                colorPatterns = TRUE,
                highlight = 0,
                patterns = c("-", alphabet(myXStringSet, baseOnly=TRUE)),
                colors = substring(rainbow(length(patterns),
                                 v=0.8, start=0.9, end=0.7), 1, 7),
                colWidth = Inf,
                ...)
}
\arguments{
  \item{myXStringSet}{
A \code{XStringSet} object of sequences.
}
  \item{htmlFile}{
Character string giving the location where the html file should be written.
}
  \item{colorPatterns}{
Logical specifying whether to color matched \code{patterns}, or an integer vector providing pairs of start and stop boundaries for coloring.
}
  \item{highlight}{
Numeric specifying which sequence in the set to use for comparison or \code{0} to color all sequences (default).
}
  \item{patterns}{
Character vector containing regular expressions to be colored in the XStringSet.  Regular expressions are searched sequentially with multiple matches allowed, even within other previously matched \code{patterns}.  (See details section below.)
}
  \item{colors}{
Character vector providing the color for each of the matched \code{patterns}.  Typically a character vector with elements of 7 characters, ``#'' followed by the red, blue, green values in hexadecimal (after rescaling to 0 ... 255).  Positions with background color have white font.
}
  \item{colWidth}{
Integer giving the maximum number of nucleotides wide the display can be before starting a new page.  Must be a multiple of \code{20} (e.g., \code{100}), or \code{Inf} (the default) to display all the sequences in one set of rows.
}
  \item{\dots}{
Additional arguments to adjust the appearance of the consensus sequence at the base of the display.  Passed directly to \code{ConsensusSequence} for an \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet}, or to \code{consensusString} for a \code{BStringSet}.
}
}
\details{
Some web browsers cannot quickly display a large amount colored text, so it is recommended to use \code{color = FALSE} when viewing a large \code{XStringSet}.

Patterns are not matched across column breaks, so multi-character \code{patterns} should be carefully considered when \code{colWidth} is less than the maximum sequence length.  Patterns are matched sequentially in the order provided, so it is feasible to use nested \code{patterns} such as \code{c("ACCTG", "CC")}.  In this case the ``CC'' could be colored differently inside the previously colored ``ACCTG''.  Note that \code{patterns} overlapping the boundaries of a previously matched pattern will not be matched.  For example, ``ACCTG'' would not be matched if \code{patterns=c("CC", "ACCTG")}.

Column positions identical to the \code{highlight}ed sequence are replaced with a ``.'' character.  Patterns are not matched in ``.'' regions.
}
\value{
Creates an html file containing sequence data and opens it in a web browser for easy viewing.  The viewer has the sequence name on the left, position legend on the top, cumulative number of nucleotides on the right, and consensus sequence on the bottom.

Returns \code{htmlFile} if the html file was written successfully.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{BrowseDB}}, \code{\link{ConsensusSequence}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db)
BrowseSequences(dna, colWidth=100, highlight=1)

# color bases in alternating groups with a different color scheme
BrowseSequences(dna[1:5],
	colorPatterns=seq(1, width(dna)[1], 10),
	patterns=c("A", "C", "G", "T", "-"),
	colors=c("#1E90FF", "#32CD32", "#9400D3", "#000000", "#EE3300"))

# color all restriction sites
data(RESTRICTION_ENZYMES)
sites <- RESTRICTION_ENZYMES
sites <- sites[order(nchar(sites))] # match shorter sites first
# convert all restriction sites into regular expressions
for (degeneracy in names(IUPAC_CODE_MAP)[5:15]) {
    sites <- gsub(degeneracy,
        paste("[", degeneracy, "|",
            paste(strsplit(IUPAC_CODE_MAP[degeneracy], "")[[1]],
                collapse="|"),
            "]",
            sep=""),
        sites)
}
dna <- SearchDB(db, remove="all") # unaligned sequences
BrowseSequences(dna, patterns=sites)
}