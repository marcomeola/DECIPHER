\name{BrowseSequences}
\alias{BrowseSequences}
\title{
View Sequences In A Web Browser
}
\description{
Opens an html file in a web browser to show the sequences in a \code{DNAStringSet}.
}
\usage{
BrowseSequences(myDNAStringSet,
                htmlFile = paste(tempdir(), "/dna.html", sep = ""),
                colorBases=TRUE,
                highlight=0,
                ...)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object of sequences.
}
  \item{htmlFile}{
Character string giving the location where the html file should be written.
}
  \item{colorBases}{
Logical specifying whether to color each type of base (A, C, G, and T) the same color.
}
  \item{highlight}{
Numeric specifying which sequence in the set to use for comparison or \code{0} to color every sequence (default).
}
  \item{\dots}{
Additional arguments to be passed directly to \code{ConsensusSequence}.
}
}
\details{
Some web browsers cannot quickly display a large amount data, so it is recommended to use \code{color = FALSE} (the default) when viewing a large \code{DNAStringSet}.
}
\value{
Creates an html file containing sequence data and opens it in a web browser for easy viewing.  The viewer has the sequence name on the left, position legend on the top, number of characters on the right, and consensus sequence on the bottom.

Returns \code{TRUE} if the html file was written successfully.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{BrowseDB}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db)
BrowseSequences(dna[1:5], colorBases=TRUE)
}