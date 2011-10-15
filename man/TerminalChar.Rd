\name{TerminalChar}
\alias{TerminalChar}
\title{
Determine the Number of Terminal Characters
}
\description{
Counts the number of terminal characters for every sequence in a \code{DNAStringSet}.  Terminal characters are defined as a specific character repeated at the beginning and end of a sequence.
}
\usage{
TerminalChar(myDNAStringSet,
             char = "-")
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object of sequences.
}
  \item{char}{
A single character giving the terminal character to count.
}
}
\value{
A \code{matrix} containing the results for each sequence in its respective row.  The first column contains the number of leading \code{char}, the second contains the number of trailing \code{char}, and the third contains the total number of characters in between.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{IdLengths}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db)
t <- TerminalChar(dna)
}