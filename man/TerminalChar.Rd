\name{TerminalChar}
\alias{TerminalChar}
\title{
Determine the Number of Terminal Characters
}
\description{
Counts the number of terminal characters for every sequence in an \code{XStringSet}.  Terminal characters are defined as a specific character repeated at the beginning and end of a sequence.
}
\usage{
TerminalChar(myXStringSet,
             char = "")
}
\arguments{
  \item{myXStringSet}{
An \code{XStringSet} object of sequences.
}
  \item{char}{
A single character giving the terminal character to count, or an empty character ("") indicating to count both gap ("-") and unknown (".") characters.
}
}
\value{
A \code{matrix} containing the results for each sequence in its respective row.  The first column contains the number of leading \code{char}, the second contains the number of trailing \code{char}, and the third contains the total number of characters in-between.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{IdLengths}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db)
t <- TerminalChar(dna)
}