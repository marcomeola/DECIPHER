\name{IdLengths}
\alias{IdLengths}
\title{
Determine the Number of Bases, Nonbases, and Width of Each Sequence
}
\description{
Counts the number of bases (A, C, G, T) and ambiguities/degeneracies in each sequence.
}
\usage{
IdLengths(dbFile,
          tblName = "DNA",
          identifier = "",
          type = "DNAStringSet",
          add2tbl = FALSE,
          batchSize = 10000,
          verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table where the sequences are located.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.
}
  \item{type}{
The type of \code{XStringSet} being processed.  This should be (an abbreviation of) one of \code{"DNAStringSet"} or \code{"RNAStringSet"}.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{batchSize}{
Integer specifying the number of sequences to process at a time.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\value{
A \code{data.frame} with the number of \code{bases} (``A'', ``C'', ``G'', or ``T''), \code{nonbases}, and \code{width} of each sequence.  The \code{width} is defined as the sum of \code{bases} and \code{nonbases} in each sequence.  The \code{row.names} of the \code{data.frame} correspond to the "row_names" in the \code{tblName} of the \code{dbFile}.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{Add2DB}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
l <- IdLengths(db)
}