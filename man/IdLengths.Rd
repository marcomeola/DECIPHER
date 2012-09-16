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
          add2tbl = FALSE,
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
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\value{
A \code{data.frame} with the number of \code{bases}, \code{nonbases}, and \code{width} of each sequence.  The \code{width} is defined as the sum of \code{bases} and \code{nonbases} in each sequence.  The \code{row.names} of the \code{data.frame} correspond to the "row_names" in the \code{tblName} of the \code{dbFile}.
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