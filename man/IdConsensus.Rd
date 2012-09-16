\name{IdConsensus}
\alias{IdConsensus}
\title{
Create Consensus Sequences by Groups
}
\description{
Forms a consensus sequence representing the sequences in each group.
}
\usage{
IdConsensus(dbFile,
            tblName = "DNA",
            identifier = "",
            colName = "cluster",
            add2tbl = FALSE,
            verbose = TRUE,
            ...)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table in which to form consensus.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.
}
  \item{colName}{
Column containing the group name of each sequence.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Additional arguments to be passed directly to \code{ConsensusSequence}.
}
}
\details{
Creates a consensus sequence for each of the distinct groups defined in \code{colName}.  The resulting \code{DNAStringSet} contains as many consensus sequences as there are groups in \code{colName}.  For example, it is possible to create a set of consensus sequences with one consensus sequence for each \code{"id"} or \code{"cluster"}.
}
\value{
A \code{DNAStringSet} object containing the consensus sequence for each group.  The \code{names} of the \code{DNAStringSet} contain the number of sequences and name of each group.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{Seqs2DB}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
con <- IdConsensus(db, colName="id")
BrowseSequences(con, colorBases=TRUE)
}