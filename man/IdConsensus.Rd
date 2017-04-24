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
            tblName = "Seqs",
            identifier = "",
            type = "DNAStringSet",
            colName = "identifier",
            processors = 1,
            verbose = TRUE,
            \dots)
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
  \item{type}{
The type of \code{XStringSet} (sequences) to use in forming consensus.  This should be (an abbreviation of) one of \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"}, or \code{"BStringSet"}.
}
  \item{colName}{
Column containing the group name of each sequence.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Additional arguments to be passed directly to \code{\link{ConsensusSequence}} for an \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet}, or to \code{consensusString} for a \code{BStringSet}.
}
}
\details{
Creates a consensus sequence for each of the distinct groups defined in \code{colName}.  The resulting \code{XStringSet} contains as many consensus sequences as there are distinct groups in \code{colName}.  For example, it is possible to create a set of consensus sequences with one consensus sequence for each \code{"id"} in the \code{tblName}.
}
\value{
An \code{XStringSet} object containing the consensus sequence for each group.  The \code{names} of the \code{XStringSet} contain the number of sequences and name of each group.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{Seqs2DB}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
con <- IdConsensus(db, colName="identifier", noConsensusChar="N")
BrowseSeqs(con)
}