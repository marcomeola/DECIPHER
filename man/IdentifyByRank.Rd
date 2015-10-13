\name{IdentifyByRank}
\alias{IdentifyByRank}
\title{
Identify By Taxonomic Rank
}
\description{
Identifies sequences by a specific level of their taxonomic rank.
}
\usage{
IdentifyByRank(dbFile,
               tblName = "DNA",
               level = Inf,
               add2tbl = FALSE,
               verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table where the rank information is located.
}
  \item{level}{
Level of the taxonomic rank.  (See details section below.)
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{verbose}{
Logical indicating whether to print database queries and other information.
}
}
\details{
\code{IdentifyByRank} simply identifies a sequence by a specific level of its taxonomic rank.  Requires that \code{rank} information be present in the \code{tblName}, such as that created when importing sequences from a GenBank formatted file.

The input parameter \code{level} should be a non-zero integer giving the ``level'' of the taxonomic rank to choose as the identifier.  Negative \code{level}s are interpreted as that being many levels from the last level in each rank.

If the specified level of rank does not exist then the closest rank is chosen.  Therefore, the default \code{level} (\code{Inf}) will always select the last taxonomic level (e.g., species).
}
\value{
A \code{data.frame} with the \code{rank} and corresponding identifier as \code{"id"}.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{FormGroups}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
ids <- IdentifyByRank(db)
}