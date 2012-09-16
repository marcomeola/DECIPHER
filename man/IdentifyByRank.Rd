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
               level = 3,
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
Level of the taxonomic rank.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{verbose}{
Logical indicating whether to print database queries and other information.
}
}
\details{
Simply identifies a sequence by a specific level of its taxonomic rank.  Requires that rank information be present in the \code{tblName}, such as that created when importing sequences from a GenBank file.

If the specified level of rank does not exist then the closest rank is chosen.  This makes it possible to determine the lowest level classification (e.g., genus) by specifying \code{level = 100}.
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