\name{FormGroups}
\alias{FormGroups}
\title{
Forms Groups By Rank
}
\description{
Agglomerates sequences into groups within a certain size range based on taxonomic rank.
}
\usage{
FormGroups(dbFile,
           tblName = "Seqs",
           goalSize = 1000,
           minGroupSize = 500,
           maxGroupSize = 10000,
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
  \item{goalSize}{
Number of sequences required in each group to stop adding more sequences.
}
  \item{minGroupSize}{
Minimum number of sequences in each group required to stop trying to recombine with a larger group.
}
  \item{maxGroupSize}{
Maximum number of sequences in each group allowed to continue agglomeration.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{verbose}{
Logical indicating whether to print database queries and other information.
}
}
\details{
\code{FormGroups} uses the ``rank'' field in the \code{dbFile} table to group sequences with similar taxonomic rank.  Rank information must be present in the \code{tblName}, such as that created by default when importing sequences from a GenBank formatted file.  The rank information must not contain repeated taxonomic names belonging to different lineages.

Beginning with the least common ranks, the algorithm agglomerates groups with similar ranks until the \code{goalSize} is reached.  If the group size is below \code{minGroupSize} then further agglomeration is attempted with a larger group.  If additional agglomeration results in a group larger than \code{maxGroupSize} then the agglomeration is undone so that the group is smaller.
}
\value{
A \code{data.frame} with the \code{rank} and corresponding identifier as \code{identifier}.  Note that quotes are stripped from identifiers to prevent problems that they may cause.  The \code{origin} gives the \code{rank} preceding the \code{identifier}.  If \code{add2tbl} is not \code{FALSE} then the ``identifier'' and ``origin'' columns are updated in \code{dbFile}.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{IdentifyByRank}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
g <- FormGroups(db, goalSize=10, minGroupSize=5, maxGroupSize=20)
head(g)
}