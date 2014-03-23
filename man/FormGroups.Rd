\name{FormGroups}
\alias{FormGroups}
\title{
Forms Groups By Rank
}
\description{
Agglomerates sequences into groups in a certain size range based on taxonomic rank.
}
\usage{
FormGroups(dbFile,
           tblName = "DNA",
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
\code{FormGroups} uses the ``rank'' field in the \code{dbFile} table to group sequences with similar taxonomic rank.  Requires that rank information be present in the \code{tblName}, such as that created when importing sequences from a GenBank file.  The rank information must not contain repeated names belonging to different lineages.

Beginning with the least common ranks, the algorithm agglomerates groups with similar ranks until the \code{goalSize} is reached.  If the group size is below \code{minGroupSize} then further agglomeration is attempted with a larger group.  If additional agglomeration results in a group larger than \code{maxGroupSize} then the agglomeration is undone so that the group is smaller.
}
\value{
Returns a \code{data.frame} of rank and id for each group.  If \code{add2tbl} is not \code{FALSE} then the \code{tblName} is updated with the group as the identifier.
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
}