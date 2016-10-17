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
               tblName = "Seqs",
               level = 0,
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
\code{IdentifyByRank} simply identifies a sequence by a specific level of its taxonomic rank.  Requires that \code{rank} information be present in the \code{tblName}, such as that created by default when importing sequences from a GenBank formatted file.

The input parameter \code{level} should be an integer giving the ``level'' of the taxonomic rank to choose as the identifier.  Negative \code{level}s are interpreted as being that many levels from the last level in each rank.  The \code{level} zero selects the base level (see below).

If the specified level of rank does not exist then the closest rank is chosen.  Therefore, setting \code{level} to \code{Inf} will always select the last taxonomic level (i.e., genus).

For example, a representative ``rank'' imported from a GenBank file is:\cr
Saccharomyces cerevisiae\cr
Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;\cr
Saccharomycetales; Saccharomycetaceae; Saccharomyces.

Setting \code{level} to \code{0} would result in an \code{identifier} of ``Saccharomyces cerevisiae'', because it is on the first line.  A \code{level} of \code{2} would return ``Fungi'', and \code{-2} (second to last) would return ``Saccharomycetaceae''.  A \code{level} of \code{Inf} would find the nearest level to the end, ``Saccharomyces''.
}
\value{
A \code{data.frame} with the \code{rank} and corresponding identifier as \code{identifier}.  Note that quotes are stripped from identifiers to prevent problems that they may cause.  The \code{origin} gives the \code{rank} preceding the \code{identifier}.  If \code{add2tbl} is not \code{FALSE} then the ``identifier'' column is updated in \code{dbFile}.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{FormGroups}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
ids <- IdentifyByRank(db, level=Inf)
head(ids)
}