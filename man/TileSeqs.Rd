\name{TileSeqs}
\alias{TileSeqs}
\title{
Form a Set of Tiles for Each Group of Sequences.
}
\description{
Creates a set of tiles that represent each group of sequences in the database for downstream applications.
}
\usage{
TileSeqs(dbFile,
         tblName = "DNA",
         identifier = "",
         minLength = 26,
         maxLength = 27,
         maxTilePermutations = 10,
         minCoverage = 0.9,
         add2tbl = FALSE,
         verbose = TRUE,
         ...)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table of sequences to use for forming tiles.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.
}
  \item{minLength}{
Integer providing the minimum number of nucleotides in each tile.  Typically the same  or slightly less than \code{maxLength}.
}
  \item{maxLength}{
Integer providing the maximum number of nucleotides in each tile.  Tiles are designed primarily for this length, which should ideally be slightly greater than the maximum length of oligos used in downstream functions.
}
  \item{maxTilePermutations}{
Integer specifying the maximum number of tiles in each target site.
}
  \item{minCoverage}{
Numeric providing the fraction of coverage that is desired for each target site in the group.  For example, a \code{minCoverage} of 0.9 request that additional tiles are added until 90% of the group is represented by the tile permutations.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Additional arguments to be passed directly to \code{SearchDB}.
}
}
\details{
This function will create a set of overlapping tiles representing each target site in an alignment of sequences.  The most common tile permutations are added until the minimum group coverage is obtained.

Target sites with one more more tiles not meeting a set of requirements are marked with \code{misprime} equals \code{TRUE}.  Requirements are a minimum group coverage, minimum length, and a maximum length.  Additionally, tiles are required not to contain more than four runs of a single base or four di-nucleotide repeats.
}
\value{
A \code{data.frame} with a row for each tile, and multiple columns of information.  The \code{row_names} column gives the row number.  The \code{start}, \code{end}, \code{start_aligned}, and \code{end_aligned} columns provide positioning of the tile in a consensus sequence formed from the group.  The column \code{misprime} is a logical specifying whether the tile meets the specified constraints.  The columns \code{width} and \code{id} indicate the tiles length and group of origin, respectively.

The \code{coverage} field gives the fraction of sequences containing the tile in the group that encompass the tiles start and end positions in the alignment, whereas the \code{groupCoverage} contains the fraction of all sequences in the group containing a tile at their respective target site.  For example, if 10% of sequences have a certain tile permutation in the first alignment position, but only a single sequence includes this position in the alignment then \code{coverage} would be 100%, while \code{groupCoverage} would be 10%.

The final column, \code{target_site}, provides the sequence of the tile.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\note{
If \code{add2tbl} is \code{TRUE} then the tiles will be added to the database table that currently contains the sequences used for tiling.  The added tiles may cause interference when querying a table of sequences.  Therefore, it is recommended to add the tiles to their own table, for example, by using \code{add2tbl="Tiles"}.
}
\seealso{
\code{\link{DesignPrimers}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
tiles <- TileSeqs(db, identifier="Pseudomonas")
}