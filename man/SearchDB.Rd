\name{SearchDB}
\alias{SearchDB}
\title{
Obtain Specific Sequences from A Database
}
\description{
Returns the set of sequences meeting the search criteria.
}
\usage{
SearchDB(dbFile,
         tblName = "DNA",
         identifier = "",
         type = "DNAStringSet",
         limit = -1,
         replaceChar = "-",
         nameBy = "row_names",
         orderBy = "row_names",
         countOnly = FALSE,
         removeGaps = "none",
         verbose = TRUE,
         ...)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table where the sequences are located.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" (the default) then all identifiers are selected.
}
  \item{type}{
The type of \code{XStringSet} (sequences) to return.  This should be (an unambiguous abbreviation of) one of \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"}, \code{"BStringSet"}, \code{"QualityScaledDNAStringSet"}, \code{"QualityScaledRNAStringSet"}, \code{"QualityScaledAAStringSet"}, or \code{"QualityScaledBStringSet"}.
}
  \item{limit}{
Number of results to display.  The default (\code{-1}) does not limit the number of results.
}
  \item{replaceChar}{
Optional character used to replace any characters of the sequence that are not present in the \code{XStringSet}'s alphabet.  Not applicable if \code{type=="BStringSet"}.  (See details section below.)
}
  \item{nameBy}{
Character string giving the column name for naming the \code{XStringSet}.
}
  \item{orderBy}{
Character string giving the column name for sorting the results.  Defaults to the order of entries in the database.  Optionally can be followed by \code{" ASC"} or \code{" DESC"} to specify ascending (the default) or descending order.
}
  \item{countOnly}{
Logical specifying whether to return only the number of sequences.
}
  \item{removeGaps}{
Determines how gaps are removed in the sequences.  This should be (an unambiguous abbreviation of) one of \code{"none"}, \code{"all"} or \code{"common"}.
}
  \item{verbose}{
Logical indicating whether to display queries as they are sent to the database.
}
  \item{\dots}{
Additional expressions to add as part of a where clause in the query.  Further arguments provided in \dots will be added to the query separated by \code{" and "} as part of the where clause.
}
}
\details{
If \code{seqtype} is \code{"DNA"} then all U's are converted to T's before creating the \code{DNAStringSet}, and vise-versa if \code{seqtype} is \code{"RNAStringSet"}.  All remaining characters not in the \code{XStringSet}'s alphabet are converted to \code{replaceChar}.  Quality information is interpreted as \code{PredQuality} scores.
}
\value{
An \code{XStringSet} or \code{QualityScaledXStringSet} with the sequences that meet the specified criteria.  The \code{names} of the object correspond to the value in the \code{nameBy} column of the database.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{Seqs2DB}}, \code{\link{DB2FASTA}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db)
}