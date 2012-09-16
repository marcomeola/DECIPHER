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
         limit = -1,
         replaceChar = "-",
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
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.
}
  \item{limit}{
Number of results to display.  The default (\code{-1}) does not limit the number of results.
}
  \item{replaceChar}{
Optional character used to replace any characters of the sequence that are not present in the \code{DNA_ALPHABET}.
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
If RNA is present in the database then all U's are converted to T's before creating the \code{DNAStringSet}.
}
\value{
A \code{DNAStringSet} with the sequences that meet the specified criteria.  The \code{names} of the \code{DNAStringSet} correspond to the value in the "row_names" column of the database.
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