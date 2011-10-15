\name{DB2FASTA}
\alias{DB2FASTA}
\title{
Export Database Sequences to FASTA File
}
\description{
Exports a database containing sequences to a FASTA formatted file of sequences.
}
\usage{
DB2FASTA(file,
         dbFile,
         tblName = "DNA",
         identifier = "",
         limit = -1,
         replaceChar = NULL,
         orderBy = "row_names",
         append = FALSE,
         comments = TRUE,
         removeGaps = "none",
         verbose = TRUE,
         ...)
}
\arguments{
  \item{file}{
Character string giving the location where the FASTA file should be written.
}
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table in which to extract the data.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.
}
  \item{limit}{
Number of results to display.  The default (\code{-1}) does not limit the number of results.
}
  \item{replaceChar}{
Optional character used to replace any sequence characters not present in the \code{DNA_ALPHABET}.  If \code{NULL} (the default) then no replacement occurs and the sequences are exported identical to how they were upon import.
}
  \item{orderBy}{
Character string giving the column name for sorting the results.  Defaults to the order of entries in the database.  Optionally can be followed by \code{" ASC"} or \code{" DESC"} to specify ascending (the default) or descending order.
}
  \item{append}{
Logical indicating whether to append the results to the existing \code{file}.
}
  \item{comments}{
Logical specifying whether to add the value of any database fields into the FASTA record description separated by semicolons.
}
  \item{removeGaps}{
Determines how gaps are removed in the sequences.  This should be (an unambiguous abbreviation of) one of \code{"none"}, \code{"all"} or \code{"common"}.
}
  \item{verbose}{
Logical indicating whether to display status.
}
  \item{\dots}{
Additional expressions to add as part of a where clause in the query.  Further arguments provided in \dots will be added to the query separated by \code{" and "} as part of the where clause.
}
}
\value{
Writes a FASTA formatted file containing the sequences in the database.

Returns \code{TRUE} if the file was written successfully.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
tf <- tempfile()
DB2FASTA(tf, db, l=10)
file.show(tf)
unlink(tf)
}