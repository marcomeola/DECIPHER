\name{Seqs2DB}
\alias{Seqs2DB}
\title{
Add Sequences from Text File to Database
}
\description{
Adds sequences to a database.
}
\usage{
Seqs2DB(seqs,
        type,
        dbFile,
        identifier,
        tblName = "DNA",
        chunkSize = 99999,
        replaceTbl = FALSE,
        verbose = TRUE)
}
\arguments{
  \item{seqs}{
Either a character string specifying the file path to the file containing the sequences, or a \code{DNAStringSet} object.
}
  \item{type}{
The type of sequences being imported.  This should be (an unambiguous abbreviation of) one of \code{"FASTA"}, \code{"GenBank"}, or \code{"DNAStringSet"}.
}
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.  If the \code{dbFile} does not exist then a new database is created at this location.
}
  \item{identifier}{
Character string specifying the \code{"id"} to give the imported sequences in the database.
}
  \item{tblName}{
Character string specifying the table in which to add the sequences.
}
  \item{chunkSize}{
Number of lines of the \code{seqs} to read at a time.  For very large sequence files, using \code{1e7} results in a quicker import than the default (\code{99999}), but only if enough memory is available.
}
  \item{replaceTbl}{
Logical.  If \code{FALSE} (the default) then the sequences are appended to any already existing in the table.  If \code{TRUE} then any sequences already in the table are overwritten.
}
  \item{verbose}{
Logical indicating whether to display each query as it is sent to the database.
}
}
\details{
Sequences are imported into the database in chunks of lines specified by \code{chunkSize}.  The sequences can then be identified by searching the database for the \code{identifier} provided.  Sequences are added to the database verbatim, so that no sequence information is lost when the sequences are exported from the database.
}
\value{
The total number of sequences in the database table is returned after import.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{SearchDB}}, \code{\link{DB2FASTA}}
}
\examples{
gen <- system.file("extdata", "Bacteria_175seqs.gen", package="DECIPHER")
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(gen, "GenBank", dbConn, "Bacteria")
BrowseDB(dbConn)
dbDisconnect(dbConn)
}