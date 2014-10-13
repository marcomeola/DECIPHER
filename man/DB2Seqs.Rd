\name{DB2Seqs}
\alias{DB2Seqs}
\title{
Export Database Sequences to a FASTA or FASTQ File
}
\description{
Exports a database containing sequences to a FASTA or FASTQ formatted file of sequence records.
}
\usage{
DB2Seqs(file,
         dbFile,
         tblName = "DNA",
         identifier = "",
         type = "BStringSet",
         limit = -1,
         replaceChar = "-",
         nameBy = "description",
         orderBy = "row_names",
         removeGaps = "none",
         append = FALSE,
         width = 80,
         chunkSize = 1e5,
         clause = "",
         verbose = TRUE)
}
\arguments{
  \item{file}{
Character string giving the location where the file should be written.
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
  \item{type}{
The type of \code{XStringSet} (sequences) to export to a FASTA formatted file or \code{QualityScaledXStringSet} to export to a FASTQ formatted file.  This should be (an unambiguous abbreviation of) one of \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"}, \code{"BStringSet"}, \code{"QualityScaledDNAStringSet"}, \code{"QualityScaledRNAStringSet"}, \code{"QualityScaledAAStringSet"}, or \code{"QualityScaledBStringSet"}.  (See details section below.)
}
  \item{limit}{
Number of results to display.  The default (\code{-1}) does not limit the number of results.
}
  \item{replaceChar}{
Optional character used to replace any characters of the sequence that are not present in the \code{XStringSet}'s alphabet.  Not applicable if \code{type=="BStringSet"}.  (See details section below.)
}
  \item{nameBy}{
Character string giving the column name(s) for identifying each sequence record.  If more than one column name is provided, the information in each column is concatenated, separated by pairs of colons (``::''), in the order specified.
}
  \item{orderBy}{
Character string giving the column name for sorting the results.  Defaults to the order of entries in the database.  Optionally can be followed by \code{" ASC"} or \code{" DESC"} to specify ascending (the default) or descending order.
}
  \item{removeGaps}{
Determines how gaps ("-" or "." characters) are removed in the sequences.  This should be (an unambiguous abbreviation of) one of \code{"none"}, \code{"all"} or \code{"common"}.
}
  \item{append}{
Logical indicating whether to append the output to the existing \code{file}.
}
  \item{width}{
Integer specifying the maximum number of characters per line of sequence.  Not applicable when exporting to a FASTQ formatted file.
}
  \item{chunkSize}{
Number of lines of the \code{file} to write at a time.  Cannot be less than the total number of sequences if \code{removeGaps} is \code{"common"}.
}
  \item{clause}{
An optional character string to append to the query as a clause.
}
  \item{verbose}{
Logical indicating whether to display status.
}
}
\details{
Sequences are exported into either a FASTA or FASTQ file as determined by the \code{type} of sequences.  If \code{type} is an \code{XStringSet} then sequences are exported to FASTA format.  Quality information for \code{QualityScaledXStringSet}s are interpreted as \code{PredQuality} scores before export to FASTQ format.

If \code{type} is \code{"BStringSet"} (the default) then sequences are exported to a FASTA file exactly the same as they were when imported.  If \code{type} is \code{"DNAStringSet"} then all U's are converted to T's before export, and vise-versa if \code{type} is \code{"RNAStringSet"}.  All remaining characters not in the \code{XStringSet}'s alphabet are converted to \code{replaceChar}.
}
\value{
Writes a FASTA or FASTQ formatted file containing the sequence records in the database.

Returns the number of sequence records written to the \code{file}.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
tf <- tempfile()
DB2Seqs(tf, db, limit=10)
file.show(tf)
unlink(tf)
}