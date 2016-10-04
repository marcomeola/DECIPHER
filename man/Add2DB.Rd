\name{Add2DB}
\alias{Add2DB}
\title{
Add Data to a Database
}
\description{
Adds a \code{data.frame} to a database table by its \code{row.names}.
}
\usage{
Add2DB(myData,
       dbFile,
       tblName = "Seqs",
       clause = "",
       verbose = TRUE)
}
\arguments{
  \item{myData}{
Data frame containing information to be added to the \code{dbFile}.
}
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table in which to add the data.
}
  \item{clause}{
An optional character string to append to the query as part of a ``where clause''.
}
  \item{verbose}{
Logical indicating whether to display each query as it is sent to the database.
}
}
\details{
Data contained in \code{myData} will be added to the \code{tblName} by its respective \code{row.names}.
}
\value{
Returns \code{TRUE} if the data was added successfully, or \code{FALSE} otherwise.
}
\references{
ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{Seqs2DB}}, \code{\link{SearchDB}}, \code{\link{BrowseDB}}
}
\examples{
# Create a sequence database
gen <- system.file("extdata", "Bacteria_175seqs.gen", package="DECIPHER")
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(gen, "GenBank", dbConn, "Bacteria")

# Identify the sequence lengths
l <- IdLengths(dbConn)

# Add lengths to the database
Add2DB(l, dbConn)

# View the added lengths
BrowseDB(dbConn)
dbDisconnect(dbConn)
}