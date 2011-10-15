\name{Add2DB}
\alias{Add2DB}
\title{
Add Data To A Database
}
\description{
Adds a \code{data.frame} to a database table by \code{row.names}.
}
\usage{
Add2DB(myData,
       dbFile,
       tblName = "DNA",
       verbose = TRUE,
       ...)
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
  \item{verbose}{
Logical indicating whether to display each query as it is sent to the database.
}
  \item{\dots}{
Additional expressions to add as part of a where clause in the query.  Further arguments provided in \dots will be added to the query separated by \code{" and "} as part of the where clause.
}
}
\details{
Data contained in \code{myData} will be added to the \code{tblName} by its respective \code{row.names}.
}
\value{
Returns \code{TRUE} if the data was added successfully.
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