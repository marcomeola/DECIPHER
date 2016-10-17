\name{FindChimeras}
\alias{FindChimeras}
\title{
Find Chimeras in a Sequence Database
}
\description{
Finds chimeras present in a database of sequences.  Makes use of a reference database of (presumed to be) good quality sequences.
}
\usage{
FindChimeras(dbFile,
             tblName = "Seqs",
             identifier = "",
             dbFileReference,
             tblNameReference = "Seqs",
             batchSize = 100,
             minNumFragments = 20000,
             tb.width = 5,
             multiplier = 20,
             minLength = 70,
             minCoverage = 0.6,
             overlap = 100,
             minSuspectFragments = 6,
             showPercentCoverage = FALSE,
             add2tbl = FALSE,
             maxGroupSize = -1,
             minGroupSize = 100,
             excludeIDs = NULL,
             processors = 1,
             verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file to be checked for chimeric sequences.
}
  \item{tblName}{
Character string specifying the table in which to check for chimeras.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.
}
  \item{dbFileReference}{
A SQLite connection object or a character string specifying the path to the reference database file of (presumed to be) good quality sequences.  A 16S reference database is available from \url{http://DECIPHER.codes}.
}
  \item{tblNameReference}{
Character string specifying the table with reference sequences.
}
  \item{batchSize}{
Number sequences to tile with fragments at a time.
}
  \item{minNumFragments}{
Number of suspect fragments to accumulate before searching through other groups.
}
  \item{tb.width}{
A single integer [1..14] giving the number of nucleotides at the start of each fragment that are part of the trusted band.
}
  \item{multiplier}{
A single integer specifying the multiple of fragments found out-of-group greater than fragments found in-group in order to consider a sequence a chimera.
}
  \item{minLength}{
Minimum length of a chimeric region in order to be considered as a chimera.
}
  \item{minCoverage}{
Minimum fraction of coverage necessary in a chimeric region.
}
  \item{overlap}{
Number of nucleotides at the end of the sequence that the chimeric region must overlap in order to be considered a chimera.
}
  \item{minSuspectFragments}{
Minimum number of suspect fragments belonging to another group required to consider a sequence a chimera.
}
  \item{showPercentCoverage}{
Logical indicating whether to list the percent coverage of suspect fragments in each chimeric region in the output.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{maxGroupSize}{
Maximum number of sequences searched in a group.  A value of less than 0 means the search is unlimited.
}
  \item{minGroupSize}{
The minimum number of sequences in a group to be considered as part of the search for chimeras.  May need to be set to a small value for reference database with mostly small groups.
}
  \item{excludeIDs}{
Optional character vector of \code{identifier}(s) to exclude from database searches, or \code{NULL} (the default) to not exclude any.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{FindChimeras} works by finding suspect fragments that are uncommon in the group where the sequence belongs, but very common in another group where the sequence does not belong.  Each sequence in the \code{dbFile} is tiled into short sequence segments called fragments.  If the fragments are infrequent in their respective group in the \code{dbFileReference} then they are considered suspect.  If enough suspect fragments from a sequence meet the specified constraints then the sequence is flagged as a chimera.

The default parameters are optimized for full-length 16S sequences (> 1,000 nucleotides).  Shorter 16S sequences require two parameters that are different than the defaults:  \code{minLength = 40}, and \code{minSuspectFragments = 2}.

Groups are determined by the identifier present in each database.  For this reason, the groups in the \code{dbFile} should exist in the groups of the \code{dbFileReference}.  The reference database is assumed to contain many sequences of only good quality.

If a reference database is not present then it is feasible to create a reference database by using the input database as the reference database.  Removing chimeras from the reference database and then iteratively repeating the process can result in a clean reference database.

For non-16S sequences it may be necessary to optimize the parameters for the particular sequences.  The simplest way to perform an optimization is to experiment with different input parameters on artificial chimeras such as those created using \code{\link{CreateChimeras}}.  Adjusting input parameters until the maximum number of artificial chimeras are identified is the easiest way to determine new defaults.
}
\value{
A \code{data.frame} containing only the sequences that meet the specifications for being chimeric.  The chimera column contains information on the chimeric region and to which group it belongs.  The \code{row.names} of the \code{data.frame} correspond to those of the sequences in the \code{dbFile}.
}
\references{
ES Wright et al. (2012) "DECIPHER: A Search-Based Approach to Chimera Identification for 16S rRNA Sequences." Applied and Environmental Microbiology, doi:10.1128/AEM.06516-11.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{CreateChimeras}}, \code{\link{Add2DB}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
# It is necessary to set dbFileReference to the file path of the
# 16S reference database available from http://DECIPHER.codes
chimeras <- FindChimeras(db, dbFileReference=db)
}