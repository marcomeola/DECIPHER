\name{AlignDB}
\alias{AlignDB}
\title{
Align Two Sets of Aligned Sequences in a Sequence Database
}
\description{
Merges the two separate sequence alignments in a database.  The aligned sequences must have separate identifiers in the same table or be located in different database tables.
}
\usage{
AlignDB(dbFile,
        tblName = "Seqs",
        identifier = "",
        type = "DNAStringSet",
        add2tbl = "Seqs",
        batchSize = 10000,
        perfectMatch = 5,
        misMatch = 0,
        gapOpening = -13,
        gapExtension = -1,
        gapPower = -1,
        terminalGap = -5,
        normPower = c(1, 0),
        substitutionMatrix = NULL,
        processors = 1,
        verbose = TRUE,
        \dots)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table(s) where the sequences are located.  If two \code{tblName}s are provided then the sequences in both tables will be aligned.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.  If two \code{identifier}s are provided then the set of sequences matching each \code{identifier} will be aligned.
}
  \item{type}{
The type of \code{XStringSet} being processed.  This should be (an abbreviation of) one of \code{"AAStringSet"}, \code{"DNAStringSet"}, or \code{"RNAStringSet"}.
}
  \item{add2tbl}{
Character string specifying the table name in which to add the aligned sequences.
}
  \item{batchSize}{
Integer specifying the number of sequences to process at a time.
}
  \item{perfectMatch}{
Numeric giving the reward for aligning two matching nucleotides in the alignment.  Only used when \code{type} is \code{DNAStringSet} or \code{RNAStringSet}.
}
  \item{misMatch}{
Numeric giving the cost for aligning two mismatched nucleotides in the alignment.  Only used when \code{type} is \code{DNAStringSet} or \code{RNAStringSet}.
}
  \item{gapOpening}{
Numeric giving the cost for opening a gap in the alignment.
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap in the alignment.
}
  \item{gapPower}{
Numeric specifying the exponent to use in the gap cost function.
}
  \item{terminalGap}{
Numeric giving the cost for allowing leading and trailing gaps ("-" or "." characters) in the alignment.  Either two numbers, the first for leading gaps and the second for trailing gaps, or a single number for both.
}
  \item{normPower}{
Numeric giving the exponent that controls the degree of normalization applied to scores by column occupancy.  If two numerics are provided, the first controls the normalization power of terminal gaps, while the second controls that of internal gaps.  A \code{normPower} of \code{0} does not normalize the scores, which results in all columns of the profiles being weighted equally, and is the optimal value for aligning fragmentary sequences.  A \code{normPower} of \code{1} normalizes the score for aligning two positions by their column occupancy (1 - fraction of gaps).  A \code{normPower} greater than \code{1} more strongly discourages aligning with ``gappy'' regions of the alignment.
}
  \item{substitutionMatrix}{
Either a substitution matrix representing the substitution scores for an alignment or the name of the amino acid substitution matrix to use in alignment.  The latter may be one of the following:  ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250'', or ``MIQS''.  The default (NULL) will use the \code{perfectMatch} and \code{misMatch} penalties for DNA/RNA or ``MIQS'' for AA.  (See examples section below.)
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{Codec}}.
}
}
\details{
Sometimes it is useful to align two large sets of sequences, where each set of sequences is already aligned but the two sets are not aligned to each other.  \code{AlignDB} first builds a profile of each sequence set in increments of \code{batchSize} so that the entire sequence set is not required to fit in memory.  Next the two profiles are aligned using dynamic programming.  Finally, the new alignment is applied to all the sequences as they are incrementally added to the \code{add2tbl}.

Two \code{identifier}s or \code{tblName}s must be provided, indicating the two sets of sequences to align.  The sequences corresponding to the first \code{identifier} and \code{tblName} will be aligned to those of the second \code{identifier} or \code{tblName}.  The aligned sequences are added to \code{add2tbl} under a new identifier formed from the concatenation of the two \code{identifier}s or \code{tblName}s.  (See examples section below.)
}
\value{
Returns the number of newly aligned sequences added to the database.
}
\references{
Wright, E. S. (2015). DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC Bioinformatics, 16, 322. http://doi.org/10.1186/s12859-015-0749-z
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignProfiles}}, \code{\link{AlignSeqs}}, \code{\link{AlignTranslation}}
}
\examples{
gen <- system.file("extdata", "Bacteria_175seqs.gen", package="DECIPHER")
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")

# Align two tables and place result into a third
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(gen, "GenBank", dbConn, "Seqs1", tblName="Set1")
Seqs2DB(fas, "FASTA", dbConn, "Seqs2", tblName="Set2")
AlignDB(dbConn, tblName=c("Set1", "Set2"), add2tbl="AlignedSets")
l <- IdLengths(dbConn, "AlignedSets", add2tbl=TRUE)
BrowseDB(dbConn, tblName="AlignedSets") # all sequences have the same width
dbDisconnect(dbConn)

# Align two identifiers and place the result in the same table
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(gen, "GenBank", dbConn, "Seqs1")
Seqs2DB(fas, "FASTA", dbConn, "Seqs2")
AlignDB(dbConn, identifier=c("Seqs1", "Seqs2"))
l <- IdLengths(dbConn, add2tbl=TRUE)
BrowseDB(dbConn) # note the sequences with a new identifier
dbDisconnect(dbConn)
}
