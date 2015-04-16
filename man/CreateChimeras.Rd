\name{CreateChimeras}
\alias{CreateChimeras}
\title{
Create Artificial Chimeras
}
\description{
Creates artificial random chimeras from a set of sequences.
}
\usage{
CreateChimeras(myDNAStringSet,
               numChimeras = 10,
               numParts = 2,
               minLength = 80,
               maxLength = Inf,
               minChimericRegionLength = 30,
               randomLengths = TRUE,
               includeParents = TRUE,
               processors = NULL,
               verbose = TRUE)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object with aligned sequences.
}
  \item{numChimeras}{
Number of chimeras desired.
}
  \item{numParts}{
Number of chimeric parts from which to form a single chimeric sequence.
}
  \item{minLength}{
Minimum length of the complete chimeric sequence.
}
  \item{maxLength}{
Maximum length of the complete chimeric sequence.
}
  \item{minChimericRegionLength}{
Minimum length of the chimeric region of each sequence part.
}
  \item{randomLengths}{
Logical specifying whether to create random length chimeras in addition to random breakpoints.
}
  \item{includeParents}{
Whether to include the parents of each chimera in the output.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Forms a set of random chimeras from the input set of (typically good quality) sequences.  The chimeras are created by merging random sequences at random breakpoints.  These chimeras can be used for testing the accuracy of the \code{\link{FindChimeras}} or other chimera finding functions.
}
\value{
A \code{DNAStringSet} object containing chimeras.  The \code{names} of the chimeras are specified as "parent #1 name [chimeric region] (distance from parent to chimera), ...".

If \code{includeParents = TRUE} then the parents of the chimeras are included at the end of the result.  The parents are trimmed to the same length as the chimera if \code{randomLengths = TRUE}.  The \code{names} of the parents are specified as "parent #1 name [region] (distance to parent #2, ...)".
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{FindChimeras}}, \code{\link{Seqs2DB}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db)
chims <- CreateChimeras(dna)
BrowseSeqs(chims)
}