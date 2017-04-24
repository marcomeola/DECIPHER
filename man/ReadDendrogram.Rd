\name{ReadDendrogram}
\alias{ReadDendrogram}
\title{
Read a Dendrogram from a Newick Formatted File
}
\description{
Reads a dendrogram object from a file in Newick (also known as New Hampshire) parenthetic format.
}
\usage{
ReadDendrogram(file,
               convertBlanks = TRUE,
               internalLabels = TRUE,
               keepRoot = TRUE)
}
\arguments{
  \item{file}{
a connection object or a character string.
}
  \item{convertBlanks}{
Logical specifying whether to convert underscores in unquoted leaf labels to spaces.
}
  \item{internalLabels}{
Logical indicating whether to keep internal node labels as ``edgetext'' preceding the node in the \code{dendrogram}.
}
  \item{keepRoot}{
Logical specifying whether to keep the root node (if one is present) as a dendrogram leaf.
}
}
\details{
\code{ReadDendrogram} will create a dendrogram object from a Newick formatted tree.  Note that all edge lengths must be specified, but labels are optional.
}
\value{
An object of class \code{dendrogram}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{IdClusters}}, \code{\link{WriteDendrogram}}
}
\examples{
tf <- tempfile()
dists <- matrix(c(0, 10, 20, 10, 0, 5, 20, 5, 0),
    nrow=3,
    dimnames=list(c("dog", "elephant", "horse")))
dend1 <- IdClusters(dists, method="NJ", type="dendrogram")
WriteDendrogram(dend1, file=tf)

dend2 <- ReadDendrogram(tf)
layout(matrix(1:2))
plot(dend1, main="Dendrogram Written")
plot(dend2, main="Dendrogram Read")

# Note that the ordering information is lost
any(unlist(dend1) != unlist(dend2)) # TRUE

unlink(tf)
}
