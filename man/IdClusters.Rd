\name{IdClusters}
\alias{IdClusters}
\title{
Cluster Sequences By Distance
}
\description{
Groups the sequences represented by a distance matrix into clusters of similarity.
}
\usage{
IdClusters(myDistMatrix,
           method = "UPGMA",
           cutoff = -Inf,
           showPlot = FALSE,
           asDendrogram = FALSE,
           myDNAStringSet = NULL,
           model = MODELS,
           add2tbl = FALSE,
           dbFile = NULL,
           processors = NULL,
           verbose = TRUE)
}
\arguments{
  \item{myDistMatrix}{
A symmetric \eqn{N} x \eqn{N} distance matrix with the values of dissimilarity between \eqn{N} sequences.
}
  \item{method}{
An agglomeration method to be used.  This should be (an unambiguous abbreviation of) one of \code{"complete"}, \code{"single"}, \code{"UPGMA"}, \code{"average"}, \code{"NJ"} or \code{"ML"}.  (See details section below.)
}
  \item{cutoff}{
A vector with the maximum edge length separating the sequences in the same cluster.  If \code{asDendrogram=TRUE} then only one cutoff may be specified.
}
  \item{showPlot}{
Logical specifying whether or not to plot the resulting dendrogram.
}
  \item{asDendrogram}{
Logical.  If \code{TRUE} the object returned is of class \code{dendrogram}.
}
  \item{myDNAStringSet}{
\code{DNAStringSet} used in the creation of \code{myDistMatrix}.  Only necessary if \code{method} is \code{"ML"}.
}
  \item{model}{
One or more of the available \code{MODELS} of DNA evolution.  Only necessary if \code{method} is \code{"ML"}.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{dbFile}{
A connection to a SQLite database or character string giving the path to the database file.  Only necessary if \code{add2tbl} is not \code{FALSE}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Groups the input sequences into clusters using a set dissimilarities representing the distance between \eqn{N} sequences.  Initially a phylogenetic tree is formed using the specified \code{method}.  Then each leaf (sequence) of the tree is assigned to a cluster based on its edge lengths to the other sequences.

A number of different clustering methods are provided.  The method (\code{complete} assigns clusters using complete-linkage so that sequences in the same cluster are no more than \code{cutoff} percent apart.  The method \code{single} assigns clusters using single-linkage so that sequences in the same cluster are within \code{cutoff} of at least one other sequence in the same cluster.  \code{UPGMA} or \code{average} (the default) assigns clusters using average-linkage which is a compromise between the sensitivity of complete-linkage clustering to outliers and the tendency of single-linkage clustering to connect distant relatives that do not appear to be closely related.

\code{NJ} uses the Neighbor-Joining method proposed by Saitou and Nei that does not assume lineages evolve at the same rate (the molecular clock hypothesis).  The \code{NJ} method is typically the most phylogenetically accurate of the above distance-based methods.  \code{ML} creates a neighbor-joining tree and then iteratively maximizes the likelihood of the tree given the aligned sequences (\code{myDNAStringSet}).  This is accomplished through a combination of optimizing edge lengths with Brent's method and improving tree topology with nearest-neighbor interchanges (NNIs).  When \code{method="ML"}, one or more \code{MODELS} of DNA evolution must be specified.  If multiple \code{model}s are given their parameters are optimized and the best \code{model} is automatically chosen based on BIC using the sample size defined as the number of variable sites in the DNA sequence.

If a \code{add2tbl=TRUE} then the resulting data.frame is added/updated into column(s) of the default table "DNA" in \code{dbFile}.  If \code{add2tbl} is a character string then the result is added to the specified table name in \code{dbFile}.  The added/updated column names are printed if \code{verbose=TRUE}.
}
\value{
If \code{asDendrogram=FALSE} (the default), returns a data.frame with a column for each cutoff specified.  The row.names of the data.frame correspond to the dimnames of \code{myDistMatrix}.  Each one of \eqn{N} sequences is assigned to one of \eqn{M} clusters.
If \code{asDendrogram=TRUE}, returns an object of class \code{dendrogram} that can be used for further manipulation and plotting.  Leaves of the dendrogram are randomly colored by cluster number.
}
\references{
Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17(6)}, 368-376

Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. \emph{Molecular Biology and Evolution}, \bold{4(4)}, 406-425.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{DistanceMatrix}}, \code{\link{Add2DB}}, \code{\link{MODELS}}
}
\examples{
# using the matrix from the original paper by Saitou and Nei
m <- matrix(0,8,8)
m[2:8,1] <- c(7, 8, 11, 13, 16, 13, 17)
m[3:8,2] <- c(5, 8, 10, 13, 10, 14)
m[4:8,3] <- c(5, 7, 10, 7, 11)
m[5:8,4] <- c(8, 11, 8, 12)
m[6:8,5] <- c(5, 6, 10)
m[7:8,6] <- c(9, 13)
m[8,7] <- c(8)

# returns an object of class "dendrogram"
myClusters <- IdClusters(m, cutoff=10, method="NJ", showPlot=TRUE, asDendrogram=TRUE)

# example of specifying a cutoff
# returns a data frame
IdClusters(m, cutoff=c(2,6,10,20))
}