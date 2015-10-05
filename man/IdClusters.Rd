\name{IdClusters}
\alias{IdClusters}
\title{
Cluster Sequences By Distance or Sequence
}
\description{
Groups the sequences represented by a distance matrix into clusters of similarity.
}
\usage{
IdClusters(myDistMatrix = NULL,
           method = "UPGMA",
           cutoff = -Inf,
           showPlot = FALSE,
           asDendrogram = FALSE,
           myXStringSet = NULL,
           model = MODELS,
           processors = NULL,
           verbose = TRUE)
}
\arguments{
  \item{myDistMatrix}{
A symmetric \eqn{N} x \eqn{N} distance matrix with the values of dissimilarity between \eqn{N} sequences, or \code{NULL} if \code{method} is \code{"inexact"}.
}
  \item{method}{
An agglomeration method to be used.  This should be (an abbreviation of) one of \code{"complete"}, \code{"single"}, \code{"UPGMA"}, \code{"WPGMA"}, \code{"NJ"}, \code{"ML"}, or \code{"inexact"}.  (See details section below.)
}
  \item{cutoff}{
A vector with the maximum edge length separating the sequences in the same cluster.  Multiple cutoffs may be provided in ascending or descending order.  If \code{asDendrogram=TRUE} or \code{showPlot=TRUE} then only one cutoff may be specified.  (See details section below.)
}
  \item{showPlot}{
Logical specifying whether or not to plot the resulting dendrogram.  Not applicable if \code{method='inexact'}.
}
  \item{asDendrogram}{
Logical.  If \code{TRUE} then the object returned is of class \code{dendrogram}.  Not applicable if \code{method='inexact'}.
}
  \item{myXStringSet}{
If \code{method} is \code{"ML"}, the \code{DNAStringSet} or \code{RNAStringSet} used in the creation of \code{myDistMatrix}.  If \code{method} is \code{"inexact"}, the \code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet} to cluster.  Not applicable for other methods.
}
  \item{model}{
One or more of the available \code{MODELS} of DNA evolution.  Only applicable if \code{method} is \code{"ML"}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{IdClusters} groups the input sequences into clusters using a set dissimilarities representing the distance between \eqn{N} sequences.  Initially a phylogenetic tree is formed using the specified \code{method}.  Then each leaf (sequence) of the tree is assigned to a cluster based on its edge lengths to the other sequences.  The available clustering methods are described as follows:

Ultrametric methods:  The method \code{complete} assigns clusters using complete-linkage so that sequences in the same cluster are no more than \code{cutoff} percent apart.  The method \code{single} assigns clusters using single-linkage so that sequences in the same cluster are within \code{cutoff} of at least one other sequence in the same cluster.  \code{UPGMA} (the default) or \code{WPGMA} assign clusters using average-linkage which is a compromise between the sensitivity of complete-linkage clustering to outliers and the tendency of single-linkage clustering to connect distant relatives that do not appear to be closely related.  \code{UPGMA} produces an unweighted tree, where each leaf contributes equally to the average edge lengths, whereas \code{WPGMA} produces a weighted result.

Additive methods:  \code{NJ} uses the Neighbor-Joining method proposed by Saitou and Nei that does not assume lineages evolve at the same rate (the molecular clock hypothesis).  The \code{NJ} method is typically the most phylogenetically accurate of the above distance-based methods.  \code{ML} creates a neighbor-joining tree and then iteratively maximizes the likelihood of the tree given the aligned sequences (\code{myXStringSet}).  This is accomplished through a combination of optimizing edge lengths with Brent's method and improving tree topology with nearest-neighbor interchanges (NNIs).  When \code{method="ML"}, one or more \code{MODELS} of DNA evolution must be specified.  Model parameters are iteratively optimized to maximize likelihood, except base frequencies which are empirically determined.  If multiple \code{model}s are given, the best \code{model} is automatically chosen based on BIC calculated from the likelihood and the sample size (defined as the number of variable sites in the DNA sequence).

Sequence-only method:  \code{inexact} uses a greedy incremental algorithm to directly assign sequences to clusters without a distance matrix.  First the sequences are ordered by length and the longest sequence becomes the first cluster seed.  If the second sequence is less than \code{cutoff} percent distance then it is added to the cluster, otherwise it becomes a new cluster representative.  The remaining sequences are matched to cluster representatives using an ordered k-mer strategy, and then compared to the top hits with pairwise alignment.  This approach, finding the closest cluster representatives followed by pairwise alignment to obtain the percent identity, is repeated until all sequences belong to a cluster.  This process results in clusters with members generally separated by less than \code{cutoff} distance.

Cutoffs may be provided in ascending or descending order, except for \code{inexact} clustering, which requires descending cutoffs.  If multiple cutoffs are provided in \emph{descending} order then clustering at each new value of \code{cutoff} is continued within the prior \code{cutoff}'s clusters.  In this way clusters at lower values of \code{cutoff} are completely contained within their umbrella clusters at higher values of \code{cutoff}.  This is useful for defining taxonomy, where lower level groups (e.g., genera) are expected not to straddle multiple higher level groups (e.g., families).  If multiple cutoffs are provided in \emph{ascending} order then clustering at each level of \code{cutoff} is independent of the prior level.  This may result in fewer clusters for \code{NJ} and \code{ML} methods, but will have no impact on ultrametric methods.
}
\value{
If \code{asDendrogram=FALSE} (the default), then a data.frame is returned with a column for each cutoff specified.  This data.frame has dimensions \eqn{N*M}, where each one of \eqn{N} sequences is assigned to a cluster at the \eqn{M}-level of cutoff.  The row.names of the data.frame correspond to the dimnames of \code{myDistMatrix}.
If \code{asDendrogram=TRUE}, returns an object of class \code{dendrogram} that can be used for further manipulation and plotting.  Leaves of the dendrogram are randomly colored by cluster number.
}
\references{
Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17(6)}, 368-376.

Ghodsi, M., Liu, B., & Pop, M. (2011) DNACLUST. \emph{BMC Bioinformatics}, \bold{12(1)}, 271. doi:10.1186/1471-2105-12-271.

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