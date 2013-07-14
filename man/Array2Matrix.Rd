\name{Array2Matrix}
\alias{Array2Matrix}
\title{
Creates a Matrix Representation of a Microarray
}
\description{
Converts the output of \code{DesignArray} into the sparse matrix format used by \code{NNLS}.
}
\usage{
Array2Matrix(probes,
             verbose = TRUE)
}
\arguments{
  \item{probes}{
A set of microarray probes in the format output by \code{DesignArray}.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
A microarray can be represented by a matrix of hybridization efficiencies, where the rows represent each of the probes and the columns represent each the possible templates.  This matrix is sparse since microarray probes are designed to only target a small subset of the possible templates.
}
\value{
A list specifying the hybridization efficiency of each probe to its potential templates.
  \item{i}{
Element's row index in the sparse matrix.
}
  \item{j}{
Element's column index in the sparse matrix.
}
  \item{x}{
Non-zero elements' values representing hybridization efficiencies.
}
  \item{dimnames}{
A list of two components: the names of each probe, and the names of each template.
}
}
\references{
Coming Soon!
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{DesignArray}}, \code{\link{NNLS}}
}
\examples{
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
names(dna) <- 1:length(dna)
probes <- DesignArray(dna)
A <- Array2Matrix(probes)
}
