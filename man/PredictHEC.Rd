\name{PredictHEC}
\alias{PredictHEC}
\title{
Predict Protein Secondary Structure as Helix, Beta-Sheet, or Coil
}
\description{
Predicts 3-state protein secondary structure based on the primary (amino acid) sequence using the GOR IV method (Garnier et al., 1996).
}
\usage{
PredictHEC(myAAStringSet,
           type = "states",
           windowSize = 7,
           background = c(H = -0.12, E = -0.25, C = 0.23),
           HEC_MI1 = NULL,
           HEC_MI2 = NULL)
}
\arguments{
  \item{myAAStringSet}{
An \code{AAStringSet} object of sequences.
}
  \item{type}{
Character string indicating the type of results desired.  This should be (an unambiguous abbreviation of) one of \code{"states"}, \code{"scores"}, or \code{"probabilities"}.
}
  \item{windowSize}{
Numeric specifying the number of residues to the left or right of the center position to use in the prediction.
}
  \item{background}{
Numeric vector with the background ``scores'' for each of the three states (H, E, and C).
}
  \item{HEC_MI1}{
An array of dimensions 20 x 21 x 3 giving the mutual information for single residues.
}
  \item{HEC_MI2}{
An array of dimensions 20 x 20 x 21 x 21 x 3 giving the mutual information for pairs of residues.
}
}
\details{
The GOR (Garnier-Osguthorpe-Robson) method is an information-theory method for prediction of secondary structure based on the primary sequence of a protein.  Version IV of the method makes 3-state predictions based on the mutual information contained in single residues and pairs of residues within \code{windowSize} residues of the position being assigned.  This approach is about 65\% accurate, and is one of the most accurate methods for assigning secondary structure that only use a single sequence.  This implementation of GOR IV does not use decision constants or the number of contiguous states when assigning the final state.  Note that characters other than the standard 20 amino acids are not assigned a state.
}
\value{
If \code{type} is \code{"states"} (the default), then the output is a character vector with the secondary structure assignment ("H", "E", or "C") for each residue in \code{myAAStringSet}.

Otherwise, the output is a list with one element for each sequence in \code{myAAStringSet}.  Each list element contains a matrix of dimension 3 (H, E, or C) by the number of residues in the sequence.  If \code{type} is \code{"scores"}, then values in the matrix represent log-odds ``scores''.  If \code{type} is \code{"probabilities"} then the values represent the normalized probabilities of the three states at a position.
}
\references{
Garnier, J., Gibrat, J. F., & Robson, B. (1996). GOR method for predicting protein secondary structure from amino acid sequence. \emph{Methods in Enzymology}, \bold{266}, 540-553.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{HEC_MI1}}, \code{\link{HEC_MI2}}, \code{\link{PredictDBN}}
}
\examples{
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
aa <- translate(dna)
hec <- PredictHEC(aa)
head(hec)
}