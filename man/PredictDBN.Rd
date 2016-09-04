\name{PredictDBN}
\alias{PredictDBN}
\title{
Predict RNA Secondary Structure in Dot-Bracket Notation
}
\description{
Predicts a consensus RNA secondary structure from a multiple sequence alignment using mutual information.
}
\usage{
PredictDBN(myXStringSet,
           type = "states",
           minOccupancy = 0.5,
           impact = c(1, 1.2, 0.4, -1),
           avgProdCorr = 1,
           slope = 2,
           shift = 1.3,
           threshold = 0.5,
           pseudoknots = 1,
           weight = 1,
           processors = 1,
           verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} object containing aligned sequences.
}
  \item{type}{
Character string indicating the type of results desired.  This should be (an unambiguous abbreviation of) one of \code{"states"}, \code{"pairs"}, \code{"scores"}, or \code{"structures"}.
}
  \item{minOccupancy}{
Numeric specifying the minimum occupancy (1 - fraction of gaps) required to include a column of the alignment in the prediction.
}
  \item{impact}{
A vector with four elements giving the weights of A/U, G/C, G/U, and other pairings, respectively.  The last element of \code{impact} is the penalty for pairings that are inconsistent with two positions being paired (e.g., A/- or A/C).
}
  \item{avgProdCorr}{
Numeric specifying the weight of the average product correction (APC) term, as described in Buslje et al. (2009).
}
  \item{slope}{
Numeric giving the slope of the sigmoid used to convert mutual information values to scores ranging from zero to one.
}
  \item{shift}{
Numeric giving the relative shift of the sigmoid used to convert mutual information values to scores ranging from zero to one.
}
  \item{threshold}{
Numeric specifying the score threshold at which to consider positions for pairing.  Only applicable if \code{type} is \code{"states"} or \code{"pairs"}.
}
  \item{pseudoknots}{
Integer indicating the maximum order of pseudoknots that are acceptable.  A value of \code{0} will prevent pseudoknots in the structure, whereas \code{1} (the default) will search for first-order psuedoknots.
}
  \item{weight}{
A numeric vector of weights for each sequence, or a single number implying equal weights.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{PredictDBN} employs an extension of the method described by Freyhult et al. (2005) for determining a consensus RNA secondary structure.  It uses the mutual information (\eqn{H}) measure to find covarying positions in a multiple sequence alignment.  The original method is modified by the addition of different weights for each type of base pairing and each input sequence.  The formula for mutual information between positions \eqn{i} and \eqn{j} then becomes:

\deqn{H(i,j) = \sum_{XY \in bp}^{} \left( impact(XY) \cdot f_{i,j}(XY) \cdot \log_2 \left( \frac{f_{i,j}(XY)}{f_{i}(X) \cdot f_{j}(Y)} \right) \right)}

where, \eqn{bp} denotes the base pairings A/U, C/G, and G/U; \code{impact} is their weight; \eqn{f} is the frequency of single bases or pairs weighted by the corresponding \code{weight} of each sequence.

A penalty is then added for bases that are inconsistent with pairing:

\deqn{H_{mod}(i,j) = H(i,j) + \sum_{XY \notin bp}^{} \Big( impact(XY) \cdot f_{i,j}(XY) \Big)}

Next an average product correction (Buslje et al., 2009) is applied to the matrix \eqn{H}:

\deqn{H_{APC}(i,j) = H_{mod}(i,j) - avgProdCorr \cdot \frac{\overline{H_{mod}(i,.)} \cdot \overline{H_{mod}(.,j)}}{\overline{H_{mod}(.,.)}}}

The mutual information values are then rescaled between \code{0} and \code{1} by applying a sigmoidal transformation, which is controlled by \code{shift} and \code{slope}:

\deqn{H_{final}(i,j) = \left( 1 + \exp \left( slope \cdot log_e \left( \frac{H_{APC}(i,j)}{shift \cdot H_{APC}[n]} \right) \right) \right)^{-1}}

where, \eqn{n} is the number of positions having \code{minOccupancy} divided by two (i.e., the maximum possible number of paired positions) and \eqn{H_{APC}[n]} denotes the \eqn{n^{th}} highest value in the matrix \eqn{H_{APC}}.

If \code{type} is \code{"states"} or \code{"pairs"}, the secondary structure is determined using a variant of the Nussinov algorithm similar to that described by Venkatachalam et al. (2014).  Pairings with a score below \code{threshold} are not considered during the traceback.  If \code{psuedoknots} is greater than \code{0}, paired positions are removed from consideration and the method is applied again to find pseudoknots.

In practice the secondary structure prediction is most accurate when the input alignment is of high quality, contains a wide diversity of sequences, the number of sequences is large, no regions are completely conserved across all sequences, and most of the sequences span the entire alignment (i.e., there are few partial/incomplete sequences).
}
\value{
If \code{type} is \code{"states"} (the default), then the output is a character vector with the predicted secondary structure assignment for each position in \code{myXStringSet}.  Standard dot-bracket notation (DBN) is used, where ``.'' signifies an unpaired position, ``('' and ``)'' a paired position, and successive ``[]'', ``\{\}'', and ``<>'' indicate increasing order pseudoknots.

If \code{type} is \code{"pairs"}, then a matrix is returned with one row for each base pairing and three columns giving the positions of the paired bases and their \code{pseudoknot} order.

If \code{type} is \code{"scores"}, then a matrix of three rows is returned, where the values in a column represent the maximum score for a state in each position.  Columns sum to \code{1} if the position was above \code{minOccupancy} and \code{0} otherwise.

If \code{type} is \code{"structures"}, then the output is a list with one element for each sequence in \code{myXStringSet}.  Each list element contains a matrix of dimension 3 (each state) by the number of nucleotides in the sequence.  Columns of the matrix sum to zero where the nucleotide was located in a position that was below \code{minOccupancy}.  Otherwise, positions are considered paired if they are consistent with pairing (i.e., A/U, C/G, or G/U) in the consensus secondary structure.
}
\references{
Buslje, C., et al. (2009). Correction for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information. Bioinformatics, \bold{25(9)}, 1125-1131.

Freyhult, E., et al. (2005). Predicting RNA Structure Using Mutual Information. Applied Bioinformatics, \bold{4(1)}, 53-59.

Venkatachalam, B., et al. (2014). Faster algorithms for RNA-folding using the Four-Russians method. Algorithms for Molecular Biology : AMB, \bold{9(1)}, 1-12.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{PredictHEC}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
rna <- SearchDB(db, type="RNAStringSet")
p <- PredictDBN(rna, "states")
p

# color paired bases in the sequences
w <- which(strsplit(p, "")[[1]] != ".")
BrowseSeqs(c(BStringSet(p), BStringSet(rna)),
           colorPatterns=rep(w, each=2),
           patterns=RNA_BASES)

PredictDBN(rna, "pairs") # paired positions
}