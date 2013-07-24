\name{MODELS}
\alias{MODELS}

\title{
Available Models of DNA Evolution
}
\description{
The \code{MODELS} character vector contains the models of DNA evolution that can be used by \code{IdClusters}.
}
\usage{
MODELS
}
\details{
Six models of DNA evolution are available, with or without the discrete Gamma rates distribution.  These are described in order of increasing number of parameters as follows:

\bold{JC69 (Jukes and Cantor, 1969)}
The simplest substitution model that assumes equal base frequencies (1/4) and equal mutation rates.

\bold{K80 (Kimura, 1980)}
Assumes equal base frequencies, but distinguishes between the rate of transitions and transversions.

\bold{T92 (Tamura, 1992)}
In addition to distinguishing between transitions and transversions, a parameter is added to represent G+C content bias.

\bold{F81 (Felsenstein, 1981)}
Assumes equal mutation rates, but allows all bases to have different frequencies.

\bold{HKY85 (Hasegawa, Kishino and Yano, 1985)}
Distinguishes transitions from transversions and allows bases to have different frequencies.

\bold{TN93 (Tamura and Nei, 1993)}
Allows for unequal base frequencies and distinguishes between transversions and the two possible types of transitions (i.e., A <-> G & C <-> T).

\bold{+G (Yang, 1993)}
Specifying a model+G4 adds a single parameter to any of the above models to relax the assumption of equal rates among sites in the DNA sequence.  The single parameter specifies the shape of the Gamma Distribution.  The continuous distribution is represented with 2-10 discrete rates and their respective probabilities as determined by the Laguerre Quadrature method (Felsenstein, 2001).  For example, specifying a model+G8 would represent the continuous Gamma Distribution with eight rates and their associated probabilities.
}
\references{
Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. Journal of Molecular Evolution, \bold{17(6)}, 368-376.

Felsenstein, J. (2001). Taking Variation of Evolutionary Rates Between Sites into Account in Inferring Phylogenies. Journal of molecular evolution, \bold{53(4-5)}, 447-455.

Hasegawa, M., Kishino H., Yano T. (1985). Dating of human-ape splitting by a molecular clock of mitochondrial DNA. Journal of Molecular Evolution, \bold{22(2)}, 160-174.

Jukes, T. and Cantor C. (1969). Evolution of Protein Molecules. New York: Academic Press. pp. 21-132.

Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution, \bold{16(2)}, 111-120.

Tamura, K. (1992). Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C content biases. Molecular Biology and Evolution, \bold{9(4)}, 678-687.

Tamura, K. and Nei M. (1993). Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Molecular Biology and Evolution, \bold{10(3)}, 512-526.

Yang, Z. (1993). Maximum-likelihood estimation of phylogeny from DNA sequences when substitution rates differ over sites. Molecular Biology and Evolution, \bold{10(6)}, 1396-1401.
}
\seealso{
\code{\link{IdClusters}}
}
\examples{
MODELS
}
\keyword{data}
