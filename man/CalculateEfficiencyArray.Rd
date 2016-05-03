\name{CalculateEfficiencyArray}
\alias{CalculateEfficiencyArray}
\title{
Predict the Hybridization Efficiency of Probe/Target Sequence Pairs
}
\description{
Calculates the Gibbs free energy and hybridization efficiency of probe/target pairs at varying concentrations of the denaturant formamide.
}
\usage{
CalculateEfficiencyArray(probe,
                         target,
                         FA = 0,
                         dGini = 1.96,
                         Po = 10^-2.0021,
                         m = 0.1731,
                         temp = 42,
                         deltaGrules = NULL)
}
\arguments{
  \item{probe}{
A \code{DNAStringSet} object or character vector with pairwise-aligned probe sequences in 5' to 3' orientation.
}
  \item{target}{
A \code{DNAStringSet} object or character vector with pairwise-aligned target sequences in 5' to 3' orientation.
}
  \item{FA}{
A vector of one or more formamide concentrations (as percent v/v).
}
  \item{dGini}{
The initiation free energy.  The default is 1.96 [kcal/mol].
}
  \item{Po}{
The effective probe concentration.
}
  \item{m}{
The m-value defining the linear relationship of denaturation in the presence of formamide.
}
  \item{temp}{
Equilibrium temperature in degrees Celsius.
}
  \item{deltaGrules}{
Free energy rules for all possible base pairings in quadruplets.  If NULL, defaults to the parameters obtained using NimbleGen microarrays and a Linear Free Energy Model developed by Yilmaz \emph{et al}.
}
}
\details{
This function calculates the free energy and hybridization efficiency (HE) for a given formamide concentration ([FA]) using the linear free energy model given by:
\deqn{HE = Po*exp[-(dG_0 + m*FA)/RT]/(1+Po*exp[-(dG_0 + m*FA)/RT])}

The \code{probe} and \code{target} input sequences must be aligned in pairs, such that the first probe is aligned to the first target, second-to-second, and so on.  Ambiguity codes in the \code{IUPAC_CODE_MAP} are accepted in probe and target sequences.  Any ambiguities will default to perfect match pairings by inheriting the nucleotide in the same position on the opposite sequence whenever possible.  If the ambiguity results in a mismatch then ``T'', ``G'', ``C'', and ``A'' are substituted, in that order.  For example, if a probe nucleotide is ``S'' (``C'' or ``G'') then it will be considered a ``C'' if the target nucleotide in the same position is a ``C'', otherwise the ambiguity will be interpreted as a ``G''.

If \code{deltaGrules} is NULL then the rules defined in \code{data(deltaGrules)} will be used.  Note that \code{deltaGrules} of the same format may be customized for any application and specified as an input.
}
\value{
A \code{matrix} with the predicted Gibbs free energy (dG) and hybridization efficiency (HE) at each concentration of formamide ([FA]).
}
\references{
Yilmaz LS, Loy A, Wright ES, Wagner M, Noguera DR (2012) Modeling Formamide Denaturation of Probe-Target Hybrids for Improved Microarray Probe Design in Microbial Diagnostics. PLoS ONE 7(8): e43862. doi:10.1371/journal.pone.0043862.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{deltaGrules}}
}
\examples{
probes <- c("AAAAACGGGGAGCGGGGGGATACTG", "AAAAACTCAACCCGAGGAGCGGGGG")
targets <- c("CAACCCGGGGAGCGGGGGGATACTG", "TCGGGCTCAACCCGAGGAGCGGGGG")
result <- CalculateEfficiencyArray(probes, targets, FA=0:40)
dG0 <- result[, "dG_0"]
HE0 <- result[, "HybEff_0"]
plot(result[1, 1:40], xlab="[FA]", ylab="HE", main="Probe/Target # 1", type="l")
}