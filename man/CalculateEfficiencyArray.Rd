\name{CalculateEfficiencyArray}
\alias{CalculateEfficiencyArray}
\title{
Predicts the Hybridization Efficiency of Probe/Target Sequence Pairs
}
\description{
Calculates the Gibb's free energy and hybridization efficiency of probe/target pairs at varying concentrations of the denaturant formamide.
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

Probe and target input sequences must be entered in pairwise alignment, such as that given by \code{pairwiseAlignment}.  Only "A", "C", "G", "T", and "-" characters are permitted in the probe sequence.

If deltaGrules is NULL then the rules defined in \code{data(deltaGrules)} are used.
}
\value{
A \code{matrix} with the predicted Gibb's free energy (dG) and hybridization efficiency (HE) at each concentration of formamide ([FA]).
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