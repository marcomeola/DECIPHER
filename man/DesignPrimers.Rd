\name{DesignPrimers}
\alias{DesignPrimers}
\title{
Designs Primers Targeting a Specific Group of Sequences
}
\description{
Assists in the design of primer sets targeting a specific group of sequences while minimizing the potential to cross-amplify other groups of sequences.
}
\usage{
DesignPrimers(tiles,
              identifier = "",
              start = 1,
              end = NULL,
              minLength = 17,
              maxLength = 26,
              maxPermutations = 4,
              minCoverage = 0.9,
              minGroupCoverage = 0.2,
              annealingTemp = 64,
              P = 4e-07,
              monovalent = 0.07,
              divalent = 0.003,
              dNTPs = 8e-04,
              minEfficiency = 0.8,
              worstScore = -Inf,
              numPrimerSets = 0,
              minProductSize = 75,
              maxProductSize = 1200,
              maxSearchSize = 1500,
              batchSize = 1000,
              maxDistance = 0.4,
              primerDimer = 1e-07,
              ragged5Prime = TRUE,
              taqEfficiency = TRUE,
              induceMismatch = FALSE,
              verbose = TRUE)
}
\arguments{
  \item{tiles}{
A set of tiles representing each group of sequences, as in the format created by the function \code{TileSeqs}.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  Determines the target group(s) for which primers will be designed.  If "" then all identifiers are selected.
}
  \item{start}{
Integer specifying the starting position in the alignment where potential forward primer target sites begin.  Preferably a position that is included in most sequences in the alignment.
}
  \item{end}{
Integer specifying the ending position in the alignment where potential reverse primer target sites end.  Preferably a position that is included in most sequences in the alignment.
}
  \item{minLength}{
Integer providing the minimum length of primers to consider in the design.
}
  \item{maxLength}{
Integer providing the maximum length of primers to consider in the design, which must be less than or equal to the \code{maxLength} of tiles.
}
  \item{maxPermutations}{
Integer providing the maximum number of permutations considered as part of a forward or reverse primer set.
}
  \item{minCoverage}{
Numeric giving the minimum fraction of the target group's sequences that must be covered with the primer set.
}
  \item{minGroupCoverage}{
Numeric giving the minimum fraction of the target group that must have sequence information (not terminal gaps) in the region covered by the primer set.
}
  \item{annealingTemp}{
Numeric indicating the desired annealing temperature that will be used in the PCR experiment.
}
  \item{P}{
Numeric giving the molar concentration of primers in the reaction.
}
  \item{monovalent}{
The molar concentration of monovalent ([Na] and [K]) ions in solution that will be used to determine a sodium equivalent concentration.
}
  \item{divalent}{
The molar concentration of divalent ([Mg]) ions in solution that will be used to determine a sodium equivalent concentration.
}
  \item{dNTPs}{
Numeric giving the molar concentration of free nucleotides added to the solution that will be used to determine a sodium equivalent concentration.
}
  \item{minEfficiency}{
Numeric giving the minimum efficiency of hybridization desired for the primer set.  Note that an efficiency of 99% (0.99) will greatly lower predicted specificity of the primer set, however an efficiency of 50% (0.5) may be too low in actuality to amplify the target group due to error in melt temperature predictions.
}
  \item{worstScore}{
Numeric specifying the score cutoff to remove target sites from consideration.  For example, a \code{worstScore} of -5 will remove all primer sets scoring below -5, although this may eventually result in no primer sets meeting this criteria.
}
  \item{numPrimerSets}{
Integer giving the optimal number of primer sets (forward and reverse primer sets) to design.  If set to zero then all possible forward and reverse primers are returned, but the primer sets minimizing potential false positives are not chosen.
}
  \item{minProductSize}{
Integer giving the minimum number of nucleotides desired in the PCR product.
}
  \item{maxProductSize}{
Integer giving the maximum number of nucleotides desired in the PCR product.
}
  \item{maxSearchSize}{
Integer giving the maximum number of nucleotides to search for false priming upstream and downstream of the expected binding site.
}
  \item{batchSize}{
Integer specifying the number of primers to simulate hybridization per batch.  See the Description section below.
}
  \item{maxDistance}{
Numeric specifying the maximal fraction of mismatched base pairings on a rolling basis beginning from the 3' end of the primer.
}
  \item{primerDimer}{
Numeric giving the maximum amplification efficiency of primer-dimer products.
}
  \item{ragged5Prime}{
Logical specifying whether the 5' end or 3' end of primer permutations targeting the same site should be varying lengths.
}
  \item{taqEfficiency}{
Logical determining whether to make use of elongation efficiency and maxDistance to increase predictive accuracy for \emph{Taq} DNA Polymerase amplifying primers with mismatches near the 3' terminus.  Note that this should be set to FALSE if using a high-fidelity polymerase with 3' to 5' exonuclease activity.
}
  \item{induceMismatch}{
Logical or integer specifying whether to induce a mismatch in the primer with the template DNA.  If \code{TRUE} then a mismatch is induced at the 6th primer position.  If an integer value is provided between 2 and 6 then a mismatch is induced in that primer position.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Primers are designed for use with \emph{Taq} DNA Polymerase to maximize sensitivity and specificity for the target group of sequences.  The design makes use of \emph{Taq}'s bias against certain 3' terminal mismatch types in order to increase specificity further than can be achieve with hybridization efficiency alone.

Primers are designed from a set of tiles to target each \code{identifier} while minimizing affinity for all other tiled groups.  Arguments provide constraints that ensure the designed primer sets meet the specified criteria as well as being optimized for the particular experimental conditions.  A search is conducted through all tiles in the same alignment position to estimate the chance of cross-amplification with a non-target group.

If \code{numPrimers} is greater than or equal to one then the set of forward and reverse primers that minimizes potential false positive overlap is returned.  This will also initiate a thorough search through all target sites upstream and downstream of the expected binding sites to ensure that the primers do not bind to nearby positions.  Lowering the \code{maxSearchSize} will speed up the thorough search at the expense of potentially missing an unexpected target site.  The number of possible primer sets assessed is increased with the size of \code{numPrimers}.
}
\value{
A different \code{data.frame} will be returned depending on number of primer sets requested.  If no primer sets are required then columns contain the forward and reverse primers for every possible position scored by their potential to amplify other identified groups.  If one or more primer sets are requested then columns contain information for the optimal set of forward and reverse primers that should be used in combination to give the fewest potential false positives.
}
\references{
Coming Soon!
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\note{
The program OligoArrayAux (\url{http://mfold.rna.albany.edu/?q=DINAMelt/OligoArrayAux}) must be installed in a location accessible by the system.  For example, the following code should print the installed OligoArrayAux version when executed from the R console:

\code{system("hybrid-min -V")}
}
\seealso{
\code{\link{CalculateEfficiencyPCR}}, \code{\link{TileSeqs}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
# not run (must have OligoArrayAux installed first):
#tiles <- TileSeqs(db, identifier=c("Acinetobacter","Pseudomonas"))
#primers <- DesignPrimers(tiles, identifier="Acinetobacter", start=280, end=420,
#	minProductSize=50, numPrimerSets=1)
}