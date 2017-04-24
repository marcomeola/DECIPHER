\name{DesignProbes}
\alias{DesignProbes}
\title{
Design FISH Probes Targeting a Specific Group of Sequences
}
\description{
Assists in the design of single or dual probes targeting a specific group of sequences while minimizing the potential to cross-hybridize with other groups of sequences.
}
\usage{
DesignProbes(tiles,
             identifier = "",
             start = 1,
             end = NULL,
             minLength = 17,
             maxLength = 26,
             maxPermutations = 4,
             minCoverage = 0.9,
             minGroupCoverage = 0.2,
             hybTemp = 46,
             P = 2.5e-07,
             Na = 1,
             FA = 35,
             minEfficiency = 0.5,
             worstScore = -Inf,
             numProbeSets = 0,
             batchSize = 1000,
             target = "SSU",
             verbose = TRUE)
}
\arguments{
  \item{tiles}{
A set of tiles representing each group of sequences, as in the format created by the function \code{TileSeqs}.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  Determines the target group(s) for which probes will be designed.  If "" then all identifiers are selected.
}
  \item{start}{
Integer specifying the starting position in the alignment where potential target sites begin.  Preferably a position that is included in most sequences in the alignment.
}
  \item{end}{
Integer specifying the ending position in the alignment where potential target sites end.  Preferably a position that is included in most sequences in the alignment.
}
  \item{minLength}{
Integer providing the minimum length of probes to consider in the design.
}
  \item{maxLength}{
Integer providing the maximum length of probes to consider in the design, which must be less than or equal to the \code{maxLength} of tiles.
}
  \item{maxPermutations}{
Integer providing the maximum number of probe permutations required to reach the desired coverage of a target site.
}
  \item{minCoverage}{
Numeric giving the minimum fraction of the target group's sequences that must be covered by the designed probe(s).
}
  \item{minGroupCoverage}{
Numeric giving the minimum fraction of the target group that must have sequence information (not terminal gaps) in the target site's region.
}
  \item{hybTemp}{
Numeric specifying the hybridization temperature, typically \code{46} degrees Celsius.
}
  \item{P}{
Numeric giving the molar concentration of probes during hybridization.
}
  \item{Na}{
Numeric giving the molar sodium concentration in the hybridization buffer.  Values may range between 0.01M and 1M.  Note that salt correction from \code{1} molar is not available for the thermodynamic rules of RNA/RNA interactions.
}
  \item{FA}{
Numeric concentration (as percent v/v) of the denaturant formamide in the hybridization buffer.
}
  \item{minEfficiency}{
Numeric giving the minimum equilibrium hybridization efficiency desired for designed probe(s) at the defined experimental conditions.
}
  \item{worstScore}{
Numeric specifying the score cutoff to remove target sites from consideration.  For example, a \code{worstScore} of -5 will remove all probes scoring below -5, although this may eventually result in no probes meeting the design criteria.
}
  \item{numProbeSets}{
Integer giving the optimal number of dual probe sets to design.  If set to zero then all potential single probes are returned, and the probe sets minimizing potential false cross-hybridizations are not chosen.
}
  \item{batchSize}{
Integer specifying the number of probes to simulate hybridization per batch that is passed to \code{CalculateEfficiencyFISH}.
}
  \item{target}{
The target molecule used in the generation of tiles.  Either \code{"SSU"} for the small-subunit rRNA, \code{"LSU"} for the large-subunit rRNA, or "Other".  Used to determine the domain for \code{dG3} calculations, which is plus or minus 200 nucleotides of the target site if \code{"Other"}.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Probes are designed to maximize sensitivity and specificity to the target group(s) (\code{identifier}(s)).  If \code{numProbeSets > 0} then that many pairs of probes with minimal cross-hybridization overlap are returned, enabling increased specificity with a dual-color approach.

Probes are designed from a set of tiles to target each \code{identifier} while minimizing affinity for all other tiled groups.  Arguments provide constraints that ensure the designed probes meet the specified criteria as well as being optimized for the particular experimental conditions.  A search is conducted through all tiles in the same alignment position to estimate the chance of cross-hybridization with a non-target group.

Two models are used in design, both of which were experimentally calibrated using denaturation profiles from 5 organisms belonging to all three domains of life.  Probe lengths are chosen to meet the \code{minEfficiency} using a fast model of probe-target hybridization.  Candidate probes are then confirmed using a slower model that also takes into account probe-folding and target-folding.  Finally, probes are scored for their inability to cross-hybridize with non-target groups by using the fast model and taking into account any mismatches.
}
\value{
A different \code{data.frame} will be returned depending on number of primer sets requested.  If no probe sets are required then columns contain the designed probes for every possible position scored by their potential to cross-hybridize with other identified groups.  If one or more probe sets are requested then columns contain information for the optimal set of probes (probe one and probe two) that could be used in combination to give the fewest potential cross-hybridizations.
}
\references{
ES Wright et al. (2014) "Automated Design of Probes for rRNA-Targeted Fluorescence In Situ Hybridization Reveals the Advantages of Using Dual Probes for Accurate Identification." Applied and Environmental Microbiology, doi:10.1128/AEM.01685-14.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\note{
The program OligoArrayAux (\url{http://mfold.rna.albany.edu/?q=DINAMelt/OligoArrayAux}) must be installed in a location accessible by the system.  For example, the following code should print the installed OligoArrayAux version when executed from the R console:

\code{system("hybrid-min -V")}

To install OligoArrayAux from the downloaded source folder on Unix-like platforms, open the shell (or Terminal on Mac OS) and type:

\code{cd oligoarrayaux # change directory to the correct folder name}

\code{./configure}

\code{make}

\code{sudo make install}
}
\seealso{
\code{\link{CalculateEfficiencyFISH}}, \code{\link{TileSeqs}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
# not run (must have OligoArrayAux installed first):
\dontrun{tiles <- TileSeqs(db, identifier=c("Enterobacteriales","Pseudomonadales"))}
\dontrun{probes <- DesignProbes(tiles, identifier="Enterobacteriales", start=280, end=420)}
}