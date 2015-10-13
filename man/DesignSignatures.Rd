\name{DesignSignatures}
\alias{DesignSignatures}
\title{
Design PCR Primers for Amplifying Group-Specific Signatures
}
\description{
Aids the design of pairs of primers for amplifying a unique ``signature'' from each group of sequences.  Signatures are distinct PCR products that can be differentiated by their length, melt temperature, or sequence.
}
\usage{
DesignSignatures(dbFile,
                 tblName = "DNA",
                 identifier = "",
                 focusID = NA,
                 type = "melt",
                 resolution = 0.5,
                 levels = 10,
                 enzymes = NULL,
                 minLength = 17,
                 maxLength = 26,
                 maxPermutations = 4,
                 annealingTemp = 64,
                 P = 4e-07,
                 monovalent = 0.07,
                 divalent = 0.003,
                 dNTPs = 8e-04,
                 minEfficiency = 0.8,
                 ampEfficiency = 0.5,
                 numPrimerSets = 100,
                 minProductSize = 70,
                 maxProductSize = 400,
                 kmerSize = 8,
                 searchPrimers = 500,
                 maxDictionary = 20000,
                 primerDimer = 1e-07,
                 taqEfficiency = TRUE,
                 processors = NULL,
                 verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table where the DNA sequences are located.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  Determines the target group(s) for which primers will be designed.  If "" then all identifiers are selected.
}
  \item{focusID}{
Optional character string specifying which of the \code{identifier}s will be used in the initial step of designing primers.  If \code{NA} (the default), then the \code{identifier} with the most sequence information is used as the \code{focusID}.
}
  \item{type}{
Character string indicating the type of signature being used to differentiate the PCR products from each group.  This should be (an abbreviation of) one of \code{"melt"}, \code{"length"}, or \code{"sequence"}.
}
  \item{resolution}{
Numeric specifying the ``resolution'' of the experiment, or a vector giving the boundaries of bins.  (See details section below.)
}
  \item{levels}{
Numeric giving the number of ``levels'' that can be distinguished in each bin.  (See details section below.)
}
  \item{enzymes}{
Named character vector providing the cut sites of one or more restriction enzymes.  Cut sites must be delineated in the same format as \code{\link{RESTRICTION_ENZYMES}}.
}
  \item{minLength}{
Integer providing the minimum length of primers to consider in the design.
}
  \item{maxLength}{
Integer providing the maximum length of primers to consider in the design.
}
  \item{maxPermutations}{
Integer providing the maximum number of permutations allowed in a forward or reverse primer to attain greater coverage of sequences.
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
Numeric giving the minimum efficiency of hybridization desired for the primer set.
}
  \item{ampEfficiency}{
Numeric giving the minimum efficiency required for theoretical amplification of the primers.  Note that \code{ampEfficiency} must be less than or equal to \code{minEfficiency}.  Lower values of \code{ampEfficiency} will allow for more PCR products, although very low values are unrealistic experimentally.
}
  \item{numPrimerSets}{
Integer giving the optimal number of primer sets (forward and reverse primer sets) to design.
}
  \item{minProductSize}{
Integer giving the minimum number of nucleotides desired in the PCR product.
}
  \item{maxProductSize}{
Integer giving the maximum number of nucleotides desired in the PCR product.
}
  \item{kmerSize}{
Integer giving the size of k-mers to use in the preliminary search for potential primers.
}
  \item{searchPrimers}{
Numeric specifying the number of forward and reverse primers to use in searching for potential PCR products.  A lower value will result in a faster search, but potentially neglect some useful primers.
}
  \item{maxDictionary}{
Numeric giving the maximum number of primers to search for simultaneously in any given step.
}
  \item{primerDimer}{
Numeric giving the maximum amplification efficiency of potential primer-dimer products.
}
  \item{taqEfficiency}{
Logical determining whether to make use of elongation efficiency to increase predictive accuracy for \emph{Taq} DNA Polymerase amplifying primers with mismatches near the 3' terminus.  Note that this should be set to FALSE if using a high-fidelity polymerase with 3' to 5' exonuclease activity.
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Signatures are group-specific PCR products that can be differentiated by either their melt temperature profile, length, or sequence.  \code{DesignSignatures} assists in finding the optimal pair of forward and reverse primers for obtaining a distinguishable signature from each group of sequences.  Groups are delineated by their unique \code{identifier} in the database.  The algorithm works by progressively narrowing the search for optimal primers:  (1) the most frequent k-mers are found;  (2) these are used to design primers initially matching the \code{focusID} group;  (3) the most common forward and reverse primers are selected based on all of the groups, and ambiguity is added up to \code{maxPermutations};  (4) a final search is performed to find the optimal forward and reverse primer.  Pairs of primers are scored by the distance between the signatures generated for each group, which depends on the \code{type} of experiment.

The arguments \code{resolution} and \code{levels} control the theoretical resolving power of the experiment.  The \code{signature} for a group is discretized or grouped into ``bins'' each with a certain magnitude of the signal.  Here \code{resolution} determines the separation between distinguishable ``bins'', and \code{levels} controls the range of values in each bin.  A high-accuracy experiment would have many bins and/or many levels.  While \code{levels} is interpreted similarly for every \code{type} of experiment, \code{resolution} is treated differently depending on \code{type}.  If \code{type} is \code{"melt"}, then \code{resolution} can be either a vector of different melt temperatures, or a single number giving the change in temperatures that can be differentiated.  A high-resolution melt (HRM) assay would typically have a resolution between 0.25 and 1 degree Celsius.  If \code{type} is \code{"length"} then resolution is either the number of bins between the \code{minProductSize} and \code{maxProductSize}, or the bin boundaries.  For example, \code{resolution} can be lower (wider bins) at long lengths, and higher (narrower bins) at shorter lengths.  If \code{type} is \code{"sequence"} then \code{resolution} sets the k-mer size used in differentiating amplicons.  Oftentimes, 4 to 6-mers are used for the classification of amplicons.

The signatures can be diversified by using a restriction enzyme to digest the PCR products when \code{type} is \code{"melt"} or \code{"length"}.  If \code{enzymes} are supplied then the an additional search is made to find the best enzyme to use with each pair of primers.  In this case, the output includes all of the primer pairs, as well as any \code{enzymes} that will digest the PCR products of that primer pair.  The output is re-scored to rank the top primer pair and enzyme combination.  Note that \code{enzymes} is inapplicable when \code{type} is \code{"sequence"} because restriction enzymes do not alter the sequence of the DNA.  Also, it is recommended that only a subset of the available \code{RESTRICTION_ENZYMES} are used as input \code{enzymes} in order to accelerate the search for the best enzyme.
}
\value{
A \code{data.frame} with the top-scoring pairs of forward and reverse primers, their score, the total number of PCR products, and associated columns for the restriction enzyme (if \code{enzyme} is not \code{NULL}).
}
\references{
Coming soon!
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AmplifyDNA}}, \code{\link{CalculateEfficiencyPCR}}, \code{\link{DesignPrimers}}, \code{\link{DigestDNA}}, \code{\link{Disambiguate}}, \code{\link{MeltDNA}}, \code{\link{RESTRICTION_ENZYMES}}
}
\examples{
# below are suggested inputs for different types of experiments
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")

\dontrun{
# High Resolution Melt (HRM) assay:
primers <- DesignSignatures(db,
           resolution=seq(80, 100, 0.25), # degrees Celsius
           minProductSize=55, # base pairs
           maxProductSize=400)

# Primers for next-generation sequencing:
primers <- DesignSignatures(db,
           type="sequence",
           minProductSize=300, # base pairs
           maxProductSize=700,
           resolution=5, # 5-mers
           levels=5)

# Primers for community fingerprinting:
primers <- DesignSignatures(db,
           type="length",
           levels=2, # presence/absence
           minProductSize=200, # base pairs
           maxProductSize=1400,
           resolution=c(seq(200, 700, 3),
                        seq(705, 1000, 5),
                        seq(1010, 1400, 10)))

# Primers for restriction fragment length polymorphism (RFLP):
data(RESTRICTION_ENZYMES)
myEnzymes <- RESTRICTION_ENZYMES[c("EcoRI", "HinfI", "SalI")]
primers <- DesignSignatures(db,
           type="length",
           levels=2, # presence/absence
           minProductSize=200, # base pairs
           maxProductSize=600,
           resolution=c(seq(50, 100, 3),
                        seq(105, 200, 5),
                        seq(210, 600, 10)),
           enzymes=myEnzymes)
}
}