\name{DECIPHER-package}
\alias{DECIPHER-package}
\alias{DECIPHER}
\docType{package}
\title{
Tools for curating, analyzing, and manipulating biological sequences
}
\description{
DECIPHER is a software toolset that can be used for deciphering and managing biological sequences efficiently using the R statistical programming language.  The program is designed to be used with non-destructive workflows for importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.
}
\details{
\tabular{ll}{
Package: \tab DECIPHER\cr
Type: \tab Package\cr
Depends: \tab R (>= 2.13.0), Biostrings (>= 2.35.12), RSQLite (>= 1.0), stats, parallel\cr
Imports: \tab methods, DBI, S4Vectors, IRanges, XVector\cr
LinkingTo: \tab Biostrings, RSQLite, S4Vectors, IRanges, XVector\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
Add2DB                   Add Data to a Database
AdjustAlignment          Improve An Existing Alignment By Adjusting Gap
                         Placements
AlignDB                  Align Two Sets of Aligned Sequences in a Sequence
                         Database
AlignProfiles            Align Two Sets of Aligned Sequences
AlignSeqs                Align a Set of Unaligned Sequences
AlignSynteny             Pairwise Aligns Syntenic Blocks
AlignTranslation         Align Sequences By Their Amino Acid Translation
AmplifyDNA               Simulate Amplification of DNA by PCR
Array2Matrix             Create a Matrix Representation of a Microarray
BrowseDB                 View a Database Table in a Web Browser
BrowseSeqs               View Sequences in a Web Browser
CalculateEfficiencyArray Predict the Hybridization Efficiency of
                         Probe/Target Sequence Pairs
CalculateEfficiencyFISH  Predict Thermodynamic Parameters of Probe/Target
                         Sequence Pairs
CalculateEfficiencyPCR   Predict Amplification Efficiency of Primer Sequences
Codec                    Compression/Decompression of Character Vectors
ConsensusSequence        Create a Consensus Sequence
CorrectFrameshifts       Corrects Frameshift Errors In Protein Coding
                         Sequences
CreateChimeras           Create Artificial Chimeras
DB2Seqs                  Export Database Sequences to a FASTA or FASTQ File
deltaGrules              Free Energy of Hybridization of Probe/Target
                         Quadruplets
deltaHrules              Change in Enthalpy of Hybridization of Primer/Target
                         Quadruplets in Solution
deltaSrules              Change in Entropy of Hybridization of Primer/Target
                         Quadruplets in Solution
DesignArray              Design a set of DNA Microarray Probes for Detecting
                         Sequences
DesignPrimers            Design Primers Targeting a Specific Group of
                         Sequences
DesignProbes             Design FISH Probes Targeting a Specific Group of
                         Sequences
DesignSignatures         Design PCR Primers for Amplifying Group-Specific
                         Signatures
DigestDNA                Simulate Restriction Digestion of DNA
Disambiguate             Expand Ambiguities into All Permutations of a
                         DNAStringSet
DistanceMatrix           Calculate the Distance Between Sequences
FindChimeras             Find Chimeras in a Sequence Database
FindSynteny              Finds Synteny in a Sequence Database
FormGroups               Forms Groups By Rank
HEC_MI                   Mutual Information for Protein Secondary Structure
                         Prediction
IdClusters               Cluster Sequences By Distance or Sequence
IdConsensus              Create Consensus Sequences by Groups
IdentifyByRank           Identify By Taxonomic Rank
IdLengths                Determine the Number of Bases, Nonbases, and Width
                         of Each Sequence
MaskAlignment            Mask Highly Variable Regions of An Alignment
MeltDNA                  Simulate Melting of DNA
MIQS                     MIQS Amino Acid Substitution Matrix
MODELS                   Available Models of DNA Evolution
NNLS                     Sequential Coordinate-wise Algorithm for the
                         Non-negative Least Squares Problem
OrientNucleotides        Orient nucleotide sequences
PredictHEC               Predict Protein Secondary Structure
Read Dendrogram          Read a Dendrogram from a Newick Formatted File
RESTRICTION_ENZYMES      Common Restriction Enzyme's Cut Sites
SearchDB                 Obtain Specific Sequences from a Database
Seqs2DB                  Add Sequences from Text File to Database
StaggerAlignment         Produce a Staggered Alignment
Synteny-class            Synteny blocks and hits
TerminalChar             Determine the Number of Terminal Characters
TileSeqs                 Form a Set of Tiles for Each Group of Sequences
WriteDendrogram          Write a Dendrogram to Newick Format
}
}
\author{
Erik Wright

Maintainer: Erik Wright <DECIPHER@cae.wisc.edu>
}
\keyword{ package }
