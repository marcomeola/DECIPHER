\name{DECIPHER-package}
\alias{DECIPHER-package}
\alias{DECIPHER}
\docType{package}
\title{
Database Enabled Code for Ideal Probe Hybridization Employing R
}
\description{
Database Enabled Code for Ideal Probe Hybridization Employing R (DECIPHER) is a software toolset that can be used for deciphering and managing biological sequences efficiently using the R statistical programming language.  The program is designed to be used with non-destructive workflows that guide the user through the process of importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.  Some functionality of the program is provided online through web tools.  DECIPHER is an ongoing project at the University of Wisconsin Madison and is freely available for download.
}
\details{
\tabular{ll}{
Package: \tab DECIPHER\cr
Type: \tab Package\cr
Depends: \tab R (>= 2.13.0), Biostrings (>= 2.31.9), RSQLite (>= 0.9), stats, parallel\cr
Imports: \tab methods, IRanges, XVector\cr
LinkingTo: \tab Biostrings, RSQLite, IRanges, XVector\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
Add2DB                   Add Data To A Database
AlignProfiles            Aligns Two Sets of Aligned Sequences
AlignSeqs                Aligns A Set of Unaligned Sequences
AlignTranslation         Aligns Sequences By Their Amino Acid Translation
Array2Matrix             Creates a Matrix Representation of a Microarray
BrowseDB                 View A Database Table In A Web Browser
BrowseSequences          View Sequences In A Web Browser
CalculateEfficiencyArray Predicts the Hybridization Efficiency of
                         Probe/Target Sequence Pairs
CalculateEfficiencyFISH  Predicts Thermodynamic Parameters of Probe/Target
                         Sequence Pairs
CalculateEfficiencyPCR   Predicts Amplification Efficiency of Primer Sequences
ConsensusSequence        Create A Consensus Sequence
CreateChimeras           Creates Artificial Chimeras
DB2FASTA                 Export Database Sequences to FASTA File
deltaGrules              Free Energy of Hybridization of Probe/Target
                         Quadruplets
DesignArray              Designs a set of DNA Microarray Probes for Detecting
                         Sequences
DesignPrimers            Designs Primers Targeting a Specific Group of
                         Sequences
DesignProbes             Designs FISH Probes Targeting a Specific Group of
                         Sequences
DistanceMatrix           Calculate the Distance Between Sequences
FindChimeras             Find Chimeras In A Sequence Database
FormGroups               Forms Groups By Rank
IdClusters               Cluster Sequences By Distance or Sequence
IdConsensus              Create Consensus Sequences by Groups
IdentifyByRank           Identify By Taxonomic Rank
IdLengths                Determine the Number of Bases, Nonbases, and Width
                         of Each Sequence
MaskAlignment            Masks Highly Variable Regions of An Alignment
MODELS                   Available Models of DNA Evolution
NNLS                     Sequential Coordinate-wise Algorithm for the
                         Non-negative Least Squares Problem
RESTRICTION_ENZYMES      Common Restriction Sites Named By Restriction Enzyme
SearchDB                 Obtain Specific Sequences from A Database
Seqs2DB                  Add Sequences from Text File to Database
TerminalChar             Determine the Number of Terminal Characters
TileSeqs                 Form a Set of Tiles for Each Group of Sequences
}
}
\author{
Erik Wright

Maintainer: Erik Wright <DECIPHER@cae.wisc.edu>
}
\keyword{ package }
