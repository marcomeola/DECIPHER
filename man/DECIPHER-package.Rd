\name{DECIPHER-package}
\alias{DECIPHER-package}
\alias{DECIPHER}
\docType{package}
\title{
Database Enabled Code for Ideal Probe Hybridization Employing R
}
\description{
Database Enabled Code for Ideal Probe Hybridization Employing R (DECIPHER) is a software toolset that can be used for deciphering and managing DNA sequences efficiently using the R statistical programming language.  The program is designed to be used with non-destructive workflows that guide the user through the process of importing, maintaining, analyzing, manipulating, and exporting a massive amount of DNA sequences.  Some functionality of the program is provided online through web tools.  DECIPHER is an ongoing project in the Environmental Engineering Department at the University of Wisconsin Madison and is freely available for download.
}
\details{
\tabular{ll}{
Package: \tab DECIPHER\cr
Type: \tab Package\cr
Version: \tab 1.3.5\cr
Date: \tab 2011-07-22\cr
Depends: \tab R (>= 2.13.0), Biostrings (>= 2.16), RSQLite (>= 0.9), IRanges, stats\cr
Imports: \tab Biostrings, RSQLite, IRanges, stats\cr
LinkingTo: \tab Biostrings, RSQLite, IRanges, stats\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
Add2DB                  Add Data To A Database
BrowseDB                View A Database Table In A Web Browser
BrowseSequences         View Sequences In A Web Browser
ConsensusSequence       Create A Consensus Sequence
CreateChimeras          Creates Artificial Chimeras
DB2FASTA                Export Database to FASTA File
DistanceMatrix          Calculate the Distance Between DNA Sequences
FindChimeras            Find Chimeras In A Sequence Database
FormGroups              Forms Groups By Rank
IdClusters              Cluster Sequences By Distance
IdConsensus             Create Consensus Sequences by Groups
IdLengths               Determine the Number of Bases and Nonbases In
                        Each Sequence
IdentifyByRank          Update Identifier To Level of Taxonomic Rank
SearchDB                Obtain Specific Sequences from A Database
Seqs2DB                 Add Sequences from Text File to Database
TerminalChar            Determine the Number of Terminal Gaps
}
}
\author{
Erik Wright

Maintainer: Erik Wright <DECIPHER@cae.wisc.edu>
}
\keyword{ package }
