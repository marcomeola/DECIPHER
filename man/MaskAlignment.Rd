\name{MaskAlignment}
\alias{MaskAlignment}
\title{
Mask Highly Variable Regions of An Alignment
}
\description{
Automatically masks poorly aligned regions of an alignment based on sequence conservation and gap frequency.
}
\usage{
MaskAlignment(myXStringSet,
              windowSize = 5,
              threshold = 1,
              maxFractionGaps = 0.2,
              correction = FALSE,
              showPlot = FALSE)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of aligned sequences.
}
  \item{windowSize}{
Integer value specifying the size of the region to the left and right of the center-point to use in calculating the moving average.
}
  \item{threshold}{
Numeric giving the average entropy in bits from 0 to 2 below which a region is masked.
}
  \item{maxFractionGaps}{
Numeric specifying the maximum faction of gaps in an alignment column to be masked.
}
  \item{correction}{
Logical indicating whether to apply a small-sample size correction to columns with few letters (Yu et al., 2015).
}
  \item{showPlot}{
Logical specifying whether or not to show a plot of the positions that were kept or masked.
}
}
\details{
Poorly aligned regions of a multiple sequence alignment may lead to incorrect results in downstream analyses.  One method to mitigate their effects is to mask columns of the alignment that may be poorly aligned, such as highly-variable regions or regions with many insertions and deletions (gaps).

Highly variable regions are detected by their signature of having a low information content.  Here, information content is defined by the relative entropy of a column in the alignment (Yu et al., 2015), which is higher for conserved columns.  The relative entropy is based on the background distribution of letter-frequencies in the alignment.

A moving average of \code{windowSize} nucleotides to the left and right of the center-point is applied to smooth noise in the information content signal along the sequence.  Regions dropping below \code{threshold} bits or more than \code{maxFractionGaps} are masked in the returned alignment.
}
\value{
A \code{MultipleAlignment} object of the input type with masked columns where the input criteria are met.
}
\references{
Yu, Y.-K., et al. (2015). Log-odds sequence logos. Bioinformatics, 31(3), 324-331. http://doi.org/10.1093/bioinformatics/btu634
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignSeqs}}, \code{\link{IdClusters}}
}
\examples{
fas <- system.file("extdata", "Streptomyces_ITS_aligned.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
masked_dna <- MaskAlignment(dna, showPlot=TRUE)

# display only unmasked nucleotides for use in downstream analyses
not_masked <- as(masked_dna, "DNAStringSet")
BrowseSeqs(not_masked)

# display only masked nucleotides that are covered by the mask
masked <- masked_dna
colmask(masked, append="replace", invert=TRUE) <- colmask(masked)
masked <- as(masked, "DNAStringSet")
BrowseSeqs(masked)

# display the complete DNA sequence set including the mask
masks <- lapply(width(colmask(masked_dna)), rep, x="+")
masks <- unlist(lapply(masks, paste, collapse=""))
masked_dna <- replaceAt(dna, at=IRanges(colmask(masked_dna)), value=masks)
BrowseSeqs(masked_dna)
}