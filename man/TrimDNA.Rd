\name{TrimDNA}
\alias{TrimDNA}
\title{
Trims DNA Sequences to the High Quality Region Between Patterns
}
\description{
Aids in trimming DNA sequences to the high quality region between a set of patterns that are potentially present on the left and right sides.
}
\usage{
TrimDNA(myDNAStringSet,
        leftPatterns,
        rightPatterns,
        type = "ranges",
        quality = NULL,
        maxDistance = 0.1,
        minOverlap = 5,
        allowInternal = TRUE,
        alpha = 0.2,
        threshold = 0.01,
        minWidth = 36,
        verbose = TRUE)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object containing the sequences to be trimmed.
}
  \item{leftPatterns}{
A \code{DNAStringSet} or character vector of patterns to remove from the left side of \code{myDNAStringSet}, or \code{""} to prevent trimming patterns on the left.
}
  \item{rightPatterns}{
A \code{DNAStringSet} or character vector of patterns to remove from the right side of \code{myDNAStringSet}, or \code{""} to prevent trimming patterns on the right.
}
  \item{type}{
Character string indicating the type of results desired.  This should be (an abbreviation of) either \code{"ranges"}, \code{"sequences"} or \code{"both"}.
}
  \item{quality}{
Either \code{NULL} (the default) to skip quality trimming, or a \code{PhredQuality}, \code{SolexaQuality}, or \code{IlluminaQuality} object containing the quality scores corresponding to \code{myDNAStringSet}.
}
  \item{maxDistance}{
Numeric between zero and one giving the maximum distance of a match from the \code{leftPatterns} and \code{rightPatterns} to initiate trimming. For example, \code{0.1} (the default) would allow up to 10\% mismatches between a pattern and sequence.
}
  \item{minOverlap}{
Integer specifying the minimum number of nucleotides the \code{leftPatterns} and \code{rightPatterns} must overlap a sequence to initiate trimming.
}
  \item{allowInternal}{
Logical initiating whether to search for the \code{leftPatterns} and \code{rightPatterns} within \code{myDNAStringSet}, or (\code{FALSE} for) only overlapping the ends.
}
  \item{alpha}{
Numeric between zero and one giving the smoothing parameter for an exponential moving average that is applied to the quality scores before trimming.  Higher values result in less smoothing than lower values.
}
  \item{threshold}{
Numeric between zero and one specifying the threshold above which to trim the poor quality regions of the sequence.  Higher values allow more sequence to be preserved at the expense of a greater error rate.
}
  \item{minWidth}{
Integer giving the minimum number of nucleotides a pattern must overlap the sequence to initiate trimming.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
After a sequencing run, it is often necessary to trim the resulting sequences to the high quality region located between a set of patterns.  \code{TrimDNA} works as follows:  first left and right patterns are identified within the sequences if \code{allowInternal} is \code{TRUE} (the default).  If the patterns are not found internally, then a search is conducted at the flanking ends for patterns that partially overlap the sequence.  The region between the \code{leftPatterns} and \code{rightPatterns} is then returned, unless quality information is provided.

If \code{quality} contains quality scores, these are converted to error probabilities and an exponential moving average is applied to smooth the signal.  The longest region between the \code{leftPatterns} and \code{rightPatterns} where the average error probability is below \code{threshold} is then returned.
}
\value{
\code{TrimDNA} can return two \code{type}s of results: \code{IRanges} that can be used for trimming \code{myDNAStringSet}, or a trimmed \code{DNAStringSet} containing only those sequences over \code{minWidth} nucleotides after trimming.  Note that ambiguity codes (\code{IUPAC_CODE_MAP}) are supported in the \code{leftPatterns} and \code{rightPatterns}, but not in \code{myDNAStringSet} to prevent trivial matches (e.g., runs of N's).

If \code{type} is \code{"ranges"} (the default) the output is an \code{IRanges} object with the start, end, and width of every sequence.  This information can be accessed with the corresponding accessor function (see examples below).  Note that the start will be \code{1} and the end will be \code{0} for sequences that were not at least \code{minWidth} nucleotides after trimming.

If \code{type} is \code{"sequences"} then the trimmed sequences are returned that are at least \code{minWidth} nucleotides in length.

If \code{type} is \code{"both"} the output is a list of two components, the first containing the ranges and the second containing the sequences.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{CorrectFrameshifts}}
}
\examples{
dna <- DNAStringSet("AAAAAAAAAATTACTTCCCCCCCCCC")
qscores <- PhredQuality("0000000000AAAAAAAAAAAAAAAA")

x <- TrimDNA(dna,
            leftPatterns="AAAAAA",
            rightPatterns="CCCCCC",
            quality=qscores,
            minWidth=1,
            allowInternal=TRUE,
            type="both")

x[[1]]
start(x[[1]])
end(x[[1]])
width(x[[1]])

x[[2]]
}
