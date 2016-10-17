\name{Codec}
\alias{Codec}
\title{
Compression/Decompression of Character Vectors
}
\description{
Compresses character vectors into raw vectors, or decompresses raw vectors into character vectors using a variety of codecs.
}
\usage{
Codec(x,
      compression,
      compressRepeats = FALSE,
      processors = 1)
}
\arguments{
  \item{x}{
Either a character vector to be compressed, or a list of raw vectors to be decompressed.
}
  \item{compression}{
The type of compression algorithm to use when \code{x} is a character vector.  This should be (an unambiguous abbreviation of) one of \code{"nbit"} (for nucleotides), \code{"qbit"} (for quality scores), \code{"gzip"}, \code{"bzip2"}, or \code{"xz"}.  If \code{compression} is \code{"nbit"} or \code{"qbit"} then a second method can be provided for cases when \code{x} is incompressible.  Decompression type is determined automatically.  (See details section below.)
}
  \item{compressRepeats}{
Logical specifying whether to compress exact repeats and reverse complement repeats in a character vector input (\code{x}). Only applicable when \code{compression} is \code{"nbit"}.  Repeat compression in long DNA sequences generally increases compression by about 2\% while requiring three-fold more compression time.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
}
\details{
\code{Codec} can be used to compress/decompress character vectors using different algorithms.  The \code{"nbit"} and \code{"qbit"} methods are tailored specifically to nucleotides and quality scores, respectively.  These two methods will store the data as plain text (\code{"ASCII"} format) when it is incompressible.  In such cases, a second \code{compression} method can be given to use in lieu of plain text.  For example \code{compression = c("nbit", "gzip")} will use \code{"gzip"} compression when \code{"nbit"} compression is inappropriate.

When performing the reverse operation, decompression, the type of \code{compression} is automatically detected based on the unique signature ("magic number") added by each compression algorithm.
}
\value{
If \code{x} is a character vector to be compressed, the output is a list with one element containing a raw vector per character string.  If \code{x} is a list of raw vectors to be decompressed, then the output is a character vector with one string per list element.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\examples{
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- as.character(readDNAStringSet(fas)) # aligned sequences
object.size(dna)

# compression
system.time(x <- Codec(dna, compression="nbit"))
object.size(x)/sum(nchar(dna)) # bytes per position

system.time(g <- Codec(dna, compression="gzip"))
object.size(g)/sum(nchar(dna)) # bytes per position

# decompression
system.time(y <- Codec(x))
stopifnot(dna==y)

system.time(z <- Codec(g))
stopifnot(dna==z)
}