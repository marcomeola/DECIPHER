Codec <- function(x,
	compression="auto",
	compressRepeats=FALSE,
	processors=1) {
	
	# error checking
	if (length(compression) != 1 || !is.character(compression))
		stop("compression must be a single character string.")
	if (!(compression %in% c("auto", "nbit", "gzip", "bzip2", "xz")))
		stop("Invalid type of compression.")
	if (typeof(x)=="list") {
		if (length(x)==0)
			return(character())
		if (!all(unlist(lapply(x, is.raw))))
			stop("All elements of x must be raw vectors.")
	} else if (!typeof(x)=="character") {
		stop("Invalid type of x.")
		if (length(x)==0)
			stop("Length of x must be greater than zero.")
	}
	if (!is.logical(compressRepeats))
		stop("compressRepeats must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors)!=processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	if (typeof(x)=="character") {
		if (compression=="auto" || compression=="nbit") {
			y <- .Call("compress",
				x,
				ifelse(compression=="nbit", 1L, 0L),
				compressRepeats,
				processors,
				PACKAGE="DECIPHER")
			compression <- "gzip"
		} else {
			y <- rep(list(raw(0)),
				length(x))
		}
		
		if (compression != "nbit") {
			w <- which(unlist(lapply(y, length))==0)
			for (i in w)
				y[[i]] <- memCompress(x[i],
					compression)
		}
	} else { # x is a list
		y <- .Call("decompress",
			x,
			processors,
			PACKAGE="DECIPHER")
		
		bz_header <- as.raw(c(0x42, 0x5a, 0x68))
		xz_header <- as.raw(c(0xfd, 0x37, 0x7a))
		.decompress <- function(x) {
			# choose decompression type
			if (identical(x[1:3], bz_header)) {
				memDecompress(x,
					type="bzip2",
					asChar=TRUE)
			} else if (identical(x[1:3], xz_header)) {
				memDecompress(x,
					type="xz",
					asChar=TRUE)
			} else { # variable header
				memDecompress(x,
					type="gzip",
					asChar=TRUE)
			}
		}
		
		w <- which(unlist(lapply(y, is.na)))
		for (i in w)
			y[i] <- .decompress(x[[i]])
	}
	
	if (!is.null(names(x)))
		names(y) <- names(x)
	
	return(y)
}
