Codec <- function(x,
	compression=c("nbit", "gzip"),
	compressRepeats=FALSE,
	processors=1) {
	
	# error checking
	if (length(compression)==2) {
		if (compression[1] != "nbit" && compression[1] != "qbit")
			stop("The first element of compression must be 'nbit' or 'qbit' when two elements are provided.")
		if (!(compression[2] %in% c("gzip", "bzip2", "xz")))
			stop("The second element of compression must be  'gzip', 'bzip2', or 'xz' when two elements are provided.")
	} else if (length(compression) != 1 || !is.character(compression)) {
		stop("compression must be a character string.")
	} else if (!(compression %in% c("nbit", "qbit", "gzip", "bzip2", "xz"))) {
		stop("Invalid type of compression.")
	}
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
		if (compression[1]=="nbit") {
			y <- .Call(compression[1],
				x,
				2L - length(compression),
				compressRepeats,
				processors,
				PACKAGE="DECIPHER")
			if (length(compression)==2) {
				w <- which(lengths(y)==0)
				for (i in w)
					y[[i]] <- memCompress(x[i],
						compression[2])
			}
		} else if (compression[1]=="qbit") {
			y <- .Call(compression[1],
				x,
				2L - length(compression),
				processors,
				PACKAGE="DECIPHER")
			if (length(compression)==2) {
				w <- which(lengths(y)==0)
				for (i in w)
					y[[i]] <- memCompress(x[i],
						compression[2])
			}
		} else {
			y <- lapply(x,
				memCompress,
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
