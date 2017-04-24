OrientNucleotides <- function(myXStringSet,
	reference=which.max(width(myXStringSet)),
	type="sequences",
	orientation="all",
	threshold=0.05,
	verbose=TRUE,
	processors=1) {
	
	# error checking
	if (!(is(myXStringSet, "DNAStringSet") ||
		is(myXStringSet, "RNAStringSet")))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	l <- length(myXStringSet)
	if (l==0)
		stop("myXStringSet does not contain any sequences.")
	if (all(width(myXStringSet)==0L))
		stop("All sequences in myXStringSet are zero width.")
	reference <- unique(reference)
	if (length(reference)==0)
		stop("reference must be at least one number.")
	if (any(reference < 1))
		stop("reference must be at least one.")
	if (any(reference > l))
		stop("reference greater than the length of myXStringSet.")
	if (!all(reference==floor(reference)))
		stop("reference must be a whole number.")
	if (l <= length(reference))
		stop("myXStringSet must contain non-reference sequences.")
	TYPES <- c("sequences", "orientations", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	ORIENT <- c("all", "both", "reverse", "complement")
	orientation <- unique(pmatch(orientation, ORIENT))
	if (any(is.na(orientation)))
		stop("Invalid orientation.")
	if (any(orientation == -1))
		stop("Ambiguous orientation.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold <= 0)
		stop("threshold must be greater than 0.")
	if (threshold > 1)
		stop("threshold cannot be greater than 1.")
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
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') must be removed before orienting.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') must be removed before orienting.")
	
	if (verbose) {
		time.1 <- Sys.time()
		if (1L %in% orientation) {
			tot <- 4
		} else {
			tot <- 1
			if (2L %in% orientation)
				tot <- tot + 1
			if (3L %in% orientation)
				tot <- tot + 1
			if (4L %in% orientation)
				tot <- tot + 1
		}
		count <- 1
	} else {
		pBar <- NULL
	}
	
	result <- character(l)
	result[reference] <- NA
	
	wordSize <- ceiling(log(100*mean(width(myXStringSet)), 4))
	if (wordSize > 15)
		wordSize <- 15
	if (wordSize < 2)
		wordSize <- 2
	
	v <- .Call("enumerateSequence",
		myXStringSet,
		wordSize,
		PACKAGE="DECIPHER")
	v <- lapply(v, sort.int, method="radix")
	X <- v[reference]
	
	dists <- function(x) {
		x <- .Call("matchListsDual",
			x,
			X,
			verbose,
			pBar,
			processors,
			PACKAGE="DECIPHER")
		return(1 - apply(x, 1L, max))
	}
	
	org <- numeric(l)
	if (verbose) {
		pBar <- txtProgressBar(min=0,
			max=tot*100,
			style=3)
	}
	org[-reference] <- dists(v[-reference])
	
	w <- seq_len(l)[-reference]
	if (1L %in% orientation || 2L %in% orientation) {
		seqs <- reverseComplement(myXStringSet[w])
		v <- .Call("enumerateSequence",
			seqs,
			wordSize,
			PACKAGE="DECIPHER")
		v <- lapply(v, sort.int)
		
		if (verbose) {
			pBar <- txtProgressBar(min=-100,
				max=tot*100 - 100,
				style=3)
			count <- count + 1
		}
		v <- dists(v)
		
		s <- which(v + threshold <= org[w])
		if (length(s) > 0) {
			if (type==2L || type==3L)
				result[w[s]] <- "rc"
			if (type==1L || type==3L)
				myXStringSet[w[s]] <- seqs[s]
			w <- w[-s]
		}
	}
	
	if (length(w) > 0 &&
		(1L %in% orientation || 3L %in% orientation)) {
		w <- which(result=="")
		seqs <- reverse(myXStringSet[w])
		v <- .Call("enumerateSequence",
			seqs,
			wordSize,
			PACKAGE="DECIPHER")
		v <- lapply(v, sort.int)
		
		if (verbose) {
			pBar <- txtProgressBar(min=-100*count,
				max=tot*100 - 100*count,
				style=3)
			count <- count + 1
		}
		v <- dists(v)
		
		s <- which(v + threshold <= org[w])
		if (length(s) > 0) {
			if (type==2L || type==3L)
				result[w[s]] <- "r"
			if (type==1L || type==3L)
				myXStringSet[w[s]] <- seqs[s]
			w <- w[-s]
		}
	}
	
	if (length(w) > 0 &&
		(1L %in% orientation || 4L %in% orientation)) {
		w <- which(result=="")
		seqs <- complement(myXStringSet[w])
		v <- .Call("enumerateSequence",
			seqs,
			wordSize,
			PACKAGE="DECIPHER")
		v <- lapply(v, sort.int)
		
		if (verbose) {
			pBar <- txtProgressBar(min=-100*count,
				max=tot*100 - 100*count,
				style=3)
			count <- count + 1
		}
		v <- dists(v)
		
		s <- which(v + threshold <= org[w])
		if (length(s) > 0) {
			if (type==2L || type==3L)
				result[w[s]] <- "c"
			if (type==1L || type==3L)
				myXStringSet[w[s]] <- seqs[s]
		}
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
	}
	
	if (type==1L) {
		return(myXStringSet)
	} else if (type==2L) {
		return(result)
	} else { # type==3L
		return(list(result,
			myXStringSet))
	}
}
