CorrectFrameshifts <- function(myXStringSet,
	myAAStringSet,
	type="indels",
	acceptDistance=0.01,
	rejectDistance=0.60,
	maxComparisons=10,
	gapOpening=-13,
	gapExtension=-1,
	frameShift=-15,
	geneticCode=GENETIC_CODE,
	substitutionMatrix="MIQS",
	verbose=TRUE,
	processors=1) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	if (!is(myAAStringSet, "AAStringSet"))
		stop("myAAStringSet must be an AAStringSet.")
	if (length(myXStringSet)==0)
		stop("At least one sequence is required in myXStringSet.")
	if (length(myAAStringSet)==0)
		stop("At least one sequence is required in myAAStringSet.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("+", myXStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("-", myAAStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') in myAAStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("+", myAAStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') in myAAStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern(".", myAAStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') in myAAStringSet must be removed before correcting frameshifts.")
	myAAStringSet <- unique(myAAStringSet)
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
	TYPES <- c("indels", "sequences", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(frameShift))
		stop("frameShift must be a numeric.")
	if (!is.numeric(rejectDistance))
		stop("rejectDistance must be a numeric.")
	if (rejectDistance > 1)
		stop("rejectDistance can be at most 1.")
	if (rejectDistance <= 0)
		stop("rejectDistance must be greater than zero.")
	if (!is.numeric(acceptDistance))
		stop("acceptDistance must be a numeric.")
	if (acceptDistance > rejectDistance)
		stop("acceptDistance can be at most rejectDistance.")
	if (acceptDistance < 0)
		stop("acceptDistance must be at least zero.")
	if (!is.numeric(maxComparisons))
		stop("maxComparisons must be a numeric.")
	if (maxComparisons < 1)
		stop("maxComparisons must be at least one.")
	if (floor(maxComparisons)!=maxComparisons)
		stop("maxComparisons must be a whole number.")
	AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
	if (is.character(substitutionMatrix)) {
		if (substitutionMatrix=="MIQS") {
			subMatrix <- matrix(c(3.2,-1.3,-0.4,-0.4,1.5,-0.2,-0.4,0.4,-1.2,-1.3,-1.4,-0.7,-1,-2.3,-0.1,0.8,0.8,-3.6,-2.4,0,-6.1,-1.3,6.2,-0.1,-1.5,-2.7,1.8,-0.7,-1.9,0.9,-2.4,-2.5,3.3,-1.1,-3.3,-1.1,-0.3,-0.9,-3.8,-1.9,-2.3,-6.1,-0.4,-0.1,5.1,2.6,-1.6,0.9,0.8,0.2,1,-3.6,-3.5,0.7,-2.3,-3.5,-1.4,0.9,0,-4.5,-1.5,-2.6,-6.1,-0.4,-1.5,2.6,5.7,-3.7,0.9,2.7,-0.5,0.3,-4.5,-4.6,0.4,-3.3,-5.8,-0.3,0.3,-0.2,-5.3,-3.9,-3.5,-6.1,1.5,-2.7,-1.6,-3.7,11.7,-2.8,-3.2,-1.7,-1.2,0.2,-2.3,-3.2,0.1,-2.8,-2.8,1,0,-6.1,-0.7,1.8,-6.1,-0.2,1.8,0.9,0.9,-2.8,3.6,2.1,-1.6,1.2,-2.2,-1.9,1.7,-0.4,-2.4,-0.4,0.4,0.1,-5.4,-2.8,-1.8,-6.1,-0.4,-0.7,0.8,2.7,-3.2,2.1,4.3,-1.3,-0.2,-3.3,-2.8,1.1,-2.3,-4.1,0,0.4,-0.2,-5.8,-2.4,-2.3,-6.1,0.4,-1.9,0.2,-0.5,-1.7,-1.6,-1.3,7.6,-1.6,-5.4,-4.8,-1.7,-3.6,-4.6,-1.6,0,-1.9,-4.8,-4.5,-3.8,-6.1,-1.2,0.9,1,0.3,-1.2,1.2,-0.2,-1.6,7.5,-2.2,-1.9,0,-2.1,0,-1.5,0,-0.2,-0.3,2.1,-2.3,-6.1,-1.3,-2.4,-3.6,-4.5,0.2,-2.2,-3.3,-5.4,-2.2,4.6,3.1,-2.3,1.7,0.7,-3.7,-2.8,-0.7,-0.7,-0.8,3.3,-6.1,-1.4,-2.5,-3.5,-4.6,-2.3,-1.9,-2.8,-4.8,-1.9,3.1,4.6,-2.4,3.2,2.1,-2.8,-2.9,-1.6,-0.2,0,2,-6.1,-0.7,3.3,0.7,0.4,-3.2,1.7,1.1,-1.7,0,-2.3,-2.4,3.6,-1.1,-3.7,-0.1,0,0,-4,-2.3,-2,-6.1,-1,-1.1,-2.3,-3.3,0.1,-0.4,-2.3,-3.6,-2.1,1.7,3.2,-1.1,5.4,1.4,-2.8,-1.8,-0.8,-2.1,-0.9,1.4,-6.1,-2.3,-3.3,-3.5,-5.8,-2.8,-2.4,-4.1,-4.6,0,0.7,2.1,-3.7,1.4,7.4,-3.7,-2.6,-2.3,4.2,5.2,-0.3,-6.1,-0.1,-1.1,-1.4,-0.3,-2.8,-0.4,0,-1.6,-1.5,-3.7,-2.8,-0.1,-2.8,-3.7,8.4,-0.1,-0.5,-3.6,-4.5,-2.5,-6.1,0.8,-0.3,0.9,0.3,1,0.4,0.4,0,0,-2.8,-2.9,0,-1.8,-2.6,-0.1,3.1,1.6,-3.5,-1.5,-1.4,-6.1,0.8,-0.9,0,-0.2,0,0.1,-0.2,-1.9,-0.2,-0.7,-1.6,0,-0.8,-2.3,-0.5,1.6,3.8,-5.3,-2.1,-0.1,-6.1,-3.6,-3.8,-4.5,-5.3,-6.1,-5.4,-5.8,-4.8,-0.3,-0.7,-0.2,-4,-2.1,4.2,-3.6,-3.5,-5.3,14.8,4.9,-3.3,-6.1,-2.4,-1.9,-1.5,-3.9,-0.7,-2.8,-2.4,-4.5,2.1,-0.8,0,-2.3,-0.9,5.2,-4.5,-1.5,-2.1,4.9,8.3,-1.2,-6.1,0,-2.3,-2.6,-3.5,1.8,-1.8,-2.3,-3.8,-2.3,3.3,2,-2,1.4,-0.3,-2.5,-1.4,-0.1,-3.3,-1.2,3.5,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,1),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
	"PAM30", "PAM40", "PAM70", "PAM120", "PAM250")) {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix, envir=environment())))
		} else {
			stop("Invalid substitutionMatrix.")
		}
	} else if (is.matrix(substitutionMatrix)) {
		if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
			any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
			stop("substitutionMatrix is incomplete.")
		subMatrix <- substitutionMatrix
	} else {
		stop("Invalid substitutionMatrix.")
	}
	subMatrix <- subMatrix[AAs, AAs] + 0 # convert to numeric
	
	# de-replicate
	w <- which(!duplicated(myXStringSet))
	l <- length(myXStringSet)
	ns <- names(myXStringSet)
	if (length(w) < l) {
		dupes <- TRUE
		derep <- numeric(l)
		derep[w] <- seq_along(w)
		derep[-w] <- match(myXStringSet[-w],
			myXStringSet[w])
		myXStringSet <- myXStringSet[w]
		l <- length(w)
	} else {
		dupes <- FALSE
	}
	
	ends <- width(myXStringSet)
	if (min(ends) < 2)
		stop("All sequences in myXStringSet must be at least two nucleotides long.")
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Finding the closest reference amino acid sequences:\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(max=100, style=3)
	} else {
		pBar <- NULL
	}
	
	frames <- AAStringSet()
	for (i in 1:3) {
		end <- ends
		offset <- end - i + 1
		end <- end - offset %% 3
		end <- ifelse(end < i - 1,
			i - 1,
			end)
		
		AA <- translate(subseq(myXStringSet,
				i,
				end),
			genetic.code=geneticCode,
			if.fuzzy.codon="solve")
		frames <- c(frames,
			AA)
	}
	
	v1 <- .Call("enumerateSequenceAA",
		frames,
		7L,
		PACKAGE="DECIPHER")
	v2 <- .Call("enumerateSequenceAA",
		myAAStringSet,
		7L,
		PACKAGE="DECIPHER")
	
	v1 <- lapply(v1,
		sort.int,
		method="quick")
	v2 <- lapply(v2,
		sort.int,
		method="quick")
	
	d <- .Call("matchListsDual",
		v1,
		v2,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	d <- rowsum(d,
		rep(seq_len(l),
			3),
		na.rm=TRUE)
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		
		time.1 <- Sys.time()
		cat("\nAssessing frameshifts in nucleotide sequences:\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(max=100, style=3)
	}
	
	widths <- width(myAAStringSet)
	
	index <- matrix(nrow=nrow(d), ncol=ncol(d))
	for (i in seq_len(nrow(d)))
		index[i,] <- order(d[i,], widths, decreasing=TRUE)
	if (ncol(index) > maxComparisons) {
		index <- index[, seq_len(maxComparisons), drop=FALSE]
	} else {
		maxComparisons <- ncol(index)
	}
	
	X <- .Call("findFrameshifts",
		myAAStringSet,
		ends,
		frames,
		index,
		maxComparisons,
		gapOpening,
		gapExtension,
		frameShift,
		acceptDistance,
		rejectDistance,
		subMatrix,
		verbose,
		pBar,
		PACKAGE="DECIPHER")
	if (type != 2) {
		X <- lapply(X,
			setNames,
			c("insertions", "deletions", "distance", "index"))
	}
	
	if (type > 1) {
		pos <- lapply(X,
			function(x) {
				ins <- rle(x[[1]])
				dels <- rle(x[[2]])
				IRanges(start=c(ins[[2]],
						dels[[2]]),
					width=c(ins[[1]],
						rep(0,
							length(dels[[2]]))))
			})
		pos <- unname(pos)
		
		val <- lapply(X,
			function(x) {
				ins <- rep("",
					length(unique(x[[1]])))
				dels <- rle(x[[2]])
				types <- c("N", "NN")
				dels <- types[dels[[1]]]
				c(ins, dels)
			})
		
		if (type==2L) { # sequences
			X <- replaceAt(myXStringSet,
				RangesList(pos),
				val)
		} else { # both
			X <- list(indels=X,
				sequences=replaceAt(myXStringSet,
					RangesList(pos),
					val))
		}
	}
	
	if (dupes) { # re-replicate
		if (type==1 || type==2) {
			X <- X[derep]
		} else {
			X[[1]] <- X[[1]][derep]
			X[[2]] <- X[[2]][derep]
		}
	}
	
	if (type==1 || type==2) {
		names(X) <- ns
	} else {
		names(X[[1]]) <- names(X[[2]]) <- ns
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
	
	return(X)
}
