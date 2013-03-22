DistanceMatrix <- function(myDNAStringSet,
	includeTerminalGaps=FALSE,
	penalizeGapLetterMatches=TRUE,
	penalizeGapGapMatches=FALSE,
	removeDuplicates=FALSE,
	correction="none",
	verbose=TRUE) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	CORRECTIONS <- c("none", "Jukes-Cantor", "JC")
	correction <- pmatch(correction, CORRECTIONS)
	if (is.na(correction))
		stop("Invalid distance correction method.")
	if (correction == -1)
		stop("Ambiguous distance correction method.")
	if (correction==3)
		correction <- 2
	if (!is.logical(includeTerminalGaps))
		stop("includeTerminalGaps must be a logical.")
	if (!is.logical(penalizeGapGapMatches))
		stop("penalizeGapGapMatches must be a logical.")
	if (!is.logical(penalizeGapLetterMatches))
		stop("penalizeGapLetterMatches must be a logical.")
	if (!is.logical(removeDuplicates))
		stop("removeDuplicates must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	
	maxW <- unique(width(myDNAStringSet))
	if (length(maxW)!=1) {
		if (verbose)
			warning("\n",
				length(maxW),
				" different sequence lengths.\n",
				"Using shorter length in each comparison.\n")
	}
	numF <- length(myDNAStringSet)
	if (numF < 2) {
		stop("Too few sequences!")
	}
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
	} else {
		pBar <- NULL
	}
	
	if (removeDuplicates) {
		nRemoved <- length(myDNAStringSet)
		myDNAStringSet <- unique(myDNAStringSet)
		nRemoved <- nRemoved - length(myDNAStringSet)
	}
	
	# calculate the distance matrix
	distMatrix <- .Call("distMatrix",
		myDNAStringSet,
		includeTerminalGaps,
		penalizeGapGapMatches,
		penalizeGapLetterMatches,
		verbose,
		pBar,
		PACKAGE="DECIPHER")
	dimnames(distMatrix) <- list(names(myDNAStringSet),
		names(myDNAStringSet))
	
	# apply distance correction
	if (correction==2) { # Jukes-Cantor
		# JC func is undefined above p = .75
		w <- which(distMatrix > .75)
		if (length(w) > 0) {
			# rather than produce NaNs
			# make the distance infinite
			distMatrix[w] <- .75
		}
		distMatrix <- -3/4*log(1 - 4/3*distMatrix)
		attr(distMatrix, "correction") <- "Jukes-Cantor"
	} else {
		attr(distMatrix, "correction") <- "none"
	}
	
	if (verbose) {
		close(pBar)
		#if (removeDuplicates)
		#	cat("\nRemoved", nRemoved, "exact duplicate sequences.")
		#if (includeTerminalGaps==TRUE & penalizeGapGapMatches==FALSE)
		#	cat("\nCompared the union of internal ranges.",
		#		"\nGap-gap matches not included in distance.")
		#if (includeTerminalGaps==FALSE & penalizeGapGapMatches==FALSE)
		#	cat("\nCompared the intersection of internal ranges.",
		#		"\nGap-gap matches not included in distance.")
		#if (includeTerminalGaps==TRUE & penalizeGapGapMatches==TRUE)
		#	cat("\nCompared the entire sequence length.",
		#		"\nConsidered gap-gap matches as mis-matches.")
		#if (includeTerminalGaps==FALSE & penalizeGapGapMatches==TRUE)
		#	cat("\nCompared the intersection of internal ranges.",
		#		"\nConsidered gap-gap matches as mis-matches.")
		#if (penalizeGapLetterMatches)
		#	cat("\nConsidered gap-letter matches as mis-matches.")
		#else
		#	cat("\nGap-letter matches not included in distance.")
		#if (correction==2)
		#	cat("\nApplied Jukes-Cantor correction.")
		#else if (correction==3)
		#	cat("\nApplied Huber-Hugenholtz correction.")
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	return(distMatrix)
}