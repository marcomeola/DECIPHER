DistanceMatrix <- function(myXStringSet,
	includeTerminalGaps=FALSE,
	penalizeGapLetterMatches=TRUE,
	penalizeGapGapMatches=FALSE,
	correction="none",
	processors=1,
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
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is(myXStringSet, "XStringSet"))
		stop("myXStringSet must be an XStringSet.")
	if (is(myXStringSet, "BStringSet"))
		stop("myXStringSet cannot be a BStringSet.")
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
	
	maxW <- unique(width(myXStringSet))
	if (length(maxW)!=1) {
		if (verbose)
			warning("\n",
				length(maxW),
				" different sequence lengths.\n",
				"Using shorter length in each comparison.\n")
	}
	numF <- length(myXStringSet)
	if (numF < 2) {
		stop("Too few sequences!")
	}
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
	} else {
		pBar <- NULL
	}
	
	# calculate the distance matrix
	distMatrix <- .Call("distMatrix",
		myXStringSet,
		ifelse(is(myXStringSet, "AAStringSet"), 3L, 1L),
		includeTerminalGaps,
		penalizeGapGapMatches,
		penalizeGapLetterMatches,
		TRUE, # full matrix
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	dimnames(distMatrix) <- list(names(myXStringSet),
		names(myXStringSet))
	
	# apply distance correction
	if (correction==2) { # Jukes-Cantor
		if (is(myXStringSet, "DNAStringSet") || is(myXStringSet, "RNAStringSet")) {
			# JC func is undefined above p = .75
			w <- which(distMatrix > .75)
			if (length(w) > 0) {
				# rather than produce NaNs
				# make the distance infinite
				distMatrix[w] <- .75
			}
			distMatrix <- -3/4*log(1 - 4/3*distMatrix)
			attr(distMatrix, "correction") <- "Jukes-Cantor"
		} else if (is(myXStringSet, "AAStringSet")) {
			# JC func is undefined above p = .95
			w <- which(distMatrix > .95)
			if (length(w) > 0) {
				# rather than produce NaNs
				# make the distance infinite
				distMatrix[w] <- .95
			}
			distMatrix <- -19/20*log(1 - 20/19*distMatrix)
			attr(distMatrix, "correction") <- "Jukes-Cantor"
		} else {
			warning("Jukes-Cantor correction is not available for a BStringSet input.")
			attr(distMatrix, "correction") <- "none"
		}
	} else {
		attr(distMatrix, "correction") <- "none"
	}
	
	if (verbose) {
		close(pBar)
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
