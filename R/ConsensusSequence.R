ConsensusSequence <- function(myDNAStringSet,
	threshold=.05,
	ambiguity=TRUE,
	noConsensusChar="N",
	minInformation=.75,
	ignoreNonBases=FALSE,
	includeTerminalGaps=FALSE,
	verbose=TRUE) {
	
	# error checking
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (!is.logical(ambiguity))
		stop("ambiguity must be a logical.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold > 1)
		stop("threshold cannot be greater than 1.")
	if (threshold < 0)
		stop("threshold cannot be negative.")
	if (!is.numeric(minInformation))
		stop("minInformation must be a numeric.")
	if (minInformation > 1)
		stop("minInformation cannot be greater than 1.")
	if (minInformation < 0)
		stop("minInformation cannot be negative.")
	if (is.na(pmatch(noConsensusChar, DNA_ALPHABET)))
		stop("noConsensusChar must be a character in the DNA_ALPHABET.")
	if (!is.logical(ignoreNonBases))
		stop("ignoreNonBases must be a logical.")
	
	# initialize variables
	time.1 <- Sys.time()
	maxW <- unique(width(myDNAStringSet))
	if (length(maxW)!=1 & verbose) {
		warning("\n",
			length(maxW),
			" different sequence lengths.\n",
			"End represents consensus of less sequences.\n")
	}
	
	seq <- .Call("consensusSequence",
		myDNAStringSet,
		threshold,
		ambiguity,
		minInformation,
		ignoreNonBases,
		includeTerminalGaps,
		PACKAGE="DECIPHER")
	seq[1] <- gsub("?",
		noConsensusChar,
		seq[1],
		fixed=TRUE)
	consensusSeq <- DNAStringSet(seq[1])
	
	if (verbose) {
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=1))
		cat("\n")
	}
	
	return(consensusSeq)
}