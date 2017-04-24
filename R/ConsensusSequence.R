ConsensusSequence <- function(myXStringSet,
	threshold=0.05,
	ambiguity=TRUE,
	noConsensusChar="+",
	minInformation=1 - threshold,
	ignoreNonBases=FALSE,
	includeTerminalGaps=FALSE) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet") && !is(myXStringSet, "AAStringSet"))
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	if (!is.logical(ambiguity))
		stop("ambiguity must be a logical.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold >= 1)
		stop("threshold must be less than one.")
	if (threshold < 0)
		stop("threshold cannot be negative.")
	if (!is.numeric(minInformation))
		stop("minInformation must be a numeric.")
	if (minInformation > 1)
		stop("minInformation cannot be greater than one.")
	if (minInformation <= 0)
		stop("minInformation must be greater than zero.")
	if (is(myXStringSet, "DNAStringSet") && is.na(pmatch(noConsensusChar, DNA_ALPHABET)))
		stop("noConsensusChar must be a character in the DNA_ALPHABET.")
	if (is(myXStringSet, "RNAStringSet") && is.na(pmatch(noConsensusChar, RNA_ALPHABET)))
		stop("noConsensusChar must be a character in the RNA_ALPHABET.")
	if (is(myXStringSet, "AAStringSet") && is.na(pmatch(noConsensusChar, AA_ALPHABET)))
		stop("noConsensusChar must be a character in the AA_ALPHABET.")
	if (!is.logical(ignoreNonBases))
		stop("ignoreNonBases must be a logical.")
	
	if (is(myXStringSet, "AAStringSet")) {
		seq <- .Call("consensusSequenceAA",
			myXStringSet,
			threshold,
			ambiguity,
			minInformation,
			ignoreNonBases,
			includeTerminalGaps,
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		seq <- .Call("consensusSequence",
			myXStringSet,
			threshold,
			ambiguity,
			minInformation,
			ignoreNonBases,
			includeTerminalGaps,
			PACKAGE="DECIPHER")
		if (is(myXStringSet, "RNAStringSet"))
			seq[1] <- gsub("T",
				"U",
				seq[1],
				fixed=TRUE)
	}
	
	seq[1] <- gsub("?",
		noConsensusChar,
		seq[1],
		fixed=TRUE)
	if (is(myXStringSet, "DNAStringSet")) {
		consensusSeq <- DNAStringSet(seq[1])
	} else if (is(myXStringSet, "RNAStringSet")) {
		consensusSeq <- RNAStringSet(seq[1])
	} else { # AAStringSet
		consensusSeq <- AAStringSet(seq[1])
	}
	
	return(consensusSeq)
}
