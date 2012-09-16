TerminalChar <- function(myDNAStringSet,
	char="-") {
	
	# error checking
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (is.na(pmatch(char, DNA_ALPHABET)))
		stop("char must be a character in the DNA_ALPHABET")
	
	if (char=="-") {
		gaps <- .Call("gaps",
			myDNAStringSet,
			PACKAGE="DECIPHER")
		dimnames(gaps) <- list(names(myDNAStringSet),
			c("leadingChar","trailingChar","difference"))
	} else {
		numF <- length(myDNAStringSet)
		maxW <- max(width(myDNAStringSet))
		
		gaps <- matrix(data=0, nrow=numF, ncol=3,
		dimnames=list(names(myDNAStringSet),
			c("leadingChar","trailingChar","difference")))
		myDNAStringSetRev <- reverse(myDNAStringSet)
		
		lead <- isMatchingAt(char, myDNAStringSet,
			seq_len(maxW))
		trail <- isMatchingAt(char, myDNAStringSetRev,
			seq_len(maxW))
		
		for (i in 1:numF) {
			# leading characters
			index <- which(!lead[,i])[1L]
			
			if (is.na(index)) { # sequence is all that character
				gaps[i,1] <- width(myDNAStringSet[i])
				gaps[i,2] <- width(myDNAStringSet[i])
				gaps[i,3] <- NA
			} else {
				gaps[i,1] <- index - 1L
				
				# trailing characters
				index <- which(!trail[,i])[1L]
				gaps[i,2] <- index - 1L
				
				# difference between leading and trailing gaps
				gaps[i,3] <- width(myDNAStringSet[i]) - gaps[i,2] - gaps[i,1]
			}
		}
	}
	return(gaps)
}