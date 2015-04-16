TerminalChar <- function(myXStringSet,
	char="") {
	
	# error checking
	if (!is(myXStringSet, "XStringSet"))
		stop("myXStringSet must be an XStringSet.")
		
	
	if (char=="") {
		if (is(myXStringSet, "BStringSet"))
			stop("A single character must be specified with a BStringSet input.")
		gaps <- .Call("gaps",
			myXStringSet,
			ifelse(is(myXStringSet, "AAStringSet"), 3L, 1L),
			PACKAGE="DECIPHER")
		dimnames(gaps) <- list(names(myXStringSet),
			c("leadingChar","trailingChar","difference"))
	} else {
		numF <- length(myXStringSet)
		maxW <- max(width(myXStringSet))
		
		gaps <- matrix(data=0, nrow=numF, ncol=3,
		dimnames=list(names(myXStringSet),
			c("leadingChar","trailingChar","difference")))
		myXStringSetRev <- reverse(myXStringSet)
		
		lead <- isMatchingAt(char, myXStringSet,
			seq_len(maxW))
		trail <- isMatchingAt(char, myXStringSetRev,
			seq_len(maxW))
		
		for (i in 1:numF) {
			# leading characters
			index <- which(!lead[,i])[1L]
			
			if (is.na(index)) { # sequence is all that character
				gaps[i,1] <- width(myXStringSet[i])
				gaps[i,2] <- width(myXStringSet[i])
				gaps[i,3] <- NA
			} else {
				gaps[i,1] <- index - 1L
				
				# trailing characters
				index <- which(!trail[,i])[1L]
				gaps[i,2] <- index - 1L
				
				# difference between leading and trailing gaps
				gaps[i,3] <- width(myXStringSet[i]) - gaps[i,2] - gaps[i,1]
			}
		}
	}
	return(gaps)
}
