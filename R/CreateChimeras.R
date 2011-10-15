CreateChimeras <- function(myDNAStringSet,
	numChimeras=10,
	numParts=2,
	minLength=80,
	maxLength=Inf,
	minChimericRegionLength=30,
	randomLengths=TRUE,
	includeParents=TRUE,
	verbose=TRUE) {
	
	# error checking
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(includeParents))
		stop("includeParents must be a logical.")
	if (!is.logical(randomLengths))
		stop("randomLengths must be a logical.")
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (!is.numeric(numChimeras))
		stop("numChimeras must be a numeric.")
	if (floor(numChimeras)!=numChimeras)
		stop("numChimeras must be a whole number.")
	if (floor(numChimeras) < 1)
		stop("numChimeras must be greater than zero.")
	if (!is.numeric(numParts))
		stop("numParts must be a numeric.")
	if (floor(numParts)!=numParts)
		stop("numChimeras must be a whole number.")
	if (floor(numParts) < 2)
		stop("numParts must be greater than one.")
	if (!is.numeric(minLength))
		stop("minLength must be a numeric.")
	if (floor(minLength)!=minLength)
		stop("minLength must be a whole number.")
	if (minLength < 0)
		stop("minLength must be greater than or equal to zero.")
	if (!is.numeric(minChimericRegionLength))
		stop("minChimericRegionLength must be a numeric.")
	if (floor(minChimericRegionLength)!=minChimericRegionLength)
		stop("minChimericRegionLength must be a whole number.")
	if (minChimericRegionLength < 0)
		stop("minChimericRegionLength must be greater than or equal to zero.")
	if (!is.numeric(maxLength))
		stop("maxLength must be a numeric.")
	if (floor(maxLength)!= maxLength)
		stop("maxLength must be a whole number.")
	if (maxLength < 0)
		stop("maxLength must be greater than or equal to zero.")
	
	# initialize variables
	time.1 <- Sys.time()
	
	# calculate number of sequences
	myDNAStringSet <- unique(myDNAStringSet)
	uw <- unique(width(myDNAStringSet))
	if (length(uw) > 1)
		stop("All sequences must have the same width (be aligned).")
	numF <- length(myDNAStringSet)
	if (numF < numParts)
		stop("myDNAStringSet must contain more unique sequences.")
	if (numChimeras > numF/numParts)
		stop("Cannot make more chimeras than half of the number of unique myDNAStringSet1 sequences.")
	if (is.null(names(myDNAStringSet)))
		names(myDNAStringSet) <- 1:length(myDNAStringSet)
	
	# initialize a progress bar
	if (verbose)
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
	
	# initialize variables
	i <- 0
	chimeras <- DNAStringSet()
	if (includeParents)
		parents <- DNAStringSet()
	s <- sample(numF)
	for (p in 1:numParts) {
		expr1 <- paste("s",
			p,
			" <- s[",
			((p - 1)*floor(length(s)/numParts) + 1),
			":",
			(p*floor(length(s)/numParts)),
			"]",
			sep="")
		eval(parse(text=expr1))
	}
	for (p in 1:numParts) {
		expr1 <- paste("t",
			p,
			" <- TerminalChar(myDNAStringSet[s",
			p,
			"])",
			sep="")
		eval(parse(text=expr1))
	}
	
	# splice sequences together
	while (length(chimeras) < numChimeras &&
		(i + 1) <= length(s1)) {
		# increment counter
		i <- i + 1
		
		# generate subsequences
		for (p in 1:numParts) {
			expr1 <- paste("seq",
				p,
				" <- myDNAStringSet[[s",
				p,
				"[",
				i,
				"]]]",
				sep="")
			eval(parse(text=expr1))
		}
		
		for (p in 1:numParts) {
			expr1 <- paste("width",
				p,
				" <- length(seq",
				p,
				")",
				sep="")
			eval(parse(text=expr1))
		}
		
		bp <- sample(2:(uw - 1), numParts - 1) # break points
		bp <- bp[order(bp)]
		bp <- c(1, bp, uw + 1)
		
		for (p in 1:numParts) {
			expr1 <- paste("subseq",
				p,
				" <- subseq(seq",
				p,
				", ",
				bp[p],
				", ",
				bp[p + 1] - 1,
				")",
				sep="")
			eval(parse(text=expr1))
		}
		
		# merge subsequences together
		expr1 <- "chimera <- c(subseq1"
		for (p in 2:numParts) {
			expr1 <- paste(expr1,
				", subseq",
				p,
				sep="")
		}
		expr1 <- paste(expr1, ")", sep="")
		eval(parse(text=expr1))
		
		if (randomLengths) {
			l1 <- sample(2:bp[2] - 1, 1)
			l2 <- sample((bp[numParts] + 1):(uw - 1), 1)
			if (l1 > 0)
				chimera <- c(DNAString(paste(rep("-",l1 - 1), collapse="")),
					subseq(chimera, l1, uw))
			chimera <- c(subseq(chimera, 1, l2),
				DNAString(paste(rep("-", uw - l2), collapse="")))
		}
		
		chimera <- DNAStringSet(chimera)
		aF <- sum(alphabetFrequency(chimera)[,1:15])
		if ((aF < minLength) || (aF > maxLength))
			next
		
		expr1 <- "d <- DistanceMatrix(c(DNAStringSet(seq1"
		for (p in 2:numParts) {
			expr1 <- paste(expr1,
				"), DNAStringSet(seq",
				p,
				sep="")
		}
		expr1 <- paste(expr1, "), chimera), verbose=FALSE)", sep="")
		eval(parse(text=expr1))
		
		add2set <- TRUE
		for (p in 1:numParts) {
			if (is.na(d[numParts + 1, p]))
				add2set <- FALSE
		}
		if (!add2set)
			next
		
		for (p in 1:numParts) {
			if (d[numParts + 1, p]==0) {
				add2set <- FALSE
				break
			}
			
			expr1 <- paste("myName <- names(myDNAStringSet)[s",
				p,
				"[",
				i,
				"]]",
				sep="")
			eval(parse(text=expr1))
			if (p==1) {
				start <- ifelse(randomLengths,
					ifelse(l1==0, 1, l1), 1)
				end <- bp[p + 1] - 1
				if (start > end) {
					add2set <- FALSE
					break
				}
				regionLength <- sum(alphabetFrequency(subseq(chimera,
					start,
					end))[,1:15])
				
				if (regionLength < minChimericRegionLength) {
					add2set <- FALSE
					break
				}
				
				myNames <- paste(myName,
					" [",
					ifelse(randomLengths, ifelse(l1==0, 1, l1), bp[p]),
					"-",
					bp[p + 1] - 1,
					",",
					regionLength,
					"] (",
					round(100*d[numParts + 1, p],digits=1),
					"%)",
					sep="")
			} else {
				start <- bp[p]
				end <- ifelse((p==numParts) && randomLengths,
					l2, bp[p + 1] - 1)
				if (start > end) {
					add2set <- FALSE
					break
				}
				regionLength <- sum(alphabetFrequency(subseq(chimera,
					start,
					end))[,1:15])
				
				if (regionLength < minChimericRegionLength) {
					add2set <- FALSE
					break
				}
				
				myNames <- paste(myNames,
					", ",
					myName,
					" [",
					bp[p],
					"-",
					ifelse((p==numParts) && randomLengths,
						l2, bp[p + 1] - 1),
					",",
					regionLength,
					"] (",
					round(100*d[numParts + 1, p],digits=1),
					"%)",
					sep="")
			}
		}
		if (!add2set)
			next
		
		names(chimera) <- myNames
		
		# add another chimera to sequence set
		chimeras <- c(chimeras,
			chimera)
		
		if (includeParents) {
			for (p in 1:numParts) {
				if (randomLengths) {
					expr1 <- paste("parents",
						p,
						" <- DNAStringSet(c(DNAString(paste(rep('-', ",
						l1 - 1,
						"), collapse='')), subseq(myDNAStringSet[[s",
						p,
						"[",
						i,
						"]]], l1, l2), DNAString(paste(rep('-', ",
						uw - l2,
						"), collapse=''))))",
						sep="")
					eval(parse(text=expr1))
				} else {
					expr1 <- paste("parents",
						p,
						" <- myDNAStringSet[[s",
						p,
						"[",
						i,
						"]]]",
						sep="")
					eval(parse(text=expr1))
				}
			}
			
			if (randomLengths) {
				# need to recalculate distance matrix with shorter parents
				expr1 <- "d <- DistanceMatrix(c(parents1"
				for (p in 2:numParts) {
					expr1 <- paste(expr1,
						", parents",
						p,
						sep="")
				}
			}
			expr1 <- paste(expr1, "), verbose=FALSE)", sep="")
			eval(parse(text=expr1))
			
			for (p in 1:numParts) {
				expr1 <- paste("names(parents",
					p,
					") <- paste(names(myDNAStringSet[s",
					p,
					"[",
					i,
					"]]), '[",
					l1,
					"-",
					l2,
					"] (",
					paste(round(100*d[p,][1:numParts][-p],
							digits=1),
						collapse="%, "),
					"%)')",
					sep="")
				eval(parse(text=expr1))
			}
			expr1 <- "parents <- c(parents"
			for (p in 1:numParts) {
				expr1 <- paste(expr1,
					", parents",
					p,
					sep="")
			}
			expr1 <- paste(expr1, ")", sep="")
			eval(parse(text=expr1))
		}
		
		if (verbose)
			setTxtProgressBar(pBar,
				floor(100*length(chimeras)/numChimeras))
	}
	
	if (verbose) {
		if (length(chimeras) < numChimeras)
			warning("Unable to create the requested number of chimeras.")
		close(pBar)
		cat("\nCreated ",
			length(chimeras),
			" chimeras.",
			sep="")
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=1))
		cat("\n")
	}
	
	if (includeParents)
		chimeras <- c(chimeras, parents)
	
	return(chimeras)
}