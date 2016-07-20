TrimDNA <- function(myDNAStringSet,
	leftPatterns,
	rightPatterns,
	type="ranges",
	quality=NULL,
	maxDistance=0.1,
	minOverlap=5,
	allowInternal=TRUE,
	alpha=0.2,
	threshold=0.01,
	minWidth=36,
	verbose=TRUE) {
	
	# error checking:
	if (class(myDNAStringSet) != "DNAStringSet")
		stop("myDNAStringSet must be a DNAStringSet.")
	if (!(class(quality) %in% c("NULL",
		"PhredQuality",
		"SolexaQuality",
		"IlluminaQuality")))
		stop("quality must be a PhredQuality, SolexaQuality, or IlluminaQuality object.")
	if (!is.null(quality) &&
		any(width(myDNAStringSet) != width(quality)))
		stop("The widths of myDNAStringSet must match the widths of quality.")
	if (length(leftPatterns)==0)
		stop("leftPatterns must be at least length one.")
	if (class(leftPatterns) %in% c("DNAString", "DNAStringSet")) {
		leftPattern <- as.character(leftPatterns)
	} else if (!is.character(leftPatterns)) {
		stop("leftPatterns must be a DNAStringSet or character vector.")
	}
	if (length(rightPatterns)==0)
		stop("rightPatterns must be at least length one.")
	if (class(rightPatterns) %in% c("DNAString", "DNAStringSet")) {
		rightPattern <- as.character(rightPatterns)
	} else if (!is.character(rightPatterns)) {
		stop("rightPatterns must be a DNAStringSet or character vector.")
	}
	TYPES <- c("ranges", "sequences", "both")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (!is.numeric(maxDistance))
		stop("maxDistance must be a numeric.")
	if (maxDistance >= 1 || maxDistance < 0)
		stop("maxDistance must be between zero and one.")
	if (!is.numeric(minOverlap))
		stop("minOverlap must be a numeric.")
	if (minOverlap < 1)
		stop("minOverlap must be at least one.")
	if (!is.logical(allowInternal))
		stop("allowInternal must be a logical.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold >= 1 || threshold <= 0)
		stop("threshold must be between zero and one.")
	if (!is.numeric(alpha))
		stop("alpha must be a numeric.")
	if (alpha > 1 || alpha <= 0)
		stop("alpha must be between zero and one.")
	if (!is.numeric(minWidth))
		stop("minWidth must be a numeric.")
	if (minWidth < 0)
		stop("minWidth must be at least zero.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	
	if (verbose)
		time.1 <- Sys.time()
	
	.findPatternEnd <- function(pattern) {
		left <- numeric(length(myDNAStringSet))
		l <- nchar(pattern)
		
		if (allowInternal) {
			# find interal matches without indels
			m <- vmatchPattern(pattern,
				myDNAStringSet,
				max.mismatch=floor(l*maxDistance),
				fixed="subject")
			s <- endIndex(m)
			x <- which(lengths(s) > 0)
			if (length(x) > 0) {
				left[x] <- sapply(s[x],
					function(x) {
						if (length(x) > 1) {
							max(x)
						} else {
							x
						}
					})
			}
			
			if (verbose) {
				n <- length(x)
				cat(round(n/length(left)*100, 1),
					"% internal, ",
					sep="")
				flush.console()
			}
			
			x <- which(lengths(s)==0)
		} else {
			x <- seq_len(length(myDNAStringSet))
			n <- 0
		}
		
		if (l > minOverlap) {
			s <- substring(pattern, seq_len(l - minOverlap) + 1)
			for (i in seq_along(s)) {
				if (length(x)==0)
					break
				w <- which(isMatchingAt(s[i],
					.subset(myDNAStringSet, x),
					max.mismatch=floor((l - i + 1)*maxDistance),
					with.indels=TRUE,
					fixed="subject"))
				if (length(w) > 0) {
					left[x[w]] <- l - i
					x <- x[-w]
				}
			}
		}
		
		if (verbose) {
			n <- length(left) - length(x) - n
			cat(round(n/length(left)*100, 1),
				"% flanking",
				sep="")
			flush.console()
		}
		
		return(left)
	}
	
	# try to find the flanking left pattern
	lefts <- integer(length(myDNAStringSet))
	lefts[] <- 1L
	for (k in which(nchar(leftPatterns) >= minOverlap)) {
		if (verbose) {
			cat(ifelse(k==1, "", "\n"),
				"Finding left pattern",
				ifelse(length(leftPatterns) > 1,
					paste(" (#", k, "): ", sep=""),
					": "),
				sep="")
			flush.console()
		}
		left <- .findPatternEnd(leftPatterns[k]) + 1L
		
		w <- which(left > lefts)
		if (length(w) > 0)
			lefts[w] <- left[w]
	}
	
	# reverse the patterns and subjects
	dna <- myDNAStringSet
	myDNAStringSet <- reverse(myDNAStringSet)
	rightPatterns <- sapply(strsplit(rightPatterns,
			"",
			fixed=TRUE),
		function(x) {
			paste(rev(x), collapse="")
		})
	
	# try to find the flanking right pattern
	rights <- width(myDNAStringSet)
	for (k in which(nchar(rightPatterns) >= minOverlap)) {
		if (verbose) {
			cat("\nFinding right pattern",
				ifelse(length(rightPatterns) > 1,
					paste(" (#", k, "): ", sep=""),
					": "),
				sep="")
			flush.console()
		}
		right <- width(myDNAStringSet) - .findPatternEnd(rightPatterns[k])
		
		w <- which(right < rights)
		if (length(w) > 0)
			rights[w] <- right[w]
	}
	
	if (!is.null(quality)) {
		if (verbose) {
			cat("\nTrimming by quality score: ")
			flush.console()
		}
		bounds <- .Call("movAvg",
			quality,
			match(class(quality),
				c("PhredQuality",
					"SolexaQuality",
					"IlluminaQuality")),
			alpha,
			threshold,
			as.integer(lefts),
			as.integer(rights),
			PACKAGE="DECIPHER")
		if (verbose) {
			wl <- which(lefts != bounds[[1]])
			wr <- which(rights != bounds[[2]])
			cat(100*round(length(wl)/length(lefts), 1),
				"% left, ",
				100*round(length(wr)/length(rights), 1),
				"% right",
				sep="")
			flush.console()
		}
		lefts <- bounds[[1]]
		rights <- bounds[[2]]
	}
	
	w <- which((rights - lefts + 1L) < minWidth)
	if (length(w) > 0) {
		lefts[w] <- 1L
		rights[w] <- 0L
	}
	
	if (type==1 || type==3)
		ranges <- IRanges(start=lefts,
			end=rights,
			names=names(myDNAStringSet))
	if (type==2 || type==3) {
		w <- which(lefts <= rights)
		if (length(w) > 0) {
			myDNAStringSet <- dna[w]
			myDNAStringSet <- subseq(myDNAStringSet,
				lefts[w],
				rights[w])
		} else {
			myDNAStringSet <- DNAStringSet()
		}
	}
	
	if (verbose) {
		cat("\n\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	if (type==1) {
		return(ranges)
	} else if (type==2) {
		return(myDNAStringSet)
	} else { # type==3
		return(list(ranges,
			myDNAStringSet))
	}
}
