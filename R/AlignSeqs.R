AlignSeqs <- function(myDNAStringSet,
	perfectMatch=2,
	misMatch=-3,
	gapOpening=-6,
	gapExtension=-4,
	terminalGap=-1,
	doNotAlign=0.75,
	guideTree=NULL,
	orient=TRUE,
	verbose=TRUE) {
	
	# error checking
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (min(width(myDNAStringSet)) < 2)
		stop("All sequences in myDNAStringSet must be at least two nucleotides long.")
	if (length(myDNAStringSet) < 2)
		stop("At least two sequences are required.")
	if (!is.null(guideTree)) {
		if (!is.data.frame(guideTree))
			stop("guideTree must be a data.frame.")
		if (dim(guideTree)[1] != length(myDNAStringSet))
			stop("guideTree have as many rows as sequences in myDNAStringSet.")
		if (!all(apply(guideTree, 2, is.numeric)))
			stop("guideTree must be a data.frame of numerics.")
	}
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!all(is.numeric(terminalGap)))
		stop("terminalGap must be a numeric.")
	if (length(terminalGap) > 2 || length(terminalGap) < 1)
		stop("terminalGap must be of length 1 or 2.")
	if (length(terminalGap)==1)
		terminalGap[2] <- terminalGap[1]
	if (!all(is.numeric(doNotAlign)))
		stop("doNotAlign must be a numeric.")
	if (length(doNotAlign)!=1)
		stop("doNotAlign must be length 1.")
	if (!is.logical(orient))
		stop("orient must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- alphabetFrequency(myDNAStringSet)
	if (any(a[, "-"]) > 0)
		stop("Gaps must be removed before alignment.")
	
	if (verbose)
		time.1 <- Sys.time()
	
	wordSize <- ceiling(mean(width(myDNAStringSet))^0.25)
	if (wordSize > 10)
		wordSize <- 10
	doNotAlign <- c(1 - (1 - doNotAlign)^wordSize*mean(width(myDNAStringSet)), doNotAlign)
	if (doNotAlign[1] < 0)
		doNotAlign[1] <- doNotAlign[2]
	
	if (orient) {
		tuples <- DNAStringSet(mkAllStrings(c("A", "C", "T", "G"), width=wordSize))
		pdict <- PDict(tuples,
			algorithm="ACtree2",
			tb.end=wordSize)
		v1 <- vwhichPDict(pdict,
			myDNAStringSet,
			fixed=TRUE)
		v1 <- unlist(lapply(v1, length))
		
		tuples <- reverseComplement(tuples)
		pdict <- PDict(tuples,
			algorithm="ACtree2",
			tb.end=wordSize)
		v2 <- vwhichPDict(pdict,
			myDNAStringSet,
			fixed=TRUE)
		v2 <- unlist(lapply(v2, length))
		
		w <- which(v2*0.5 > v1) # reverseComplement >> given orientation
		if (length(w) > 0)
			myDNAStringSet[w] <- reverseComplement(myDNAStringSet[w])
	}
	
	cluster <- FALSE
	cutoffs <- 0
	if (is.null(guideTree)) {
		cluster <- TRUE
		
		tuples <- mkAllStrings(c("A", "C", "T", "G"), width=wordSize)
		pdict <- PDict(tuples,
			algorithm="ACtree2",
			tb.end=wordSize)
		v <- vwhichPDict(pdict,
			myDNAStringSet,
			fixed=TRUE)
		
		if (verbose) {
			cat("Determining distance matrix based on shared ",
				wordSize,
				"-mers:\n",
				sep="")
			flush.console()
			pBar <- txtProgressBar(max=100, style=3)
		} else {
			pBar <- NULL
		}
		
		d <- .Call("matchLists", v, verbose, pBar, PACKAGE="DECIPHER")
		if (doNotAlign[1] < 1) {
			w <- which(is.na(d))
			if (length(w) > 0)
				d[w] <- doNotAlign[1] - 0.01
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
		
		if (verbose) {
			cat("\nClustering into groups by similarity:\n")
			flush.console()
		}
		cutoffs <- seq(0, max(d, na.rm=TRUE), 0.01)
		dimnames(d) <- NULL
		guideTree <- IdClusters(d, method="UPGMA", cutoff=cutoffs, verbose=verbose)
	}
	guideTree$Top <- 1
	cutoffs <- c(cutoffs, doNotAlign[1])
	
	# calculate weights based on branch lengths
	lastSplit <- numeric(dim(guideTree)[1])
	lastSplit[] <- 1
	weights <- numeric(dim(guideTree)[1])
	numGroups <- 1
	for (i in (dim(guideTree)[2] - 1):1) {
		if (length(unique(guideTree[, i])) > numGroups) {
			numGroups <- length(unique(guideTree[, i]))
			for (j in unique(guideTree[, i + 1])) {
				w <- which(guideTree[, i + 1]==j)
				if (length(unique(guideTree[w, i])) > 1) {
					weights[w] <- weights[w] + (lastSplit[w] - cutoffs[i])/length(w)
					lastSplit[w] <- cutoffs[i]
				}
			}
		}
	}
	weights <- weights + lastSplit
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Aligning Sequences:\n")
		flush.console()
		pBar <- txtProgressBar(style=3)
	}
	
	# determine the number of steps
	dna <- myDNAStringSet
	ns <- names(myDNAStringSet)
	steps <- 0
	for (i in 1:dim(guideTree)[2]) {
		if (cutoffs[i] > doNotAlign[1])
			break
		for (u in unique(guideTree[, i])) {
			w <- which(guideTree[, i]==u) # groups
			if (length(w)==1)
				next
			if (i==1) {
				for (j in 2:length(w)) {
					if (length(which(duplicated(dna[w[1:j]]))) == (j - 1))
						next
					steps <- steps + 1
				}
			} else {
				g <- unique(guideTree[w, i - 1]) # sub-groups
				steps <- steps + length(g) - 1
			}
		}
	}
	l <- steps
	steps <- 0
	mergers <- character(l)
	
	for (i in 1:dim(guideTree)[2]) {
		if (cutoffs[i] > doNotAlign[1]) {
			dna <- xscat(dna,
				unlist(lapply(max(width(dna)) - width(dna),
					function(x) paste(rep("-", times=x), collapse=""))))
			break
		}
		for (u in unique(guideTree[, i])) {
			w <- which(guideTree[, i]==u) # groups
			if (length(w)==1)
				next
			if (i==1) {
				for (j in 2:length(w)) {
					if (length(which(duplicated(dna[w[1:j]]))) == (j - 1))
						next
					
					p.weight <- weights[w[1:(j - 1)]]
					p.weight <- p.weight/mean(p.weight)
					s.weight <- weights[w[j]]
					s.weight <- s.weight/mean(s.weight)
					
					dna[w[1:j]] <- AlignProfiles(dna[w[1:(j - 1)]],
						dna[w[j]],
						p.weight,
						s.weight,
						perfectMatch,
						misMatch,
						gapOpening,
						gapExtension,
						terminalGap)
					
					steps <- steps + 1
					x <- grep(w[1], mergers[1:steps], fixed=TRUE)
					if (length(x) > 0) {
						mergers[steps] <- paste(x[length(x)],
							w[j],
							sep="|")
					} else {
						mergers[steps] <- paste(paste(w[1:(j - 1)], collapse=","),
							w[j],
							sep="|")
					}
					
					if (verbose)
						setTxtProgressBar(pBar, steps/l)
				}
			} else {
				g <- unique(guideTree[w, i - 1]) # sub-groups
				if (length(g)==1)
					next
				g1 <- which(guideTree[w, i - 1]==g[1])
				for (j in 2:length(g)) {
					g2 <- which(guideTree[w, i - 1]==g[j])
					
					p.weight <- weights[w[g1]]
					p.weight <- p.weight/mean(p.weight)
					s.weight <- weights[w[g2]]
					s.weight <- s.weight/mean(s.weight)
					
					dna[w[c(g1, g2)]] <- AlignProfiles(dna[w[g1]],
						dna[w[g2]],
						p.weight,
						s.weight,
						perfectMatch,
						misMatch,
						gapOpening,
						gapExtension,
						terminalGap)
					
					steps <- steps + 1
					x <- grep(w[g1[1]], mergers[1:steps], fixed=TRUE)
					y <- grep(w[g2[1]], mergers[1:steps], fixed=TRUE)
					if (!(length(x) > 0 && nchar(mergers[x[length(x)]]) > 1000) &&
						!(length(y) > 0 && nchar(mergers[y[length(y)]]) > 1000)) {
						if (length(x) > 0 && length(y) > 0) {
							mergers[steps] <- paste(mergers[x[length(x)]],
								mergers[y[length(y)]],
								sep="|")
						} else if (length(x) > 0) {
							mergers[steps] <- paste(mergers[x[length(x)]],
								paste(w[g2], collapse=","),
								sep="|")
						} else if (length(y) > 0) {
							mergers[steps] <- paste(paste(w[g1], collapse=","),
								mergers[y[length(y)]],
								sep="|")
						} else {
							mergers[steps] <- paste(paste(w[g1], collapse=","),
								paste(w[g2], collapse=","),
								sep="|")
						}
					}
					
					if (verbose)
						setTxtProgressBar(pBar, steps/l)
					g1 <- c(g1, g2)
				}
			}
		}
	}
	names(dna) <- ns
	
	if (cluster) {
		if (verbose) {
			setTxtProgressBar(pBar, 1)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\nDetermining distance matrix based on alignment:\n")
			flush.console()
		}
		
		suppressWarnings(d <- DistanceMatrix(dna, verbose=verbose))
		if (doNotAlign[2] < 1) {
			w <- which(is.na(d))
			if (length(w) > 0)
				d[w] <- doNotAlign[2] - 0.01
		}
		
		if (verbose) {
			cat("Reclustering into groups by similarity:\n")
			flush.console()
		}
		orgTree <- guideTree
		cutoffs <- seq(0, max(d, na.rm=TRUE), length.out=100)
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- IdClusters(d, cutoff=cutoffs, verbose=verbose, method="UPGMA"))
		guideTree$Top <- 1
		cutoffs <- c(cutoffs, doNotAlign[2])
		w <- which(is.na(d))
		if (length(w) > 0)
			d[w] <- 1
		
		# calculate weights based on branch lengths
		lastSplit <- numeric(dim(guideTree)[1])
		lastSplit[] <- 1
		weights <- numeric(dim(guideTree)[1])
		numGroups <- 1
		for (i in (dim(guideTree)[2] - 1):1) {
			if (length(unique(guideTree[, i])) > numGroups) {
				numGroups <- length(unique(guideTree[, i]))
				for (j in unique(guideTree[, i + 1])) {
					w <- which(guideTree[, i + 1]==j)
					if (length(unique(guideTree[w, i])) > 1) {
						weights[w] <- weights[w] + (lastSplit[w] - cutoffs[i])/length(w)
						lastSplit[w] <- cutoffs[i]
					}
				}
			}
		}
		weights <- weights + lastSplit
		
		if (verbose) {
			time.1 <- Sys.time()
			cat("Realigning Sequences:\n")
			flush.console()
			pBar <- txtProgressBar(style=3)
		}
		
		# determine the number of steps
		steps <- 0
		for (i in 1:dim(guideTree)[2]) {
			if (cutoffs[i] > doNotAlign[2])
				break
			for (u in unique(guideTree[, i])) {
				w <- which(guideTree[, i]==u) # groups
				if (length(w)==1)
					next
				if (i==1) {
					for (j in 2:length(w)) {
						if (length(which(duplicated(myDNAStringSet[w[1:j]]))) == (j - 1))
							next
						steps <- steps + 1
					}
				} else {
					g <- unique(guideTree[w, i - 1]) # sub-groups
					steps <- steps + length(g) - 1
				}
			}
		}
		l <- steps
		steps <- 0
		mergers2 <- character(l)
		
		ns <- names(myDNAStringSet)
		for (i in 1:dim(guideTree)[2]) {
			if (cutoffs[i] > doNotAlign[2]) {
				myDNAStringSet <- xscat(myDNAStringSet,
					unlist(lapply(max(width(myDNAStringSet)) - width(myDNAStringSet),
						function(x) paste(rep("-", times=x), collapse=""))))
				break
			}
			for (u in unique(guideTree[, i])) {
				w <- which(guideTree[, i]==u) # groups
				if (length(w)==1)
					next
				
				if (i==1) {
					for (j in 2:length(w)) {
						if (length(which(duplicated(myDNAStringSet[w[1:j]]))) == (j - 1))
							next
						
						steps <- steps + 1
						x <- grep(w[1], mergers2[1:steps], fixed=TRUE)
						if (length(x) > 0) {
							mergers2[steps] <- paste(x[length(x)],
								w[j],
								sep="|")
						} else {
							mergers2[steps] <- paste(paste(w[1:(j - 1)], collapse=","),
								w[j],
								sep="|")
						}
						if (length(grep(mergers2[steps], mergers, fixed=TRUE) > 0)) {
							myDNAStringSet[w[1:j]] <- DNAStringSet(.Call("commonGaps",
								substr(dna[w[1:j]], 1, width(dna[w[1]])),
								PACKAGE="DECIPHER"))
							next
						}
						
						p.weight <- weights[w[1:(j - 1)]]
						p.weight <- p.weight/mean(p.weight)
						s.weight <- weights[w[j]]
						s.weight <- s.weight/mean(s.weight)
						
						myDNAStringSet[w[1:j]] <- AlignProfiles(myDNAStringSet[w[1:(j - 1)]],
							myDNAStringSet[w[j]],
							p.weight,
							s.weight,
							perfectMatch,
							misMatch,
							gapOpening,
							gapExtension,
							terminalGap)
						
						if (verbose)
							setTxtProgressBar(pBar, steps/l)
					}
				} else {
					g <- unique(guideTree[w, i - 1]) # sub-groups
					if (length(g)==1)
						next
					g1 <- which(guideTree[w, i - 1]==g[1])
					for (j in 2:length(g)) {
						g2 <- which(guideTree[w, i - 1]==g[j])
						
						steps <- steps + 1
						x <- grep(w[g1[1]], mergers2[1:steps], fixed=TRUE)
						y <- grep(w[g2[1]], mergers2[1:steps], fixed=TRUE)
						if (!(length(x) > 0 && nchar(mergers2[x[length(x)]]) > 1000) &&
							!(length(y) > 0 && nchar(mergers2[y[length(y)]]) > 1000)) {
							if (length(x) > 0 && length(y) > 0) {
								mergers2[steps] <- paste(mergers2[x[length(x)]],
									mergers2[y[length(y)]],
									sep="|")
							} else if (length(x) > 0) {
								mergers2[steps] <- paste(mergers2[x[length(x)]],
									paste(w[g2], collapse=","),
									sep="|")
							} else if (length(y) > 0) {
								mergers2[steps] <- paste(paste(w[g1], collapse=","),
									mergers2[y[length(y)]],
									sep="|")
							} else {
								mergers2[steps] <- paste(paste(w[g1], collapse=","),
									paste(w[g2], collapse=","),
									sep="|")
							}
							if (length(grep(mergers2[steps], mergers, fixed=TRUE) > 0)) {
								myDNAStringSet[c(w[g1], w[g2])] <- DNAStringSet(.Call("commonGaps",
									substr(dna[c(w[g1], w[g2])], 1, width(dna[w[g1[1]]])),
									PACKAGE="DECIPHER"))
								if (verbose)
									setTxtProgressBar(pBar, steps/l)
								g1 <- c(g1, g2)
								next
							}
						}
						
						p.weight <- weights[w[g1]]
						p.weight <- p.weight/mean(p.weight)
						s.weight <- weights[w[g2]]
						s.weight <- s.weight/mean(s.weight)
						
						myDNAStringSet[w[c(g1, g2)]] <- AlignProfiles(myDNAStringSet[w[g1]],
								myDNAStringSet[w[g2]],
								p.weight,
								s.weight,
								perfectMatch,
								misMatch,
								gapOpening,
								gapExtension,
								terminalGap)
						
						if (verbose)
							setTxtProgressBar(pBar, steps/l)
						g1 <- c(g1, g2)
					}
				}
			}
		}
		names(myDNAStringSet) <- ns
	} else {
		myDNAStringSet <- dna
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(myDNAStringSet)
}
