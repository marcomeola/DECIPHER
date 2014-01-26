AlignSeqs <- function(myXStringSet,
	guideTree=NULL,
	orient=FALSE,
	processors=NULL,
	verbose=TRUE,
	...) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet") && !is(myXStringSet, "AAStringSet"))
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	if (length(myXStringSet) < 2)
		stop("At least two sequences are required.")
	if (!is.null(guideTree)) {
		if (!is.data.frame(guideTree))
			stop("guideTree must be a data.frame.")
		if (dim(guideTree)[1] != length(myXStringSet))
			stop("guideTree must have as many rows as sequences in myXStringSet.")
		if (!all(apply(guideTree, 2, is.numeric)))
			stop("guideTree must be a data.frame of numerics.")
	}
	if (!is.logical(orient))
		stop("orient must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a) > 0)
		stop("Gaps must be removed before alignment.")
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
	
	if (verbose)
		time.1 <- Sys.time()
	
	if (is(myXStringSet, "AAStringSet")) {
		wordSize <- ceiling(mean(width(myXStringSet))^0.2) # always >= 1
		if (wordSize > 7)
			wordSize <- 7
	} else {
		wordSize <- ceiling(mean(width(myXStringSet))^0.25) # always >= 1
		if (wordSize > 15)
			wordSize <- 15
	}
	
	if (orient && !is(myXStringSet, "AAStringSet")) {
		w <- which.max(width(myXStringSet))
		v1 <- .Call("enumerateSequence",
			myXStringSet,
			wordSize,
			PACKAGE="DECIPHER")
		X <- v1[[w]]
		
		f <- function(x) {
			m <- match(x, X)
			return(length(which(!is.na(m))))
		}
		v1 <- unlist(lapply(v1, f))
		
		v2 <- .Call("enumerateSequence",
			reverseComplement(myXStringSet),
			wordSize,
			PACKAGE="DECIPHER")
		v2 <- unlist(lapply(v2, f))
		
		v3 <- .Call("enumerateSequence",
			reverse(myXStringSet),
			wordSize,
			PACKAGE="DECIPHER")
		v3 <- unlist(lapply(v3, f))
		
		v4 <- .Call("enumerateSequence",
			complement(myXStringSet),
			wordSize,
			PACKAGE="DECIPHER")
		v4 <- unlist(lapply(v4, f))
		
		w <- which(v2*0.5 > v1 & v2 > v3 & v2 > v4)
		if (length(w) > 0) # reverseComplement >> given orientation
			myXStringSet[w] <- reverseComplement(myXStringSet[w])
		w <- which(v3*0.5 > v1 & v3 > v2 & v3 > v4)
		if (length(w) > 0) # reverse >> given orientation
			myXStringSet[w] <- reverse(myXStringSet[w])
		w <- which(v4*0.5 > v1 & v4 > v2 & v4 > v3)
		if (length(w) > 0) # complement >> given orientation
			myXStringSet[w] <- complement(myXStringSet[w])
	} else if (orient && is(myXStringSet, "AAStringSet")) {
		warning("orient is ignored when myXStringSet is an AAStringSet.")
	}
	
	cluster <- FALSE
	cutoffs <- 0
	if (is.null(guideTree)) {
		cluster <- TRUE
		
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
		
		if (is(myXStringSet, "AAStringSet")) {
			v <- .Call("enumerateSequenceAA",
				myXStringSet,
				wordSize,
				PACKAGE="DECIPHER")
		} else { # DNAStringSet or RNAStringSet
			v <- .Call("enumerateSequence",
				myXStringSet,
				wordSize,
				PACKAGE="DECIPHER")
		}
		
		d <- .Call("matchOrder", v, verbose, pBar, processors, PACKAGE="DECIPHER")
		
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
		guideTree <- suppressWarnings(IdClusters(d, method="UPGMA", cutoff=cutoffs, verbose=verbose, processors=processors))
	} else {
		cutoffs <- rep(1, dim(guideTree)[2])
	}
	guideTree$Top <- 1
	cutoffs <- c(cutoffs, 1)
	
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
	seqs <- myXStringSet
	ns <- names(myXStringSet)
	steps <- 0
	for (i in 1:dim(guideTree)[2]) {
		for (u in unique(guideTree[, i])) {
			w <- which(guideTree[, i]==u) # groups
			if (length(w)==1)
				next
			if (i==1) {
				for (j in 2:length(w)) {
					if (length(which(duplicated(seqs[w[1:j]]))) == (j - 1))
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
		subMatrix <- ifelse(cutoffs[i] < 0.01,
			"BLOSUM100",
				ifelse(cutoffs[i] < 0.2,
					"BLOSUM80",
					ifelse(cutoffs[i] < 0.4,
						"BLOSUM62",
						ifelse(cutoffs[i] < 0.5,
							"BLOSUM50",
							"BLOSUM45"))))
		for (u in unique(guideTree[, i])) {
			w <- which(guideTree[, i]==u) # groups
			if (length(w)==1)
				next
			if (i==1) {
				for (j in 2:length(w)) {
					if (length(which(duplicated(seqs[w[1:j]]))) == (j - 1))
						next
					
					p.weight <- weights[w[1:(j - 1)]]
					p.weight <- p.weight/mean(p.weight)
					s.weight <- weights[w[j]]
					s.weight <- s.weight/mean(s.weight)
					
					seqs[w[1:j]] <- AlignProfiles(pattern=seqs[w[1:(j - 1)]],
						subject=seqs[w[j]],
						p.weight=p.weight,
						s.weight=s.weight,
						substitutionMatrix="BLOSUM100",
						processors=processors,
						...)
					
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
					
					seqs[w[c(g1, g2)]] <- AlignProfiles(pattern=seqs[w[g1]],
						subject=seqs[w[g2]],
						p.weight=p.weight,
						s.weight=s.weight,
						substitutionMatrix=subMatrix,
						processors=processors,
						...)
					
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
	names(seqs) <- ns
	
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
		
		suppressWarnings(d <- DistanceMatrix(seqs, verbose=verbose, processors=processors))
		
		if (verbose) {
			cat("Reclustering into groups by similarity:\n")
			flush.console()
		}
		orgTree <- guideTree
		cutoffs <- seq(0, max(d, na.rm=TRUE), length.out=100)
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- suppressWarnings(IdClusters(d, cutoff=cutoffs, verbose=verbose, method="UPGMA", processors=processors)))
		guideTree$Top <- 1
		cutoffs <- c(cutoffs, 1)
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
			for (u in unique(guideTree[, i])) {
				w <- which(guideTree[, i]==u) # groups
				if (length(w)==1)
					next
				if (i==1) {
					for (j in 2:length(w)) {
						if (length(which(duplicated(myXStringSet[w[1:j]]))) == (j - 1))
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
		
		ns <- names(myXStringSet)
		for (i in 1:dim(guideTree)[2]) {
			subMatrix <- ifelse(cutoffs[i] < 0.01,
				"BLOSUM100",
					ifelse(cutoffs[i] < 0.2,
						"BLOSUM80",
						ifelse(cutoffs[i] < 0.4,
							"BLOSUM62",
							ifelse(cutoffs[i] < 0.5,
								"BLOSUM50",
								"BLOSUM45"))))
			for (u in unique(guideTree[, i])) {
				w <- which(guideTree[, i]==u) # groups
				if (length(w)==1)
					next
				
				if (i==1) {
					for (j in 2:length(w)) {
						if (length(which(duplicated(myXStringSet[w[1:j]]))) == (j - 1))
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
							if (is(myXStringSet, "AAStringSet")) {
								myXStringSet[w[1:j]] <- AAStringSet(.Call("commonGaps",
									substr(seqs[w[1:j]], 1, width(seqs[w[1]])),
									PACKAGE="DECIPHER"))
							} else if (is(myXStringSet, "DNAStringSet")) {
								myXStringSet[w[1:j]] <- DNAStringSet(.Call("commonGaps",
									substr(seqs[w[1:j]], 1, width(seqs[w[1]])),
									PACKAGE="DECIPHER"))
							} else { # RNAStringSet
								myXStringSet[w[1:j]] <- RNAStringSet(.Call("commonGaps",
									substr(seqs[w[1:j]], 1, width(seqs[w[1]])),
									PACKAGE="DECIPHER"))
							}
							next
						}
						
						p.weight <- weights[w[1:(j - 1)]]
						p.weight <- p.weight/mean(p.weight)
						s.weight <- weights[w[j]]
						s.weight <- s.weight/mean(s.weight)
						
						myXStringSet[w[1:j]] <- AlignProfiles(pattern=myXStringSet[w[1:(j - 1)]],
							subject=myXStringSet[w[j]],
							p.weight=p.weight,
							s.weight=s.weight,
							substitutionMatrix="BLOSUM100",
							processors=processors,
							...)
						
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
								if (is(myXStringSet, "AAStringSet")) {
									myXStringSet[c(w[g1], w[g2])] <- AAStringSet(.Call("commonGaps",
										substr(seqs[c(w[g1], w[g2])], 1, width(seqs[w[g1[1]]])),
										PACKAGE="DECIPHER"))
								} else if (is(myXStringSet, "DNAStringSet")) {
									myXStringSet[c(w[g1], w[g2])] <- DNAStringSet(.Call("commonGaps",
										substr(seqs[c(w[g1], w[g2])], 1, width(seqs[w[g1[1]]])),
										PACKAGE="DECIPHER"))
								} else { # RNAStringSet
									myXStringSet[c(w[g1], w[g2])] <- RNAStringSet(.Call("commonGaps",
										substr(seqs[c(w[g1], w[g2])], 1, width(seqs[w[g1[1]]])),
										PACKAGE="DECIPHER"))
								}
								
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
						
						myXStringSet[w[c(g1, g2)]] <- AlignProfiles(pattern=myXStringSet[w[g1]],
								subject=myXStringSet[w[g2]],
								p.weight=p.weight,
								s.weight=s.weight,
								substitutionMatrix=subMatrix,
								processors=processors,
								...)
						
						if (verbose)
							setTxtProgressBar(pBar, steps/l)
						g1 <- c(g1, g2)
					}
				}
			}
		}
		names(myXStringSet) <- ns
	} else {
		myXStringSet <- seqs
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
	
	return(myXStringSet)
}
