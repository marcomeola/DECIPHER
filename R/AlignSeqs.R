AlignSeqs <- function(myXStringSet,
	guideTree=NULL,
	iterations=1,
	refinements=1,
	gapOpening=c(-16, -12),
	gapExtension=c(-2, -1),
	structures=NULL,
	FUN=AdjustAlignment,
	levels=c(0.95, 0.7, 10, 5),
	processors=NULL,
	verbose=TRUE,
	...) {
	
	# error checking
	type <- switch(class(myXStringSet),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet."))
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
	if (!is.numeric(iterations))
		stop("iterations must be a numeric.")
	if (iterations != floor(iterations))
		stop("iterations must be a whole number.")
	if (iterations < 0)
		stop("iterations must be at least zero.")
	if (!is.numeric(refinements))
		stop("refinements must be a numeric.")
	if (refinements != floor(refinements))
		stop("refinements must be a whole number.")
	if (refinements < 0)
		stop("refinements must be at least zero.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') must be removed before alignment.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') must be removed before alignment.")
	FUN <- match.fun(FUN)
	if (!is.numeric(levels))
		stop("levels must be a numeric vector.")
	if (length(levels) != 4)
		stop("levels must be length 3.")
	if (levels[1] < 0 || levels[1] > 1)
		stop("levels[1] must be between zero and one (inclusive).")
	if (levels[2] < 0 || levels[2] > 1)
		stop("levels[2] must be between zero and one (inclusive).")
	if (levels[3] <= 0 || is.infinite(levels[3]))
		stop("levels[3] must be between zero and Inf (exclusive).")
	if (levels[3] != floor(levels[3]))
		stop("levels[3] must be a whole number.")
	if (levels[4] != floor(levels[4]))
		stop("levels[4] must be a whole number.")
	if (levels[4] < 0)
		stop("levels[4] must be at least zero.")
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
	
	args <- list(...)
	n <- names(args)
	m <- character(length(n))
	for (i in seq_along(n)) {
		m[i] <- match.arg(n[i],
			names(formals(AlignProfiles)))
	}
	
	if (length(gapOpening)==2) {
		gapOpeningMin <- gapOpening[1]
		gapOpeningMax <- gapOpening[2]
	} else if (length(gapOpening)==1) {
		gapOpeningMin <- gapOpeningMax <- gapOpening
	} else {
		stop("gapOpening must be a vector of length 1 or 2.")
	}
	if (length(gapExtension)==2) {
		gapExtensionMin <- gapExtension[1]
		gapExtensionMax <- gapExtension[2]
	} else if (length(gapExtension)==1) {
		gapExtensionMin <- gapExtensionMax <- gapExtension
	} else {
		stop("gapExtension must be a vector of length 1 or 2.")
	}
	if (!is.numeric(gapOpeningMax))
		stop("gapOpening must be a numeric.")
	if (gapOpeningMin > gapOpeningMax)
		stop("gapOpening[1] must be less than or equal to gapOpening[2].")
	if (!is.numeric(gapExtensionMax))
		stop("gapExtension must be a numeric.")
	if (gapExtensionMin > gapExtensionMax)
		stop("gapExtension[1] must be less than or equal to gapExtension[2].")
	gapOpeningSlope <- gapOpeningMax - gapOpeningMin
	gapExtensionSlope <- gapExtensionMax - gapExtensionMin
	
	if (type==3L) {
		# structures are always used with AA sequences
		if (is.null(structures)) {
			structures <- PredictHEC(myXStringSet,
				type="probabilities")
		} else {
			if (length(structures) != length(myXStringSet))
				stop("structures is not the same length as myXStringSet.")
			if (typeof(structures) != "list")
				stop("structures must be a list.")
		}
		if (refinements > 0) {
			w <- which(m=="structureMatrix")
			if (length(w) > 0) {
				structureMatrix <- args[[w]]
				# assume structures and matrix are ordered the same
				if (!is.double(structureMatrix))
					stop("structureMatrix must be contain numerics.")
				if (!is.matrix(structureMatrix))
					stop("structureMatrix must be a matrix.")
				if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
					stop("structureMatrix is not square.")
			} else {
				structureMatrix <- matrix(c(5, 0, -2, 0, 9, -1, -2, -1, 2),
					nrow=3) # order is H, E, C
			}
			if (dim(structureMatrix)[1] != dim(structures[[1]])[1])
				stop("Dimensions of structureMatrix are incompatible with the structures.")
		}
		
		subM <- TRUE
		w <- which(m=="substitutionMatrix")
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (length(w)==1L) {
			sM <- args[[w]]
			args[w] <- NULL
			if (is.character(sM)) {
				if (!(sM %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
			"PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")))
					stop("Invalid substitutionMatrix.")
			}
			AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
				"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
			if (is.matrix(sM)) {
				if (any(!(AAs %in% dimnames(sM)[[1]])) ||
					any(!(AAs %in% dimnames(sM)[[2]])))
					stop("substitutionMatrix is incomplete.")
			} else {
				sM <- eval(parse(text=data(list=sM, envir=environment())))
			}
			sM <- sM[AAs, AAs]
			sM <- sM + 0 # convert to numeric matrix
		} else {
			# use MIQS
			sM <- matrix(c(3.2,-1.3,-0.4,-0.4,1.5,-0.2,-0.4,0.4,-1.2,-1.3,-1.4,-0.7,-1,-2.3,-0.1,0.8,0.8,-3.6,-2.4,0,-6.1,-1.3,6.2,-0.1,-1.5,-2.7,1.8,-0.7,-1.9,0.9,-2.4,-2.5,3.3,-1.1,-3.3,-1.1,-0.3,-0.9,-3.8,-1.9,-2.3,-6.1,-0.4,-0.1,5.1,2.6,-1.6,0.9,0.8,0.2,1,-3.6,-3.5,0.7,-2.3,-3.5,-1.4,0.9,0,-4.5,-1.5,-2.6,-6.1,-0.4,-1.5,2.6,5.7,-3.7,0.9,2.7,-0.5,0.3,-4.5,-4.6,0.4,-3.3,-5.8,-0.3,0.3,-0.2,-5.3,-3.9,-3.5,-6.1,1.5,-2.7,-1.6,-3.7,11.7,-2.8,-3.2,-1.7,-1.2,0.2,-2.3,-3.2,0.1,-2.8,-2.8,1,0,-6.1,-0.7,1.8,-6.1,-0.2,1.8,0.9,0.9,-2.8,3.6,2.1,-1.6,1.2,-2.2,-1.9,1.7,-0.4,-2.4,-0.4,0.4,0.1,-5.4,-2.8,-1.8,-6.1,-0.4,-0.7,0.8,2.7,-3.2,2.1,4.3,-1.3,-0.2,-3.3,-2.8,1.1,-2.3,-4.1,0,0.4,-0.2,-5.8,-2.4,-2.3,-6.1,0.4,-1.9,0.2,-0.5,-1.7,-1.6,-1.3,7.6,-1.6,-5.4,-4.8,-1.7,-3.6,-4.6,-1.6,0,-1.9,-4.8,-4.5,-3.8,-6.1,-1.2,0.9,1,0.3,-1.2,1.2,-0.2,-1.6,7.5,-2.2,-1.9,0,-2.1,0,-1.5,0,-0.2,-0.3,2.1,-2.3,-6.1,-1.3,-2.4,-3.6,-4.5,0.2,-2.2,-3.3,-5.4,-2.2,4.6,3.1,-2.3,1.7,0.7,-3.7,-2.8,-0.7,-0.7,-0.8,3.3,-6.1,-1.4,-2.5,-3.5,-4.6,-2.3,-1.9,-2.8,-4.8,-1.9,3.1,4.6,-2.4,3.2,2.1,-2.8,-2.9,-1.6,-0.2,0,2,-6.1,-0.7,3.3,0.7,0.4,-3.2,1.7,1.1,-1.7,0,-2.3,-2.4,3.6,-1.1,-3.7,-0.1,0,0,-4,-2.3,-2,-6.1,-1,-1.1,-2.3,-3.3,0.1,-0.4,-2.3,-3.6,-2.1,1.7,3.2,-1.1,5.4,1.4,-2.8,-1.8,-0.8,-2.1,-0.9,1.4,-6.1,-2.3,-3.3,-3.5,-5.8,-2.8,-2.4,-4.1,-4.6,0,0.7,2.1,-3.7,1.4,7.4,-3.7,-2.6,-2.3,4.2,5.2,-0.3,-6.1,-0.1,-1.1,-1.4,-0.3,-2.8,-0.4,0,-1.6,-1.5,-3.7,-2.8,-0.1,-2.8,-3.7,8.4,-0.1,-0.5,-3.6,-4.5,-2.5,-6.1,0.8,-0.3,0.9,0.3,1,0.4,0.4,0,0,-2.8,-2.9,0,-1.8,-2.6,-0.1,3.1,1.6,-3.5,-1.5,-1.4,-6.1,0.8,-0.9,0,-0.2,0,0.1,-0.2,-1.9,-0.2,-0.7,-1.6,0,-0.8,-2.3,-0.5,1.6,3.8,-5.3,-2.1,-0.1,-6.1,-3.6,-3.8,-4.5,-5.3,-6.1,-5.4,-5.8,-4.8,-0.3,-0.7,-0.2,-4,-2.1,4.2,-3.6,-3.5,-5.3,14.8,4.9,-3.3,-6.1,-2.4,-1.9,-1.5,-3.9,-0.7,-2.8,-2.4,-4.5,2.1,-0.8,0,-2.3,-0.9,5.2,-4.5,-1.5,-2.1,4.9,8.3,-1.2,-6.1,0,-2.3,-2.6,-3.5,1.8,-1.8,-2.3,-3.8,-2.3,3.3,2,-2,1.4,-0.3,-2.5,-1.4,-0.1,-3.3,-1.2,3.5,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,1),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		}
	} else {
		if (is.null(structures)) {
			subM <- FALSE
			w <- which(m=="structureMatrix")
			if (length(w) > 0)
				stop("structureMatrix provided without structures.")
			structureMatrix <- numeric() # needed for colScores
		} else {
			subM <- TRUE
			if (length(structures) != length(myXStringSet))
				stop("structures is not the same length as myXStringSet.")
			if (typeof(structures) != "list")
				stop("structures must be a list.")
			
			w <- which(m=="structureMatrix")
			if (length(w) > 0) {
				structureMatrix <- args[[w]]
				# assume structures and matrix are ordered the same
				if (!is.double(structureMatrix))
					stop("structureMatrix must be contain numerics.")
				if (!is.matrix(structureMatrix))
					stop("structureMatrix must be a matrix.")
				if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
					stop("structureMatrix is not square.")
			} else {
				stop("structureMatrix must be specified when structures are provided.")
			}
			if (dim(structureMatrix)[1] != dim(structures[[1]])[1])
				stop("Dimensions of structureMatrix are incompatible with the structures.")
		}
		
		bases <- c("A", "C", "G",
			ifelse(type==2L, "U", "T"))
		w <- which(m=="substitutionMatrix")
		if (length(w) > 0) {
			sM <- args[[w]]
			args[w] <- NULL
			if (is.matrix(sM)) {
				if (any(!(bases %in% dimnames(sM)[[1]])) ||
					any(!(bases %in% dimnames(sM)[[2]])))
					stop("substitutionMatrix is incomplete.")
				sM <- sM[bases, bases]
				sM <- sM + 0 # convert to numeric matrix
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		} else {
			w <- which(m=="perfectMatch")
			if (length(w) > 0) {
				PM <- args[[w]]
				if (!is.numeric(PM))
					stop("perfectMatch must be a numeric.")
			} else {
				PM <- 5
			}
			w <- which(m=="misMatch")
			if (length(w) > 0) {
				MM <- args[[w]]
				if (!is.numeric(MM))
					stop("misMatch must be a numeric.")
			} else {
				MM <- 0
			}
			sM <- matrix(as.numeric(MM),
				nrow=4,
				ncol=4,
				dimnames=list(bases, bases))
			diag(sM) <- PM
		}
	}
	if (length(args)==0)
		args <- NULL
	
	if (verbose)
		time.1 <- Sys.time()
	
	if (type==3L) {
		wordSize <- ceiling(mean(width(myXStringSet))^0.2)
		if (wordSize > 7)
			wordSize <- 7
	} else {
		wordSize <- ceiling(mean(width(myXStringSet))^0.25)
		if (wordSize > 15)
			wordSize <- 15
		if (wordSize < 8)
			wordSize <- 8
	}
	
	cluster <- FALSE
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
		
		if (type==3L) {
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
		cutoffs <- sort(unique(c(seq(0, 1, 0.01),
			round(quantile(d, seq(0, 1, 0.01), na.rm=TRUE), 5))))
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- IdClusters(d,
			method="single",
			cutoff=cutoffs,
			verbose=verbose,
			processors=processors))
	} else {
		cutoffs <- rep(1, dim(guideTree)[2])
	}
	guideTree$Top <- 1
	guideTree <- as.matrix(guideTree)
	guideTree[, 1] <- 1:dim(guideTree)[1]
	cutoffs <- c(cutoffs, 1)
	
	# calculate weights based on branch lengths
	lastSplit <- numeric(dim(guideTree)[1])
	weights <- numeric(dim(guideTree)[1])
	numGroups <- 1
	for (i in (dim(guideTree)[2] - 1):1) {
		if (length(unique(guideTree[, i])) > numGroups) {
			if (numGroups==1) { # mid-point root of tree
				lastSplit[] <- cutoffs[i]
			} else {
				for (j in unique(guideTree[, i + 1])) {
					w <- which(guideTree[, i + 1]==j)
					if (length(unique(guideTree[w, i])) > 1) {
						weights[w] <- weights[w] + (lastSplit[w] - cutoffs[i])/length(w)
						lastSplit[w] <- cutoffs[i]
					}
				}
			}
			numGroups <- length(unique(guideTree[, i]))
		}
	}
	weights <- weights + lastSplit
	weights <- ifelse(weights==0, 1, weights)
	weights <- weights/mean(weights)
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Aligning Sequences:\n")
		flush.console()
		pBar <- txtProgressBar(style=3, max=100)
		percentComplete <- before <- 0L
	}
	
	# determine the number of steps
	seqs <- seqs_original <- myXStringSet
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
	
	GO <- gapOpeningMin + gapOpeningSlope*cutoffs
	GE <- gapExtensionMin + gapExtensionSlope*cutoffs
	
	for (i in 1:dim(guideTree)[2]) {
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
					
					if (subM) {
						pattern <- .Call("subsetXStringSet",
							seqs,
							w[1:(j - 1)],
							type,
							processors)
						subject <- .Call("subsetXStringSet",
							seqs,
							w[j],
							type,
							processors)
						seqs[w[1:j]] <- do.call(AlignProfiles,
							args=c(list(pattern=pattern,
									subject=subject,
									p.weight=p.weight,
									s.weight=s.weight,
									p.struct=structures[w[1:(j - 1)]],
									s.struct=structures[w[j]],
									substitutionMatrix=sM,
									processors=processors,
									gapOpening=GO[i],
									gapExtension=GE[i]),
								args))
					} else { # do not need to specify substitutionMatrix
						pattern <- .Call("subsetXStringSet",
							seqs,
							w[1:(j - 1)],
							type,
							processors)
						subject <- .Call("subsetXStringSet",
							seqs,
							w[j],
							type,
							processors)
						seqs[w[1:j]] <- AlignProfiles(pattern=pattern,
							subject=subject,
							p.weight=p.weight,
							s.weight=s.weight,
							processors=processors,
							gapOpening=GO[i],
							gapExtension=GE[i],
							...)
					}
					
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
					
					if (verbose) {
						percentComplete <- as.integer(100*steps/l)
						if (percentComplete > before) {
							setTxtProgressBar(pBar, percentComplete)
							before <- percentComplete
						}
					}
				}
			} else {
				g <- unique(guideTree[w, i - 1]) # sub-groups
				if (length(g)==1)
					next
				g <- g[order(width(myXStringSet)[match(g, guideTree[w, i - 1])],
					decreasing=TRUE)] # sort sub-groups by increasing width
				g1 <- which(guideTree[w, i - 1]==g[1])
				for (j in 2:length(g)) {
					g2 <- which(guideTree[w, i - 1]==g[j])
					
					p.weight <- weights[w[g1]]
					p.weight <- p.weight/mean(p.weight)
					s.weight <- weights[w[g2]]
					s.weight <- s.weight/mean(s.weight)
					
					if (subM) {
						pattern <- .Call("subsetXStringSet",
							seqs,
							w[g1],
							type,
							processors)
						subject <- .Call("subsetXStringSet",
							seqs,
							w[g2],
							type,
							processors)
						temp <- do.call(AlignProfiles,
							args=c(list(pattern=pattern,
									subject=subject,
									p.weight=p.weight,
									s.weight=s.weight,
									p.struct=structures[w[g1]],
									s.struct=structures[w[g2]],
									substitutionMatrix=sM,
									processors=processors,
									gapOpening=GO[i],
									gapExtension=GE[i]),
								args))
					} else { # do not need to specify substitutionMatrix
						pattern <- .Call("subsetXStringSet",
							seqs,
							w[g1],
							type,
							processors)
						subject <- .Call("subsetXStringSet",
							seqs,
							w[g2],
							type,
							processors)
						temp <- AlignProfiles(pattern=pattern,
							subject=subject,
							p.weight=p.weight,
							s.weight=s.weight,
							processors=processors,
							gapOpening=GO[i],
							gapExtension=GE[i],
							...)
					}
					if (cutoffs[i] <= levels[1] ||
						(length(g1) + length(g2)) < levels[4]) {
						seqs[w[c(g1, g2)]] <- temp
					} else {
						if (type==3L) { # AAStringSet
							seqs[w[c(g1, g2)]] <- FUN(temp,
								substitutionMatrix=sM,
								weight=weights[w[c(g1, g2)]]/mean(weights[w[c(g1, g2)]]),
								processors=processors)
						} else { # DNAStringSet or RNAStringSet
							seqs[w[c(g1, g2)]] <- FUN(temp,
								weight=weights[w[c(g1, g2)]]/mean(weights[w[c(g1, g2)]]),
								processors=processors)
						}
					}
					
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
					
					if (verbose) {
						percentComplete <- as.integer(100*steps/l)
						if (percentComplete > before) {
							setTxtProgressBar(pBar, percentComplete)
							before <- percentComplete
						}
					}
					g1 <- c(g1, g2)
				}
			}
		}
	}
	
	if (refinements==0 && iterations==0) {
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			flush.console()
		}
		
		names(seqs) <- ns
		return(seqs)
	}
	
	for (it in seq_len(iterations)) {
		myXStringSet <- seqs_original
		
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\nDetermining distance matrix based on alignment")
			if (iterations > 1)
				 cat(" (iteration ",
					it,
					" of ",
					iterations,
					")",
					sep="")
			cat(":\n")
			flush.console()
		}
		
		suppressWarnings(d <- DistanceMatrix(seqs,
			verbose=verbose,
			processors=processors,
			includeTerminalGaps=TRUE))
		
		if (verbose) {
			cat("Reclustering into groups by similarity")
			if (iterations > 1)
				 cat(" (iteration ",
					it,
					" of ",
					iterations,
					")",
					sep="")
			cat(":\n")
			flush.console()
		}
		orgTree <- guideTree
		cutoffs <- sort(unique(c(seq(0, 1, 0.01),
			round(quantile(d, seq(0, 1, 0.01), na.rm=TRUE), 5))))
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- IdClusters(d,
			cutoff=cutoffs,
			verbose=verbose,
			method="UPGMA",
			processors=processors))
		guideTree$Top <- 1
		guideTree <- as.matrix(guideTree)
		guideTree[, 1] <- 1:dim(guideTree)[1]
		cutoffs <- c(cutoffs, 1)
		w <- which(is.na(d))
		if (length(w) > 0)
			d[w] <- 1
		
		# calculate weights based on branch lengths
		lastSplit <- numeric(dim(guideTree)[1])
		weights <- numeric(dim(guideTree)[1])
		numGroups <- 1
		for (i in (dim(guideTree)[2] - 1):1) {
			if (length(unique(guideTree[, i])) > numGroups) {
				if (numGroups==1) { # mid-point root of tree
					lastSplit[] <- cutoffs[i]
				} else {
					for (j in unique(guideTree[, i + 1])) {
						w <- which(guideTree[, i + 1]==j)
						if (length(unique(guideTree[w, i])) > 1) {
							weights[w] <- weights[w] + (lastSplit[w] - cutoffs[i])/length(w)
							lastSplit[w] <- cutoffs[i]
						}
					}
				}
				numGroups <- length(unique(guideTree[, i]))
			}
		}
		weights <- weights + lastSplit
		weights <- ifelse(weights==0, 1, weights)
		weights <- weights/mean(weights)
		
		if (verbose) {
			time.1 <- Sys.time()
			cat("Realigning Sequences")
			if (iterations > 1)
				 cat(" (iteration ",
					it,
					" of ",
					iterations,
					")",
					sep="")
			cat(":\n")
			flush.console()
			pBar <- txtProgressBar(style=3, max=100)
			percentComplete <- before <- 0L
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
		
		GO <- gapOpeningMin + gapOpeningSlope*cutoffs
		GE <- gapExtensionMin + gapExtensionSlope*cutoffs
		
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
						if (length(grep(mergers2[steps], mergers, fixed=TRUE)) > 0) {
							myXStringSet[w[1:j]] <- .Call("removeCommonGaps",
								seqs[w[1:j]],
								type,
								processors,
								PACKAGE="DECIPHER")
							next
						}
						
						p.weight <- weights[w[1:(j - 1)]]
						p.weight <- p.weight/mean(p.weight)
						s.weight <- weights[w[j]]
						s.weight <- s.weight/mean(s.weight)
						
						if (subM) {
							pattern <- .Call("subsetXStringSet",
								myXStringSet,
								w[1:(j - 1)],
								type,
								processors)
							subject <- .Call("subsetXStringSet",
								myXStringSet,
								w[j],
								type,
								processors)
							myXStringSet[w[1:j]] <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										p.struct=structures[w[1:(j - 1)]],
										s.struct=structures[w[j]],
										substitutionMatrix=sM,
										processors=processors,
										gapOpening=GO[i],
										gapExtension=GE[i]),
									args))
						} else { # do not need to specify substitutionMatrix
							pattern <- .Call("subsetXStringSet",
								myXStringSet,
								w[1:(j - 1)],
								type,
								processors)
							subject <- .Call("subsetXStringSet",
								myXStringSet,
								w[j],
								type,
								processors)
							myXStringSet[w[1:j]] <- AlignProfiles(pattern=pattern,
								subject=subject,
								p.weight=p.weight,
								s.weight=s.weight,
								processors=processors,
								gapOpening=GO[i],
								gapExtension=GE[i],
								...)
						}
						
						if (verbose) {
							percentComplete <- as.integer(100*steps/l)
							if (percentComplete > before) {
								setTxtProgressBar(pBar, percentComplete)
								before <- percentComplete
							}
						}
					}
				} else {
					g <- unique(guideTree[w, i - 1]) # sub-groups
					if (length(g)==1)
						next
					g <- g[order(width(myXStringSet)[match(g, guideTree[w, i - 1])],
						decreasing=TRUE)] # sort sub-groups by increasing width
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
							if (length(grep(mergers2[steps], mergers, fixed=TRUE)) > 0) {
								myXStringSet[c(w[g1], w[g2])] <- .Call("removeCommonGaps",
									seqs[c(w[g1], w[g2])],
									type,
									processors,
									PACKAGE="DECIPHER")
								
								if (verbose) {
									percentComplete <- as.integer(100*steps/l)
									if (percentComplete > before) {
										setTxtProgressBar(pBar, percentComplete)
										before <- percentComplete
									}
								}
								g1 <- c(g1, g2)
								next
							}
						}
						
						p.weight <- weights[w[g1]]
						p.weight <- p.weight/mean(p.weight)
						s.weight <- weights[w[g2]]
						s.weight <- s.weight/mean(s.weight)
						
						if (subM) {
							pattern <- .Call("subsetXStringSet",
								myXStringSet,
								w[g1],
								type,
								processors)
							subject <- .Call("subsetXStringSet",
								myXStringSet,
								w[g2],
								type,
								processors)
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										p.struct=structures[w[g1]],
										s.struct=structures[w[g2]],
										substitutionMatrix=sM,
										processors=processors,
										gapOpening=GO[i],
										gapExtension=GE[i]),
									args))
						} else { # do not need to specify substitutionMatrix
							pattern <- .Call("subsetXStringSet",
								myXStringSet,
								w[g1],
								type,
								processors)
							subject <- .Call("subsetXStringSet",
								myXStringSet,
								w[g2],
								type,
								processors)
							temp <- AlignProfiles(pattern=pattern,
								subject=subject,
								p.weight=p.weight,
								s.weight=s.weight,
								processors=processors,
								gapOpening=GO[i],
								gapExtension=GE[i],
								...)
						}
						if (cutoffs[i] <= levels[2] ||
							(length(g1) + length(g2)) < levels[4]) {
							myXStringSet[w[c(g1, g2)]] <- temp
						} else {
							if (type==3L) { # AAStringSet
								myXStringSet[w[c(g1, g2)]] <- FUN(temp,
									substitutionMatrix=sM,
									weight=weights[w[c(g1, g2)]]/mean(weights[w[c(g1, g2)]]),
									processors=processors)
							} else { # DNAStringSet or RNAStringSet
								myXStringSet[w[c(g1, g2)]] <- FUN(temp,
									weight=weights[w[c(g1, g2)]]/mean(weights[w[c(g1, g2)]]),
									processors=processors)
							}
						}
						
						if (verbose) {
							percentComplete <- as.integer(100*steps/l)
							if (percentComplete > before) {
								setTxtProgressBar(pBar, percentComplete)
								before <- percentComplete
							}
						}
						g1 <- c(g1, g2)
					}
				}
			}
		}
		
		if (it < iterations && all(seqs==myXStringSet))
			break
		
		seqs <- myXStringSet
		mergers <- mergers2
	}
	
	myXStringSet <- seqs
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		
		if (iterations > 0 && it < iterations) {
			cat("\nAlignment converged - skipping remaining",
				ifelse(iterations - it > 1,
					"iterations.\n",
					"iteration.\n"))
			if (refinements==0)
				cat("\n")
		}
	}
	
	if (refinements > 0) {
		if (type==3L) {
			functionCall <- "colScoresAA"
		} else {
			functionCall <- "colScores"
		}
		GO <- gapOpeningMax/2
		colScores <- function(seqs, structs, weights) {
			scores <- .Call(functionCall,
				seqs,
				sM,
				GO,
				gapExtensionMax,
				weights/mean(weights),
				structs,
				structureMatrix)
			return(sum(scores))
		}
		
		# define trustworthy groups
		if (cluster || iterations > 0) {
			if (type==3L) { # AAStringSet
				w <- which(cutoffs < ifelse(iterations > 0,
					0.7, # fraction identical
					0.9)) # fraction ordered k-mers
			} else { # DNAStringSet or RNAStringSet
				w <- which(cutoffs < ifelse(iterations > 0,
					0.4, # fraction identical
					0.7)) # fraction ordered k-mers
			}
			groups <- guideTree[, w[length(w)]]
		} else { # the only guideTree was provided as input
			groups <- seq_along(myXStringSet)
		}
		
		# refinement
		u <- unique(groups)
		if (verbose) {
			time.1 <- Sys.time()
			cat("\nRefining the alignment:\n")
			flush.console()
			pBar <- txtProgressBar(style=3)
		}
		if (length(u) > 2) { # more than 2 groups
			score <- colScores(myXStringSet, structures, weights)
			for (ref in seq_len(refinements)) {
				org_score <- score
				count <- 0L
				for (i in seq_along(u)) {
					w <- which(groups==u[i])
					seqs1 <- .Call("subsetXStringSet",
						myXStringSet,
						seq_len(length(myXStringSet))[-w],
						type,
						processors)
					seqs1 <- .Call("removeCommonGaps",
						seqs1,
						type,
						processors,
						PACKAGE="DECIPHER")
					seqs2 <- .Call("subsetXStringSet",
						myXStringSet,
						w,
						type,
						processors)
					seqs2 <- .Call("removeCommonGaps",
						seqs2,
						type,
						processors,
						PACKAGE="DECIPHER")
					
					p.weight <- weights[-w]
					p.weight <- p.weight/mean(p.weight)
					s.weight <- weights[w]
					s.weight <- s.weight/mean(s.weight)
					
					if (subM) {
						temp <- do.call(AlignProfiles,
							args=c(list(pattern=seqs1,
									subject=seqs2,
									p.weight=p.weight,
									s.weight=s.weight,
									p.struct=structures[seq_len(length(myXStringSet))[-w]],
									s.struct=structures[w],
									substitutionMatrix=sM,
									processors=processors,
									gapOpening=gapOpeningMax,
									gapExtension=gapExtensionMax),
								args))
					} else {
						temp <- AlignProfiles(pattern=seqs1,
							subject=seqs2,
							p.weight=p.weight,
							s.weight=s.weight,
							processors=processors,
							gapOpening=gapOpeningMax,
							gapExtension=gapExtensionMax,
							...)
					}
					
					if (is.null(structures)) {
						temp_score <- colScores(temp, structures, c(weights[-w], weights[w]))
					} else {
						temp_score <- colScores(temp, structures[c((1:length(myXStringSet))[-w], w)], c(weights[-w], weights[w]))
					}
					if (temp_score > score) {
						score <- temp_score
						myXStringSet <- .Call("subsetXStringSet",
							temp,
							order(c((1:length(myXStringSet))[-w], w)),
							type,
							processors)
						
						count <- count + 1L
						if (count %% levels[3] ||
							length(myXStringSet) < levels[4])
							next # refine every nth change
						
						if (type==3L) { # AAStringSet
							temp <- FUN(myXStringSet,
								substitutionMatrix=sM,
								weight=weights,
								processors=processors)
							temp_score <- colScores(temp, structures, weights)
						} else { # DNAStringSet or RNAStringSet
							temp <- FUN(myXStringSet,
								weight=weights,
								processors=processors)
							temp_score <- colScores(temp, structures, weights)
						}
						
						if (temp_score > score) {
							score <- temp_score
							myXStringSet <- temp
						}
					}
					
					if (verbose)
						setTxtProgressBar(pBar,
							(i + length(u)*(ref - 1))/(length(u)*refinements))
				}
				
				if (org_score==score) # no changes
					break
			}
		} else {
			ref <- 1
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
			
			if (refinements > 0 && ref < refinements)
				cat("Alignment converged - skipping remaining",
					ifelse(refinements - ref > 1,
						"refinements.\n\n",
						"refinement.\n\n"))
		}
	}
	
	if (length(myXStringSet) >= levels[4]) {
		# apply a final adjustment
		if (type==3L) { # AAStringSet
			myXStringSet <- FUN(myXStringSet,
				substitutionMatrix=sM,
				weight=weights,
				processors=processors)
		} else { # DNAStringSet or RNAStringSet
			myXStringSet <- FUN(myXStringSet,
				weight=weights,
				processors=processors)
		}
	}
	
	names(myXStringSet) <- ns
	
	return(myXStringSet)
}
