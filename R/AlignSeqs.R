AlignSeqs <- function(myXStringSet,
	guideTree=NULL,
	iterations=1,
	refinements=1,
	gapOpening=c(-16, -12),
	gapExtension=c(-2, -1),
	useStructures=TRUE,
	structures=NULL,
	FUN=AdjustAlignment,
	levels=c(0.9, 0.7, 0.7, 0.4, 10, 5, 5, 2),
	processors=1,
	verbose=TRUE,
	...) {
	
	# error checking
	type <- switch(class(myXStringSet),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet."))
	l <- length(myXStringSet)
	if (l < 2)
		stop("At least two sequences are required.")
	if (!is.null(guideTree)) {
		if (class(guideTree)=="dendrogram") {
			if (attr(guideTree, "height") > 0.501)
				stop("Total height of guideTree can be at most 0.5.")
			heights <- rapply(guideTree,
				function(x)
					attr(x, "height"))
			if (any(heights > 0.001) || any(heights < -0.001))
				stop("All tips of guideTree must have a height of zero.")
			labels <- rapply(guideTree,
				function(x)
					attr(x, "label"))
			if (any(duplicated(labels)))
				stop("Leaf labels in guideTree must be unique.")
			labels <- match(labels, names(myXStringSet))
			if (any(is.na(labels)))
				stop("Leaf labels in guideTree must match names of myXStringSet.")
			guideTree <- dendrapply(guideTree,
				function(x) {
					if (is.leaf(x))
						x[1] <- match(attr(x, "label"),
							names(myXStringSet))
					return(x)
				})
		} else {
			stop("guideTree must be a dendrogram object.")
		}
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
	if (!is.logical(useStructures))
		stop("useStructures must be a logical.")
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
	if (length(levels) != 8)
		stop("levels must be length 8.")
	if (levels[1] < 0 || levels[1] > 1)
		stop("levels[1] must be between zero and one (inclusive).")
	if (levels[2] < 0 || levels[2] > 1)
		stop("levels[2] must be between zero and one (inclusive).")
	if (levels[3] < 0 || levels[3] > 1)
		stop("levels[3] must be between zero and one (inclusive).")
	if (levels[4] < 0 || levels[4] > 1)
		stop("levels[4] must be between zero and one (inclusive).")
	if (levels[5] <= 0 || is.infinite(levels[5]))
		stop("levels[5] must be between zero and Inf (exclusive).")
	if (levels[5] != floor(levels[5]))
		stop("levels[5] must be a whole number.")
	if (levels[6] != floor(levels[6]))
		stop("levels[6] must be a whole number.")
	if (levels[6] < 0)
		stop("levels[6] must be at least zero.")
	if (levels[7] <= 0)
		stop("levels[7] must be at least zero.")
	if (levels[8] <= 0)
		stop("levels[8] must be at least zero.")
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
	
	# prepare structures and structure matrix
	if (useStructures) {
		if (type==3L) {
			if (is.null(structures)) {
				structures <- PredictHEC(myXStringSet,
					type="probabilities")
			} else {
				if (length(structures) != l)
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
						stop("structureMatrix must contain numerics.")
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
		} else {
			if (!is.null(structures)) {
				if (length(structures) != l)
					stop("structures is not the same length as myXStringSet.")
				if (typeof(structures) != "list")
					stop("structures must be a list.")
				
				w <- which(m=="structureMatrix")
				if (length(w) > 0) {
					structureMatrix <- args[[w]]
					# assume structures and matrix are ordered the same
					if (!is.double(structureMatrix))
						stop("structureMatrix must contain numerics.")
					if (!is.matrix(structureMatrix))
						stop("structureMatrix must be a matrix.")
					if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
						stop("structureMatrix is not square.")
				} else {
					stop("structureMatrix must be specified when structures are provided.")
				}
				if (dim(structureMatrix)[1] != dim(structures[[1]])[1])
					stop("Dimensions of structureMatrix are incompatible with the structures.")
				replace <- TRUE
			} else if (type==2L) {
				w <- which(m=="structureMatrix")
				if (length(w) > 0) {
					structureMatrix <- args[[w]]
					# assume thre matrix is in the correct order
					if (!is.double(structureMatrix))
						stop("structureMatrix must contain numerics.")
					if (!is.matrix(structureMatrix))
						stop("structureMatrix must be a matrix.")
					if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
						stop("structureMatrix is not square.")
					if (dim(structureMatrix)[1] != 3)
						stop("structureMatrix must be 3 x 3 when structures is NULL.")
				} else { # use the default structureMatrix
					structureMatrix <- matrix(c(0, 1, 1, 1, 10, -5, 1, -5, 10),
						nrow=3) # order is ., (, )
				}
				replace <- FALSE
			} else {
				structureMatrix <- numeric()
				useStructures <- FALSE
			}
		}
	} else {
		w <- which(m=="structureMatrix")
		if (length(w) > 0)
			stop("structureMatrix provided when useStructures is FALSE.")
		structureMatrix <- numeric() # needed for colScores
		if (!is.null(structures))
			stop("structures provided when useStructures is FALSE.")
	}
	
	# prepare substitution matrix
	if (type==3L) {
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
		bases <- c("A", "C", "G",
			ifelse(type==2L, "U", "T"))
		w <- which(m=="substitutionMatrix")
		if (length(w) > 0) {
			subM <- TRUE
			sM <- args[[w]]
			args[w] <- NULL
			if (is.matrix(sM)) {
				if (any(!(bases %in% dimnames(sM)[[1]])) ||
					any(!(bases %in% dimnames(sM)[[2]])))
					stop("substitutionMatrix is incomplete.")
				sM <- sM[bases, bases]
				sM <- sM + 0 # convert to numeric matrix
				if (type==2L) { # RNAStringSet
					sM2 <- sM
					dimnames(sM2) <- list(DNA_BASES, DNA_BASES)
				}
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		} else if (type==2L &&
			!("misMatch" %in% m) &&
			!("perfectMatch" %in% m)) {
			subM <- TRUE
			sM <- matrix(c(11, 3, 5, 4, 3, 12, 3, 6, 5, 3, 12, 3, 4, 6, 3, 10),
				nrow=4,
				ncol=4,
				dimnames=list(bases, bases))
			sM2 <- matrix(c(11, 3, 5, 4, 3, 12, 3, 6, 5, 3, 12, 3, 4, 6, 3, 10),
				nrow=4,
				ncol=4,
				dimnames=list(DNA_BASES, DNA_BASES))
		} else {
			subM <- FALSE
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
	
	# need to maximize the recursion depth temporarily
	org.options <- options(expressions=5e5)
	on.exit(options(org.options))
	
	if (length(args)==0)
		args <- NULL
	
	if (verbose)
		time.1 <- Sys.time()
	
	w.x <- width(myXStringSet)
	if (type==3L) {
		wordSize <- ceiling(log(100*mean(w.x), 9)) - 1
		if (wordSize > 10)
			wordSize <- 10
		if (wordSize < 1)
			wordSize <- 1
	} else {
		wordSize <- ceiling(log(100*mean(w.x), 4)) - 1
		if (wordSize > 15)
			wordSize <- 15
		if (wordSize < 5)
			wordSize <- 5
	}
	
	if (is.null(guideTree)) {
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
			v <- .Call("enumerateSequenceReducedAA",
				myXStringSet,
				wordSize,
				c(A=2L, R=7L, N=6L, D=8L, C=3L, Q=4L, E=8L,
					G=5L, H=4L, I=1L, L=1L, K=7L, M=0L, F=0L,
					P=4L, S=6L, T=6L, W=4L, Y=4L, V=1L),
				PACKAGE="DECIPHER")
		} else { # DNAStringSet or RNAStringSet
			v <- .Call("enumerateSequence",
				myXStringSet,
				wordSize,
				PACKAGE="DECIPHER")
		}
		
		d <- .Call("matchOrder",
			v,
			verbose,
			pBar,
			processors,
			PACKAGE="DECIPHER")
		
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
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- IdClusters(d,
			method="single",
			type="dendrogram",
			verbose=verbose,
			processors=processors))
	}
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Aligning Sequences:\n")
		flush.console()
		pBar <- txtProgressBar(style=3, max=100)
		before <- steps <- 0L
		nsteps <- l - 1L
	}
	
	.align <- function(dend) {
		l <- length(dend)
		treeLengths <- numeric(l)
		inherit <- attr(dend, "inherit")
		if (!is.null(inherit)) {
			if (inherit) {
				u <- unlist(dend)
				seqs <<- .replace(seqs,
					.Call("removeCommonGaps",
						.subset(seqs_prev, u),
						type,
						processors,
						PACKAGE="DECIPHER"),
					u)
				
				if (verbose) {
					steps <<- steps + 1L
					percentComplete <- as.integer(100L*steps/nsteps)
					if (percentComplete > before) {
						setTxtProgressBar(pBar, percentComplete)
						before <<- percentComplete
					}
				}
			}
			
			if (l > 1) {
				for (i in seq_len(l))
					treeLengths[i] <- .align(dend[[i]])
				
				h <- attr(dend, "height")
				for (i in seq_len(l)) {
					m <- unlist(dend[i])
					treeLengths[i] <- treeLengths[i] + h - heights[m[1]]
					weights[m] <<- weights[m] + (h - heights[m])/length(m)
					heights[m] <<- h
				}
			} else if (is.leaf(dend)) {
				heights[unlist(dend)] <- attr(dend, "height")
			} else {
				treeLengths[1] <- .align(dend[[1]])
			}
		} else if (l > 1) {
			for (i in seq_len(l))
				treeLengths[i] <- .align(dend[[i]])
			
			h <- attr(dend, "height")
			members <- vector("list", l)
			for (i in seq_len(l)) {
				m <- unlist(dend[i])
				treeLengths[i] <- treeLengths[i] + h - heights[m[1]]
				weights[m] <<- weights[m] + (h - heights[m])/length(m)
				heights[m] <<- h
				members[[i]] <- m
			}
			
			h <- h*2 # total length of sub-tree
			GO <- h*gapOpeningSlope + gapOpeningMin
			GE <- h*gapExtensionSlope + gapExtensionMin
			
			for (i in 2:length(dend)) {
				x <- unlist(members[1:(i - 1)])
				y <- members[[i]]
				combo <- c(x, y)
				
				p.weight <- weights[x]
				w <- which(p.weight <= 0)
				if (length(w) > 0)
					p.weight[w] <- 1
				p.weight <- p.weight/mean(p.weight)
				s.weight <- weights[y]
				w <- which(s.weight <= 0)
				if (length(w) > 0)
					s.weight[w] <- 1
				s.weight <- s.weight/mean(s.weight)
				
				pattern <- .subset(seqs, x)
				subject <- .subset(seqs, y)
				
				if (subM) {
					if (useStructures) {
						if (is.null(structures)) {
							if (type==2L && h < LEVEL) {
								# align as DNAStringSet (faster because no pairs matrix)
								temp <- .switch(do.call(AlignProfiles,
									args=c(list(pattern=.switch(pattern),
											subject=.switch(subject),
											p.weight=p.weight,
											s.weight=s.weight,
											substitutionMatrix=sM2,
											processors=processors,
											gapOpening=GO,
											gapExtension=GE),
										args)))
							} else {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											substitutionMatrix=sM,
											processors=processors,
											gapOpening=GO,
											gapExtension=GE),
										args))
							}
						} else {
							if (type==2L && h < LEVEL) {
								# align as DNAStringSet (faster because no pairs matrix)
								temp <- .switch(do.call(AlignProfiles,
									args=c(list(pattern=.switch(pattern),
											subject=.switch(subject),
											p.weight=p.weight,
											s.weight=s.weight,
											p.struct=structures[x],
											s.struct=structures[y],
											substitutionMatrix=sM2,
											processors=processors,
											gapOpening=GO,
											gapExtension=GE),
										args)))
							} else {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											p.struct=structures[x],
											s.struct=structures[y],
											substitutionMatrix=sM,
											processors=processors,
											gapOpening=GO,
											gapExtension=GE),
										args))
							}
						}
					} else {
						if (type==2L && h < LEVEL) {
							# align as DNAStringSet (faster because no pairs matrix)
							temp <- .switch(do.call(AlignProfiles,
								args=c(list(pattern=.switch(pattern),
										subject=.switch(subject),
										p.weight=p.weight,
										s.weight=s.weight,
										substitutionMatrix=sM2,
										processors=processors,
										gapOpening=GO,
										gapExtension=GE),
									args)))
						} else {
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										substitutionMatrix=sM,
										processors=processors,
										gapOpening=GO,
										gapExtension=GE),
									args))
						}
					}
				} else {
					if (useStructures) {
						temp <- do.call(AlignProfiles,
							args=c(list(pattern=pattern,
									subject=subject,
									p.weight=p.weight,
									s.weight=s.weight,
									p.struct=structures[x],
									s.struct=structures[y],
									processors=processors,
									gapOpening=GO,
									gapExtension=GE),
								args))
					} else {
						temp <- do.call(AlignProfiles,
							args=c(list(pattern=pattern,
									subject=subject,
									p.weight=p.weight,
									s.weight=s.weight,
									processors=processors,
									gapOpening=GO,
									gapExtension=GE),
								args))
					}
				}
				
				if (h > LEVEL &&
					length(temp) >= levels[6]) {
					weight <- weights[combo]
					w <- which(weight <= 0)
					if (length(w) > 0)
						weight[w] <- 1
					weight <- weight/mean(weight)
					
					if (subM) {
						temp <- FUN(temp,
							substitutionMatrix=sM,
							weight=weight,
							processors=processors)
					} else {
						temp <- FUN(temp,
							weight=weight,
							processors=processors)
					}
				}
				
				seqs <<- .replace(seqs,
					temp,
					combo)
				
				if (verbose) {
					steps <<- steps + 1L
					percentComplete <- as.integer(100L*steps/nsteps)
					if (percentComplete > before) {
						setTxtProgressBar(pBar, percentComplete)
						before <<- percentComplete
					}
				}
			}
		} else if (is.leaf(dend)) {
			heights[unlist(dend)] <- attr(dend, "height")
		} else {
			treeLengths[1] <- .align(dend[[1]])
		}
		
		return(sum(treeLengths))
	}
	
	.reorder <- function(dend) {
		l <- length(dend)
		if (l > 1) {
			for (i in seq_len(l))
				dend[[i]] <- .reorder(dend[[i]])
			
			members <- lapply(dend, unlist)
			# sort tree by descending width
			o <- order(sapply(members,
					function(x)
						max(w.x[x])),
				lengths(members),
				sapply(members, min),
				decreasing=TRUE)
			dend[] <- dend[o]
		} else if (!is.leaf(dend)) {
			dend[[1]] <- .reorder(dend[[1]])
		}
		return(dend)
	}
	guideTree <- .reorder(guideTree)
	
	ns <- names(myXStringSet)
	seqs <- myXStringSet
	weights <- heights <- numeric(l)
	if (type==3L) {
		LEVEL <- levels[1]
	} else {
		LEVEL <- levels[3]
	}
	treeLength <- .align(guideTree)
	
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
			cat("\n")
		}
		
		names(seqs) <- ns
		return(seqs)
	}
	
	.RNAStructures <- function(seqs, weights) {
		if (verbose) {
			cat("\nComputing RNA Secondary Structures:\n")
			flush.console()
		}
		
		w <- which(weights <= 0)
		if (length(w) > 0)
			weights[w] <- 1
		weights <- weights/mean(weights)
		
		if (replace) {
			structureMatrix <- matrix(c(0, 1, 1, 1, 10, -5, 1, -5, 10),
				nrow=3) # order is ., (, )
			replace <- FALSE
		}
		
		PredictDBN(seqs,
			type="structures",
			weight=weights,
			processors=processors,
			verbose=verbose)
	}
	
	minTreeLength <- levels[7]
	
	if (iterations > 0) {
		if (type==3L) {
			LEVEL <- levels[2]
		} else {
			LEVEL <- levels[4]
		}
		
		.order <- function(dend) {
			l <- length(dend)
			if (l > 1) {
				j <<- j + 1L
				orders[[j]] <<- unlist(dend)
				
				for (i in seq_len(l))
					.order(dend[[i]])
			} else if (!is.leaf(dend)) {
				.order(dend[[1]])
			}
		}
		
		# record the alignment order in the original tree
		j <- 0L
		orders <- vector("list", l - 1L)
		.order(guideTree)
		
		.compare <- function(dend, found=FALSE) {
			l <- length(dend)
			if (found) {
				if (l > 1) {
					if (verbose)
						nsteps <<- nsteps - 1L
					
					j <<- j + 1L
					orders[[j]] <<- unlist(dend)
					
					attr(dend, "inherit") <- FALSE
					
					for (i in seq_len(l))
						dend[[i]] <- .compare(dend[[i]], found)
				} else if (!is.leaf(dend)) {
					dend[[1]] <- .compare(dend[[1]], found)
				}
			} else if (l > 1) {
				o <- unlist(dend)
				j <<- j + 1L
				orders[[j]] <<- o
				
				found <- FALSE
				w <- which(ls==length(o))
				for (i in seq_along(w)) {
					if (all(orders_prev[[w[i]]]==o)) {
						found <- TRUE
						break
					}
				}
				
				if (found)
					attr(dend, "inherit") <- TRUE
				
				for (i in seq_len(l))
					dend[[i]] <- .compare(dend[[i]], found)
			} else if (!is.leaf(dend)) {
				dend[[i]] <- .compare(dend[[1]], found)
			}
			
			return(dend)
		}
	}
	
	for (it in seq_len(iterations)) {
		seqs_prev <- seqs
		ls <- lengths(orders)
		orders_prev <- orders
		
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			if (iterations > 1)
				 cat("\nIteration ",
					it,
					" of ",
					iterations,
					":\n",
					sep="")
		}
		
		if (type==2L &&
			treeLength >= minTreeLength &&
			useStructures) {
			structures <- .RNAStructures(seqs, weights)
		} else if (verbose) {
			cat("\n")
		}
		minTreeLength <- levels[8]
		
		if (verbose) {
			cat("Determining distance matrix based on alignment:\n")
			flush.console()
		}
		
		d <- DistanceMatrix(seqs,
			verbose=verbose,
			processors=processors,
			includeTerminalGaps=TRUE)
		
		if (verbose) {
			cat("Reclustering into groups by similarity:\n")
			flush.console()
		}
		
		orgTree <- guideTree
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- IdClusters(d,
			method="UPGMA",
			type="dendrogram",
			collapse=0,
			verbose=verbose,
			processors=processors))
		
		if (verbose) {
			time.1 <- Sys.time()
			cat("Realigning Sequences:\n")
			flush.console()
			pBar <- txtProgressBar(style=3, max=100)
			before <- steps <- 0L
			nsteps <- l - 1L
		}
		
		j <- 0L
		orders <- vector("list", l - 1L)
		guideTree <- .reorder(guideTree)
		guideTree <- .compare(guideTree)
		
		seqs <- myXStringSet
		weights <- heights <- numeric(l)
		treeLength <- .align(guideTree)
		
		if (it < iterations && all(seqs==seqs_prev))
			break
	}
	
	myXStringSet <- seqs
	w <- which(weights <= 0)
	if (length(w) > 0)
		weights[w] <- 1
	weights <- weights/mean(weights)
	
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
		}
	}
	
	if (refinements > 0) {
		if (type==3L) {
			functionCall <- "colScoresAA"
		} else {
			functionCall <- "colScores"
		}
		
		GO <- gapOpeningMax/2 # applied at both ends
		colScores <- function(seqs, structs, weights) {
			scores <- .Call(functionCall,
				seqs,
				sM,
				GO,
				gapExtensionMax,
				weights,
				structs,
				structureMatrix)
			return(sum(scores))
		}
		
		# define trustworthy groups
		if (type==3L) { # AAStringSet
			cutoff <- ifelse(iterations > 0,
				levels[2]/2, # (fraction identical)/2
				levels[1]/2) # (fraction ordered k-mers)/2
		} else { # DNAStringSet or RNAStringSet
			cutoff <- ifelse(iterations > 0,
				levels[4]/2, # (fraction identical)/2
				levels[3]/2) # (fraction ordered k-mers)/2
		}
		guideTree <- cut(guideTree, cutoff)$lower
		
		# refinement
		n <- length(guideTree)
		if (n > 2) { # more than 2 groups
			if (type==2L &&
				treeLength >= minTreeLength &&
				(iterations==0 || (iterations > 0 && !all(seqs==seqs_prev))) &&
				useStructures) {
				structures <- .RNAStructures(seqs, weights)
			} else if (verbose) {
				cat("\n")
			}
			
			if (verbose) {
				time.1 <- Sys.time()
				cat("Refining the alignment:\n")
				flush.console()
				pBar <- txtProgressBar(style=3)
			}
			
			score <- colScores(myXStringSet, structures, weights)
			vec <- seq_along(myXStringSet)
			
			for (ref in seq_len(refinements)) {
				org_score <- score
				count <- 0L
				for (i in seq_len(n)) {
					x <- unlist(guideTree[[i]])
					y <- vec[-x]
					o <- c(x, y)
					
					pattern <- .subset(myXStringSet, x)
					pattern <- .Call("removeCommonGaps",
						pattern,
						type,
						processors,
						PACKAGE="DECIPHER")
					subject <- .subset(myXStringSet, y)
					subject <- .Call("removeCommonGaps",
						subject,
						type,
						processors,
						PACKAGE="DECIPHER")
					
					p.weight <- weights[x]
					p.weight <- p.weight/mean(p.weight)
					s.weight <- weights[y]
					s.weight <- s.weight/mean(s.weight)
					
					if (subM) {
						if (useStructures) {
							if (is.null(structures)) {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											substitutionMatrix=sM,
											processors=processors,
											gapOpening=gapOpeningMax,
											gapExtension=gapExtensionMax),
										args))
							} else {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											p.struct=structures[x],
											s.struct=structures[y],
											substitutionMatrix=sM,
											processors=processors,
											gapOpening=gapOpeningMax,
											gapExtension=gapExtensionMax),
										args))
							}
						} else {
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										substitutionMatrix=sM,
										processors=processors,
										gapOpening=gapOpeningMax,
										gapExtension=gapExtensionMax),
									args))
						}
					} else {
						if (useStructures) {
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										p.struct=structures[x],
										s.struct=structures[y],
										processors=processors,
										gapOpening=gapOpeningMax,
										gapExtension=gapExtensionMax),
									args))
						} else {
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										processors=processors,
										gapOpening=gapOpeningMax,
										gapExtension=gapExtensionMax),
									args))
						}
					}
					
					if (useStructures) {
						temp_score <- colScores(temp, structures[o], weights[o])
					} else {
						temp_score <- colScores(temp, NULL, weights[o])
					}
					
					if (temp_score > score) {
						score <- temp_score
						myXStringSet <- .subset(temp, order(o))
						
						count <- count + 1L
						if (count %% levels[5] ||
							l < levels[6])
							next # refine every nth change
						
						if (subM) {
							temp <- FUN(myXStringSet,
								substitutionMatrix=sM,
								weight=weights,
								processors=processors)
						} else {
							temp <- FUN(myXStringSet,
								weight=weights,
								processors=processors)
						}
						temp_score <- colScores(temp, structures, weights)
						
						if (temp_score > score) {
							score <- temp_score
							myXStringSet <- temp
						}
					}
					
					if (verbose)
						setTxtProgressBar(pBar,
							(i + n*(ref - 1))/(n*refinements))
				}
				
				if (org_score==score) # no changes
					break
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
		} else if (verbose) {
			cat("\n")
		}
	} else if (verbose) {
		cat("\n")
	}
	
	if (l >= levels[6]) {
		# apply a final adjustment
		if (subM) {
			myXStringSet <- FUN(myXStringSet,
				substitutionMatrix=sM,
				weight=weights,
				processors=processors)
		} else {
			myXStringSet <- FUN(myXStringSet,
				weight=weights,
				processors=processors)
		}
	}
	
	names(myXStringSet) <- ns
	
	return(myXStringSet)
}
