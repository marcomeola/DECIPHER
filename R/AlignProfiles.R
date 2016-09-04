AlignProfiles <- function(pattern,
	subject,
	p.weight=1,
	s.weight=1,
	p.struct=NULL,
	s.struct=NULL,
	perfectMatch=5,
	misMatch=0,
	gapOpening=-13,
	gapExtension=-1,
	gapPower=-1,
	terminalGap=-5,
	restrict=-1000,
	anchor=0.7,
	normPower=1,
	substitutionMatrix=NULL,
	structureMatrix=NULL,
	processors=1) {
	
	# error checking
	type <- switch(class(pattern),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet."))
	if (class(pattern) != class(subject))
		stop("pattern and subject must be of the same class.")
	if (length(subject) < 1)
		stop("At least one sequence is required in the subject.")
	w.p <- unique(width(pattern))
	if (length(w.p)!=1)
		stop("Sequences in pattern must be the same width (aligned).")
	w.s <- unique(width(subject))
	if (length(w.s)!=1)
		stop("Sequences in subject must be the same width (aligned).")
	if (!is.numeric(p.weight))
		stop("p.weight must be a numeric.")
	if (length(p.weight)!=1 && length(p.weight)!=length(pattern))
		stop("Length of p.weight must equal one or the length of the pattern.")
	if (length(p.weight)==1) {
		p.weight <- rep(1, length(pattern))
	} else {
		if (!isTRUE(all.equal(1, mean(p.weight))))
			stop("The mean of p.weight must be 1.")
	}
	if (!is.numeric(s.weight))
		stop("s.weight must be a numeric.")
	if (length(s.weight)!=1 && length(s.weight)!=length(subject))
		stop("Length of s.weight must equal one or the length of the subject.")
	if (length(s.weight)==1) {
		s.weight <- rep(1, length(subject))
	} else {
		if (!isTRUE(all.equal(1, mean(s.weight))))
			stop("The mean of s.weight must be 1.")
	}
	if (is.null(p.struct) != is.null(s.struct))
		stop("Both p.struct and s.struct must be specified.")
	if (!is.null(p.struct)) {
		if (is.matrix(p.struct)) {
			if (dim(p.struct)[2] != w.p)
				stop("The number of columns in p.struct does not match the width of the pattern.")
		} else if (is.list(p.struct)) {
			if (length(p.struct) != length(pattern))
				stop("p.struct is not the same length as the pattern.")
		} else {
			stop("p.struct must be a matrix or list.")
		}
	}
	if (!is.null(s.struct)) {
		if (is.matrix(s.struct)) {
			if (dim(s.struct)[2] != w.s)
				stop("The number of columns in s.struct does not match the width of the subject.")
		} else if (is.list(s.struct)) {
			if (length(s.struct) != length(subject))
				stop("s.struct is not the same length as the subject.")
		} else {
			stop("s.struct must be a matrix or list.")
		}
	}
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	gapOpening <- gapOpening/2 # split into gap opening and closing
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(gapPower))
		stop("gapPower must be a numeric.")
	if (!all(is.numeric(terminalGap)))
		stop("terminalGap must be a numeric.")
	if (length(terminalGap) > 2 || length(terminalGap) < 1)
		stop("Length of terminalGap must be 1 or 2.")
	if (any(is.infinite(terminalGap)))
		stop("terminalGap must be finite.")
	if (length(terminalGap)==1)
		terminalGap[2] <- terminalGap[1]
	if (!all(is.numeric(restrict)))
		stop("restrict must be a numeric.")
	if (restrict >= 0)
		stop("restrict must be less than zero.")
	if (!is.numeric(anchor) && !is.na(anchor))
		stop("anchor must be numeric.")
	if (is.matrix(anchor)) {
		if (any(is.na(anchor)))
			stop("The matrix anchor must not contain NA values.")
		if (dim(anchor)[1] != 4)
			stop("The matrix anchor must have four rows.")
		if (is.unsorted(anchor[1:2,], strictly=TRUE))
			stop("Anchors must be ascending order.")
		if (is.unsorted(anchor[3:4,], strictly=TRUE))
			stop("Anchors must be ascending order.")
		if (dim(anchor)[2] > 0) {
			if (anchor[1, 1] < 1)
				stop("The first anchor is outside the range of the pattern.")
			if (anchor[3, 1] < 1)
				stop("The first anchor is outside the range of the subject.")
			if (anchor[2, dim(anchor)[2]] > w.p)
				stop("The last anchor is outside the range of the pattern.")
			if (anchor[4, dim(anchor)[2]] > w.s)
				stop("The last anchor is outside the range of the subject.")
		}
	} else {
		if (is.numeric(anchor) && anchor <= 0)
			stop("anchor must be greater than zero.")
		if (is.numeric(anchor) && anchor > 1)
			stop("anchor must be less than or equal to one.")
	}
	if (!all(is.numeric(normPower)))
		stop("normPower must be a numeric.")
	if (normPower < 0)
		stop("normPower must be at least zero.")
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
	if (length(pattern) > 2147483647)
		stop(paste("Length of pattern (",
			length(pattern),
			") longer than the maximum allowable length (2,147,483,647).",
			sep=""))
	if (length(subject) > 2147483647)
		stop(paste("Length of subject (",
			length(subject),
			") longer than the maximum allowable length (2,147,483,647).",
			sep=""))
	
	if (type==3L) { # AAStringSet
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.null(substitutionMatrix)) {
			# use MIQS
			substitutionMatrix <- matrix(c(3.2,-1.3,-0.4,-0.4,1.5,-0.2,-0.4,0.4,-1.2,-1.3,-1.4,-0.7,-1,-2.3,-0.1,0.8,0.8,-3.6,-2.4,0,-6.1,-1.3,6.2,-0.1,-1.5,-2.7,1.8,-0.7,-1.9,0.9,-2.4,-2.5,3.3,-1.1,-3.3,-1.1,-0.3,-0.9,-3.8,-1.9,-2.3,-6.1,-0.4,-0.1,5.1,2.6,-1.6,0.9,0.8,0.2,1,-3.6,-3.5,0.7,-2.3,-3.5,-1.4,0.9,0,-4.5,-1.5,-2.6,-6.1,-0.4,-1.5,2.6,5.7,-3.7,0.9,2.7,-0.5,0.3,-4.5,-4.6,0.4,-3.3,-5.8,-0.3,0.3,-0.2,-5.3,-3.9,-3.5,-6.1,1.5,-2.7,-1.6,-3.7,11.7,-2.8,-3.2,-1.7,-1.2,0.2,-2.3,-3.2,0.1,-2.8,-2.8,1,0,-6.1,-0.7,1.8,-6.1,-0.2,1.8,0.9,0.9,-2.8,3.6,2.1,-1.6,1.2,-2.2,-1.9,1.7,-0.4,-2.4,-0.4,0.4,0.1,-5.4,-2.8,-1.8,-6.1,-0.4,-0.7,0.8,2.7,-3.2,2.1,4.3,-1.3,-0.2,-3.3,-2.8,1.1,-2.3,-4.1,0,0.4,-0.2,-5.8,-2.4,-2.3,-6.1,0.4,-1.9,0.2,-0.5,-1.7,-1.6,-1.3,7.6,-1.6,-5.4,-4.8,-1.7,-3.6,-4.6,-1.6,0,-1.9,-4.8,-4.5,-3.8,-6.1,-1.2,0.9,1,0.3,-1.2,1.2,-0.2,-1.6,7.5,-2.2,-1.9,0,-2.1,0,-1.5,0,-0.2,-0.3,2.1,-2.3,-6.1,-1.3,-2.4,-3.6,-4.5,0.2,-2.2,-3.3,-5.4,-2.2,4.6,3.1,-2.3,1.7,0.7,-3.7,-2.8,-0.7,-0.7,-0.8,3.3,-6.1,-1.4,-2.5,-3.5,-4.6,-2.3,-1.9,-2.8,-4.8,-1.9,3.1,4.6,-2.4,3.2,2.1,-2.8,-2.9,-1.6,-0.2,0,2,-6.1,-0.7,3.3,0.7,0.4,-3.2,1.7,1.1,-1.7,0,-2.3,-2.4,3.6,-1.1,-3.7,-0.1,0,0,-4,-2.3,-2,-6.1,-1,-1.1,-2.3,-3.3,0.1,-0.4,-2.3,-3.6,-2.1,1.7,3.2,-1.1,5.4,1.4,-2.8,-1.8,-0.8,-2.1,-0.9,1.4,-6.1,-2.3,-3.3,-3.5,-5.8,-2.8,-2.4,-4.1,-4.6,0,0.7,2.1,-3.7,1.4,7.4,-3.7,-2.6,-2.3,4.2,5.2,-0.3,-6.1,-0.1,-1.1,-1.4,-0.3,-2.8,-0.4,0,-1.6,-1.5,-3.7,-2.8,-0.1,-2.8,-3.7,8.4,-0.1,-0.5,-3.6,-4.5,-2.5,-6.1,0.8,-0.3,0.9,0.3,1,0.4,0.4,0,0,-2.8,-2.9,0,-1.8,-2.6,-0.1,3.1,1.6,-3.5,-1.5,-1.4,-6.1,0.8,-0.9,0,-0.2,0,0.1,-0.2,-1.9,-0.2,-0.7,-1.6,0,-0.8,-2.3,-0.5,1.6,3.8,-5.3,-2.1,-0.1,-6.1,-3.6,-3.8,-4.5,-5.3,-6.1,-5.4,-5.8,-4.8,-0.3,-0.7,-0.2,-4,-2.1,4.2,-3.6,-3.5,-5.3,14.8,4.9,-3.3,-6.1,-2.4,-1.9,-1.5,-3.9,-0.7,-2.8,-2.4,-4.5,2.1,-0.8,0,-2.3,-0.9,5.2,-4.5,-1.5,-2.1,4.9,8.3,-1.2,-6.1,0,-2.3,-2.6,-3.5,1.8,-1.8,-2.3,-3.8,-2.3,3.3,2,-2,1.4,-0.3,-2.5,-1.4,-0.1,-3.3,-1.2,3.5,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,1),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (is.character(substitutionMatrix)) {
			if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
		"PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")))
				stop("Invalid substitutionMatrix.")
		}
		if (is.matrix(substitutionMatrix)) {
			if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix, envir=environment())))
		}
		subMatrix <- subMatrix[AAs, AAs]
		subMatrix <- as.numeric(subMatrix)
	} else {
		if (!is.null(substitutionMatrix)) {
			if (is.matrix(substitutionMatrix)) {
				bases <- c("A", "C", "G",
					ifelse(type==2L, "U", "T"))
				if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
					any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
					stop("substitutionMatrix is incomplete.")
				substitutionMatrix <- substitutionMatrix[bases, bases]
				substitutionMatrix <- as.numeric(substitutionMatrix)
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		} else if (type==2L && missing(perfectMatch) && missing(misMatch)) {
			sM <- matrix(c(11, 3, 5, 4, 3, 12, 3, 6, 5, 3, 12, 3, 4, 6, 3, 10),
				nrow=4,
				ncol=4,
				dimnames=list(bases, bases))
		}
	}
	
	if (type==3L) {
		consensusProfile <- "consensusProfileAA"
	} else {
		consensusProfile <- "consensusProfile"
	}
	if (is.null(p.struct)) {
		p.profile <- .Call(consensusProfile,
			pattern,
			p.weight,
			NULL,
			PACKAGE="DECIPHER")
		s.profile <- .Call(consensusProfile,
			subject,
			s.weight,
			NULL,
			PACKAGE="DECIPHER")
	} else {
		if (is.null(structureMatrix)) {
			if (type==3L) {
				# assume structures from PredictHEC
				structureMatrix <- matrix(c(5, 0, -2, 0, 9, -1, -2, -1, 2),
					nrow=3) # order is H, E, C
			} else {
				structureMatrix <- matrix(c(0, 1, 1, 1, 10, -5, 1, -5, 10),
					nrow=3) # order is ., (, )
			}
		} else {
			# assume structures and matrix are ordered the same
			if (!is.double(structureMatrix))
				stop("structureMatrix must be contain numerics.")
			if (!is.matrix(structureMatrix))
				stop("structureMatrix must be a matrix.")
			if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
				stop("structureMatrix is not square.")
		}
		
		if (is.list(p.struct)) {
			if (dim(structureMatrix)[1] != dim(p.struct[[1]])[1])
				stop("Dimensions of structureMatrix are incompatible with p.struct.")
			
			p.profile <- .Call(consensusProfile,
				pattern,
				p.weight,
				p.struct,
				PACKAGE="DECIPHER")
		} else { # p.struct is a matrix
			if (dim(structureMatrix)[1] != dim(p.struct)[1])
				stop("Dimensions of structureMatrix are incompatible with p.struct.")
			
			p.profile <- .Call(consensusProfile,
				pattern,
				p.weight,
				NULL,
				PACKAGE="DECIPHER")
			
			p.profile <- rbind(p.profile, p.struct)
		}
		if (is.list(s.struct)) {
			if (dim(structureMatrix)[1] != dim(s.struct[[1]])[1])
				stop("Dimensions of structureMatrix are incompatible with s.struct.")
			
			s.profile <- .Call(consensusProfile,
				subject,
				s.weight,
				s.struct,
				PACKAGE="DECIPHER")
		} else { # s.struct is a matrix
			if (dim(structureMatrix)[1] != dim(s.struct)[1])
				stop("Dimensions of structureMatrix are incompatible with s.struct.")
			
			s.profile <- .Call(consensusProfile,
				subject,
				s.weight,
				NULL,
				PACKAGE="DECIPHER")
			
			s.profile <- rbind(s.profile, s.struct)
		}
	}
	
	f <- function(p.profile, s.profile, tGaps=terminalGap) {
		p.d <- dim(p.profile)[2]
		s.d <- dim(s.profile)[2]
		size <- as.numeric(p.d)*as.numeric(s.d)
		if (size > 2147483647) # maximum when indexing by signed integer
			stop(paste("Alignment larger (",
				format(size, big.mark=","),
				") than the maximum allowable size (2,147,483,647).",
				sep=""))
		if (p.d < 100 || s.d < 100) {
			res <- -1e10 # do not restrict
		} else {
			res <- restrict
		}
		
		if (type==3) { # AAStringSet
			if (is.null(p.struct)) {
				t <- .Call("alignProfilesAA",
					p.profile,
					s.profile,
					subMatrix,
					numeric(),
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					res,
					processors,
					PACKAGE="DECIPHER")
			} else {
				t <- .Call("alignProfilesAA",
					p.profile,
					s.profile,
					subMatrix,
					structureMatrix,
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					res,
					processors,
					PACKAGE="DECIPHER")
			}
		} else { # DNAStringSet or RNAStringSet
			if (is.null(p.struct)) {
				t <- .Call("alignProfiles",
					p.profile,
					s.profile,
					type,
					substitutionMatrix,
					numeric(),
					perfectMatch,
					misMatch,
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					res,
					processors,
					PACKAGE="DECIPHER")
			} else {
				t <- .Call("alignProfiles",
					p.profile,
					s.profile,
					type,
					substitutionMatrix,
					structureMatrix,
					perfectMatch,
					misMatch,
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					res,
					processors,
					PACKAGE="DECIPHER")
			}
		}
	}
	
	if (is.na(anchor[1])) { # don't use anchors
		inserts <- f(p.profile, s.profile)
	} else { # use anchors
		if (is.matrix(anchor)) {
			anchors <- anchor
		} else { # find anchors
			if (type==3L) { # AAStringSet
				wordSize <- 7
			} else {
				wordSize <- 15
			}
			l <- min(length(pattern), length(subject))
			o.p <- order(p.weight, decreasing=TRUE)
			o.s <- order(s.weight, decreasing=TRUE)
			if (type==3L) { # AAStringSet
				num.p <- .Call("enumerateGappedSequenceAA",
					pattern,
					wordSize,
					o.p[1:l],
					PACKAGE="DECIPHER")
				num.s <- .Call("enumerateGappedSequenceAA",
					subject,
					wordSize,
					o.s[1:l],
					PACKAGE="DECIPHER")
			} else {
				num.p <- .Call("enumerateGappedSequence",
					pattern,
					wordSize,
					o.p[1:l],
					PACKAGE="DECIPHER")
				num.s <- .Call("enumerateGappedSequence",
					subject,
					wordSize,
					o.s[1:l],
					PACKAGE="DECIPHER")
			}
			
			anchors <- .Call("matchRanges",
				num.p,
				num.s,
				wordSize,
				w.p,
				anchor,
				PACKAGE="DECIPHER")
		}
		
		numAnchors <- dim(anchors)[2]
		if (numAnchors==0) {
			inserts <- f(p.profile, s.profile)
		} else {
			max.p <- which.max(p.weight)
			max.s <- which.max(s.weight)
			
			if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				1L, anchors[1, 1],
				1L, anchors[3, 1],
				max.p,
				max.s,
				PACKAGE="DECIPHER")) {
				inserts <- f(p.profile[, 1L:anchors[1, 1], drop=FALSE],
					s.profile[, 1L:anchors[3, 1], drop=FALSE],
					tGaps=c(terminalGap[1], -1e9))
			} else {
				inserts <- list(integer(), integer(),
					integer(), integer())
			}
			
			n <- 2L
			while (n <= numAnchors) { # align regions between anchors
				if (!.Call("firstSeqsEqual",
					pattern,
					subject,
					anchors[2, n - 1], anchors[1, n],
					anchors[4, n - 1], anchors[3, n],
					max.p,
					max.s,
					PACKAGE="DECIPHER")) {
					temp <- f(p.profile[, anchors[2, n - 1]:anchors[1, n], drop=FALSE],
						s.profile[, anchors[4, n - 1]:anchors[3, n], drop=FALSE],
						tGaps=c(-1e9, -1e9))
					inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, n - 1] - 1L)
					inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, n - 1] - 1L)
					inserts[[2]] <- c(inserts[[2]], temp[[2]])
					inserts[[4]] <- c(inserts[[4]], temp[[4]])
				}
				n <- n + 1L
			}
			
			n <- 1L
			while (n <= numAnchors) { # align anchor regions
				if (!.Call("firstSeqsGapsEqual",
					pattern,
					subject,
					anchors[1, n], anchors[2, n],
					anchors[3, n], anchors[4, n],
					type,
					max.p,
					max.s,
					PACKAGE="DECIPHER")) {
					temp <- .Call("firstSeqsPosEqual",
						pattern,
						subject,
						anchors[1, n], anchors[2, n],
						anchors[3, n], anchors[4, n],
						type,
						max.p,
						max.s,
						PACKAGE="DECIPHER")
					if (length(temp)==4) {
						inserts[[1]] <- c(inserts[[1]], temp[[1]])
						inserts[[3]] <- c(inserts[[3]], temp[[3]])
					} else { # number of sites differs
						temp <- f(p.profile[, anchors[1, n]:anchors[2, n], drop=FALSE],
							s.profile[, anchors[3, n]:anchors[4, n], drop=FALSE],
							tGaps=c(-1e9, -1e9))
						inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[1, n] - 1L)
						inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[3, n] - 1L)
					}
					inserts[[2]] <- c(inserts[[2]], temp[[2]])
					inserts[[4]] <- c(inserts[[4]], temp[[4]])
				}
				n <- n + 1L
			}
			
			end.p <- anchors[2, numAnchors] == w.p
			end.s <- anchors[4, numAnchors] == w.s
			if (end.p && !end.s) { # need to add gaps after pattern
				inserts[[1]] <- c(inserts[[1]], w.p + 1L)
				inserts[[2]] <- c(inserts[[2]],
					w.s - anchors[4, numAnchors])
			} else if (end.s && !end.p) { # need to add gaps after subject
				inserts[[3]] <- c(inserts[[3]], w.s + 1L)
				inserts[[4]] <- c(inserts[[4]],
					w.p - anchors[2, numAnchors])
			} else if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				anchors[2, numAnchors], w.p,
				anchors[4, numAnchors], w.s,
				max.p,
				max.s,
				PACKAGE="DECIPHER")) { # need to align
				temp <- f(p.profile[, anchors[2, numAnchors]:w.p, drop=FALSE],
					s.profile[, anchors[4, numAnchors]:w.s, drop=FALSE],
					tGaps=c(-1e9, terminalGap[2]))
				inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, numAnchors] - 1L)
				inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, numAnchors] - 1L)
				inserts[[2]] <- c(inserts[[2]], temp[[2]])
				inserts[[4]] <- c(inserts[[4]], temp[[4]])
			} # else don't do anything
		}
	}
	
	ns.p <- names(pattern)
	ns.s <- names(subject)
	if (!is.null(ns.p) && !is.null(ns.s)) {
		ns <- c(ns.p, ns.s)
	} else if (!is.null(ns.p)) {
		ns.s <- rep(NA, length(subject))
	} else if (!is.null(ns.s)) {
		ns.p <- rep(NA, length(pattern))
	}
	ns <- c(ns.p, ns.s)
	
	if (length(inserts[[1]]) > 0) {
		o <- order(inserts[[1]])
		pattern <- .Call("insertGaps",
			pattern,
			as.integer(inserts[[1]][o]),
			as.integer(inserts[[2]][o]),
			type,
			processors,
			PACKAGE="DECIPHER")
	}
	if (length(inserts[[3]]) > 0) {
		o <- order(inserts[[3]])
		subject <- .Call("insertGaps",
			subject,
			as.integer(inserts[[3]][o]),
			as.integer(inserts[[4]][o]),
			type,
			processors,
			PACKAGE="DECIPHER")
	}
	
	result <- .append(pattern, subject)
	names(result) <- ns
	
	return(result)
}
