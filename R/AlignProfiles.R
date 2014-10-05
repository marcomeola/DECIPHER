.align <- function(t, lp, ls) {
	end.p <- t[1]
	end.s <- t[2]
	start.p <- t[3]
	start.s <- t[4]
	count <- t[5] - 4
	t <- t[-1:-5]
	
	p.inserts <- integer()
	s.inserts <- integer()
	p.ats <- integer()
	s.ats <- integer()
	
	if (start.s < start.p) {
		s.inserts <- start.p - 1
		s.ats <- 1
	} else if (start.p < start.s) {
		p.inserts <- start.s - 1
		p.ats <- 1
	}
	
	i <- start.p
	j <- start.s
	
	if (count < length(t)) {
		k <- count + 1
		while (k <= length(t)) {
			l <- length(.Call("multiMatch", t, t[k], as.integer(k), PACKAGE="DECIPHER"))
			if (t[k] == 0) {
				i <- i + l
				j <- j + l
			} else if (t[k] > 0) {
				p.inserts <- c(p.inserts, l)
				p.ats <- c(p.ats, i)
				j <- j + l
			} else {
				s.inserts <- c(s.inserts, l)
				s.ats <- c(s.ats, j)
				i <- i + l
			}
			k <- k + l
		}
	}
	
	if (ls - end.s > 0) {
		p.inserts <- c(p.inserts, ls - end.s)
		p.ats <- c(p.ats, lp + 1)
	}
	
	if (lp - end.p > 0) {
		s.inserts <- c(s.inserts, lp - end.p)
		s.ats <- c(s.ats, ls + 1)
	}
	
	if (length(p.inserts) > 0) {
		u <- unique(p.inserts)
		r <- sapply(u,
			function(x) {
				paste(rep("-", x), collapse="")
			})
		p.inserts <- r[match(p.inserts, u)]
	} else {
		p.inserts <- character()
	}
	if (length(s.inserts) > 0) {
		u <- unique(s.inserts)
		r <- sapply(u,
			function(x) {
				paste(rep("-", x), collapse="")
			})
		s.inserts <- r[match(s.inserts, u)]
	} else {
		s.inserts <- character()
	}
	
	return(list(p.ats, p.inserts, s.ats, s.inserts))
}

AlignProfiles <- function(pattern,
	subject,
	p.weight=1,
	s.weight=1,
	perfectMatch=NULL,
	misMatch=NULL,
	gapOpening=NULL,
	gapExtension=NULL,
	terminalGap=-1,
	restrict=-1000,
	anchor=0.7,
	substitutionMatrix=NULL,
	processors=NULL) {
	
	# error checking
	if (!is(pattern, "DNAStringSet") && !is(pattern, "RNAStringSet") && !is(pattern, "AAStringSet"))
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	if (length(pattern) < 1)
		stop("At least one sequence is required in the pattern.")
	if (!is(subject, "DNAStringSet") && !is(subject, "RNAStringSet") && !is(subject, "AAStringSet"))
		stop("subject must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	if (typeof(pattern) != typeof(subject))
		stop("pattern and subject must be of the same type.")
	if (length(subject) < 1)
		stop("At least one sequence is required in the subject.")
	if (length(unique(width(pattern)))!=1)
		stop("Sequences in pattern must be the same width (aligned).")
	if (length(unique(width(subject)))!=1)
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
	if (is(pattern, "DNAStringSet")) {
		if (is.null(perfectMatch))
			perfectMatch <- 6
		if (is.null(misMatch))
			misMatch <- -3
		if (is.null(gapOpening))
			gapOpening <- -11
		if (is.null(gapExtension))
			gapExtension <- -3
	} else if (is(pattern, "RNAStringSet")) {
		if (is.null(perfectMatch))
			perfectMatch <- 4
		if (is.null(misMatch))
			misMatch <- 0
		if (is.null(gapOpening))
			gapOpening <- -4
		if (is.null(gapExtension))
			gapExtension <- -1
	} else { # AAStringSet
		if (is.null(perfectMatch))
			perfectMatch <- 4
		if (is.null(misMatch))
			misMatch <- 0
		if (is.null(gapOpening))
			gapOpening <- -3
		if (is.null(gapExtension))
			gapExtension <- -5
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
		stop("anchor must be a numeric.")
	if (is.numeric(anchor) && anchor <= 0)
		stop("anchor must be greater than zero.")
	if (is.numeric(anchor) && anchor > 1)
		stop("anchor must be less than or equal to one.")
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
	
	if (is(pattern, "AAStringSet")) {
		if (is.null(substitutionMatrix)) {
			substitutionMatrix <- "BLOSUM62"
		} else if (is.character(substitutionMatrix)) {
			if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
		"PAM30", "PAM40", "PAM70", "PAM120", "PAM250")))
				stop("Invalid substitutionMatrix.")
		}
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.matrix(substitutionMatrix)) {
			if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix)))
		}
		subMatrix <- subMatrix[AAs, AAs]
	} else {
		if (!is.null(substitutionMatrix)) {
			if (is.matrix(substitutionMatrix)) {
				bases <- c("A", "C", "G", "T")
				if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
					any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
					stop("substitutionMatrix is incomplete.")
				substitutionMatrix <- substitutionMatrix[bases, bases]
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		}
	}
	
	f <- function(pattern, subject, tGaps=terminalGap) {
		if (width(pattern)[1]*width(subject)[1] > 2147483647) # maximum when indexing by signed integer
			stop(paste("Alignment larger (",
				width(pattern)[1]*width(subject)[1],
				") than the maximum allowable size (2,147,483,647).",
				sep=""))
		if (is(pattern, "AAStringSet")) {
			p.profile <- .Call("consensusProfileAA",
				pattern,
				p.weight,
				PACKAGE="DECIPHER")
			s.profile <- .Call("consensusProfileAA",
				subject,
				s.weight,
				PACKAGE="DECIPHER")
			
			t <- .Call("alignProfilesAA",
				p.profile,
				s.profile,
				as.numeric(subMatrix),
				perfectMatch,
				misMatch,
				gapOpening,
				gapExtension,
				tGaps[1],
				tGaps[2],
				restrict,
				processors,
				PACKAGE="DECIPHER")
		} else { # DNAStringSet or RNAStringSet
			p.profile <- .Call("consensusProfile",
				pattern,
				p.weight,
				PACKAGE="DECIPHER")
			s.profile <- .Call("consensusProfile",
				subject,
				s.weight,
				PACKAGE="DECIPHER")
			
			t <- .Call("alignProfiles",
				p.profile,
				s.profile,
				as.numeric(substitutionMatrix),
				perfectMatch,
				misMatch,
				gapOpening,
				gapExtension,
				tGaps[1],
				tGaps[2],
				restrict,
				processors,
				PACKAGE="DECIPHER")
		}
		lp <- dim(p.profile)[2]
		ls <- dim(s.profile)[2]
		
		.align(t, lp, ls)
	}
	
	if (is.na(anchor)) { # don't use anchors
		inserts <- f(pattern, subject)
	} else { # use anchors
		w.p <- width(pattern[1])
		w.s <- width(subject[1])
		
		if (is(pattern, "AAStringSet")) {
			wordSize <- 7
		} else {
			wordSize <- 15
		}
		l <- min(length(pattern), length(subject))
		o.p <- order(p.weight, decreasing=TRUE)
		o.s <- order(s.weight, decreasing=TRUE)
		if (is(pattern, "AAStringSet")) {
			num.p <- .Call("enumerateGappedSequenceAA",
				pattern[o.p[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
			num.s <- .Call("enumerateGappedSequenceAA",
				subject[o.s[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
		} else {
			num.p <- .Call("enumerateGappedSequence",
				pattern[o.p[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
			num.s <- .Call("enumerateGappedSequence",
				subject[o.s[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
		}
		
		anchors <- .Call("matchRanges",
			num.p,
			num.s,
			wordSize,
			w.p,
			anchor,
			PACKAGE="DECIPHER")
		
		numAnchors <- dim(anchors)[2]
		if (numAnchors==0) {
			inserts <- f(pattern, subject)
		} else {
			if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				1L, anchors[1, 1],
				1L, anchors[3, 1],
				PACKAGE="DECIPHER")) {
				inserts <- f(subseq(pattern, 1L, anchors[1, 1]),
					subseq(subject, 1L, anchors[3, 1]),
					tGaps=c(terminalGap[1], -100))
			} else {
				inserts <- list(integer(), character(),
					integer(), character())
			}
			
			n <- 2L
			while (n <= numAnchors) { # align regions between anchors
				if (!.Call("firstSeqsEqual",
					pattern,
					subject,
					anchors[2, n - 1], anchors[1, n],
					anchors[4, n - 1], anchors[3, n],
					PACKAGE="DECIPHER")) {
					temp <- f(subseq(pattern, anchors[2, n - 1], anchors[1, n]),
						subseq(subject, anchors[4, n - 1], anchors[3, n]),
						tGaps=c(-100, -100))
					inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, n - 1] - 1L)
					inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, n - 1] - 1L)
					inserts[[2]] <- c(inserts[[2]], temp[[2]])
					inserts[[4]] <- c(inserts[[4]], temp[[4]])
				}
				n <- n + 1L
			}
			
			n <- 1L
			while (n <= numAnchors) { # align anchor regions
				if (!.Call("firstSeqsEqual",
					pattern,
					subject,
					anchors[1, n], anchors[2, n],
					anchors[3, n], anchors[4, n],
					PACKAGE="DECIPHER")) {
					temp <- f(subseq(pattern, anchors[1, n], anchors[2, n]),
						subseq(subject, anchors[3, n], anchors[4, n]),
						tGaps=c(-100, -100))
					inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[1, n] - 1L)
					inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[3, n] - 1L)
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
					paste(rep("-", w.s - anchors[4, numAnchors]), collapse=""))
			} else if (end.s && !end.p) { # need to add gaps after subject
				inserts[[3]] <- c(inserts[[3]], w.s + 1L)
				inserts[[4]] <- c(inserts[[4]],
					paste(rep("-", w.p - anchors[2, numAnchors]), collapse=""))
			} else if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				anchors[2, numAnchors], w.p,
				anchors[4, numAnchors], w.s,
				PACKAGE="DECIPHER")) { # need to align
				temp <- f(subseq(pattern, anchors[2, numAnchors], w.p),
					subseq(subject, anchors[4, numAnchors], w.s),
					tGaps=c(-100, terminalGap[2]))
				inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, numAnchors] - 1L)
				inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, numAnchors] - 1L)
				inserts[[2]] <- c(inserts[[2]], temp[[2]])
				inserts[[4]] <- c(inserts[[4]], temp[[4]])
			} # else don't do anything
		}
	}
	
	if (length(inserts[[1]]) > 0)
		pattern <- replaceAt(pattern, at=inserts[1], value=inserts[[2]])
	if (length(inserts[[3]]) > 0)
		subject <- replaceAt(subject, at=inserts[3], value=inserts[[4]])
	
	return(append(pattern, subject))
}
