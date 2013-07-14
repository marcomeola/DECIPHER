AlignProfiles <- function(pattern,
	subject,
	p.weight=1,
	s.weight=1,
	perfectMatch=2,
	misMatch=-3,
	gapOpening=-6,
	gapExtension=-4,
	terminalGap=-1) {
	
	# error checking
	if (!is(pattern, "DNAStringSet"))
		stop("pattern must be a DNAStringSet.")
	if (length(pattern) < 1)
		stop("At least one sequence is required in the pattern.")
	if (min(width(pattern)) < 2)
		stop("All sequences in the pattern must be at least two nucleotides long.")
	if (!is(subject, "DNAStringSet"))
		stop("subject must be a DNAStringSet.")
	if (length(subject) < 1)
		stop("At least one sequence is required in the subject.")
	if (min(width(subject)) < 2)
		stop("All sequences in the subject must be at least two nucleotides long.")
	if (!is.numeric(p.weight))
		stop("p.weight must be a numeric.")
	if (length(p.weight)!=1 && length(p.weight)!=length(pattern))
		stop("Length of p.weight must equal one or the length of the pattern.")
	if (length(p.weight)==1)
		p.weight <- rep(1, length(pattern))
	if (!is.numeric(s.weight))
		stop("s.weight must be a numeric.")
	if (length(s.weight)!=1 && length(s.weight)!=length(subject))
		stop("Length of p.weight must equal one or the length of the subject.")
	if (length(s.weight)==1)
		s.weight <- rep(1, length(subject))
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
	
	p.profile <- .Call("consensusProfile", pattern, p.weight, PACKAGE="DECIPHER")
	s.profile <- .Call("consensusProfile", subject, s.weight, PACKAGE="DECIPHER")
	
	lp <- dim(p.profile)[2]
	ls <- dim(s.profile)[2]
	
	t <- .Call("alignProfiles",
		p.profile,
		s.profile,
		perfectMatch,
		misMatch,
		gapOpening,
		gapExtension,
		terminalGap[1],
		terminalGap[2],
		PACKAGE="DECIPHER")
	
	end.p <- t[1]
	end.s <- t[2]
	start.p <- t[3]
	start.s <- t[4]
	count <- t[5] - 4
	t <- t[-1:-5]
	
	p.inserts <- character()
	s.inserts <- character()
	p.ats <- integer()
	s.ats <- integer()
	
	if (start.s < start.p) {
		s.inserts <- paste(rep("-", start.p - 1), collapse="")
		s.ats <- 1
	} else if (start.p < start.s) {
		p.inserts <- paste(rep("-", start.s - 1), collapse="")
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
				p.inserts <- c(p.inserts, paste(rep("-", l), collapse=""))
				p.ats <- c(p.ats, i)
				j <- j + l
			} else {
				s.inserts <- c(s.inserts, paste(rep("-", l), collapse=""))
				s.ats <- c(s.ats, j)
				i <- i + l
			}
			k <- k + l
		}
	}
	
	if (ls - end.s > 0) {
		p.inserts <- c(p.inserts, paste(rep("-", ls - end.s), collapse=""))
		p.ats <- c(p.ats, lp + 1)
	}
	
	if (lp - end.p > 0) {
		s.inserts <- c(s.inserts, paste(rep("-", lp - end.p), collapse=""))
		s.ats <- c(s.ats, ls + 1)
	}
	
	if (length(p.ats) > 0)
		pattern <- replaceAt(pattern, at=list(p.ats), value=p.inserts)
	if (length(s.ats) > 0)
		subject <- replaceAt(subject, at=list(s.ats), value=s.inserts)
	
	return(append(pattern, subject))
}
