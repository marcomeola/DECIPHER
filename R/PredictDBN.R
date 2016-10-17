PredictDBN <- function(myXStringSet,
	type="states",
	minOccupancy=0.5,
	impact=c(1, 1.2, 0.4, -1),
	avgProdCorr=1,
	slope=2,
	shift=1.3,
	threshold=0.5,
	pseudoknots=1,
	weight=1,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!(is(myXStringSet, "DNAStringSet") ||
		is(myXStringSet, "RNAStringSet")))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	l <- length(myXStringSet)
	if (l < 2)
		stop("myXStringSet must contain at least two sequences.")
	u <- unique(width(myXStringSet))
	if (length(u)!=1)
		stop("Sequences in myXStringSet must be the same width (aligned).")
	TYPES <- c("states", "pairs", "scores", "structures")
	if (length(type)==0)
		stop("No type specified.")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (u==0) { # zero width sequences
		if (type==1L) {
			return("")
		} else if (type==2L) {
			return(matrix(nrow=0, ncol=3, dimnames=list(NULL, c("(", ")", "order"))))
		} else if (type==3L) {
			return(matrix(nrow=3, ncol=0, dimnames=list(c(".", "(", ")"), NULL)))
		} else {
			x <- matrix(0, nrow=3, ncol=0, dimnames=list(c(".", "(", ")"), NULL))
			x <- lapply(seq_len(l),
				function(...)
					return(x))
			return(x)
		}
	}
	if (!is.numeric(minOccupancy))
		stop("minOccupancy must be a numeric.")
	if (minOccupancy < 0 || minOccupancy > 1)
		stop("minOccupancy must be between zero and one.")
	if (!is.numeric(impact) || length(impact) != 4)
		stop("impact must be a numeric vector of length four.")
	if (impact[4] > 0)
		stop("The impact of inconsistent pairs (element four) can be at most zero.")
	if (any(impact[1:2] < 0))
		stop("The impact of A/T and G/C base pairs must be at least zero.")
	if (!is.numeric(avgProdCorr))
		stop("avgProdCorr must be a numeric.")
	if (avgProdCorr < 0)
		stop("avgProdCorr must be at least zero.")
	if (!is.numeric(slope))
		stop("slope must be a numeric.")
	if (slope <= 0)
		stop("slope must be greater than zero.")
	if (!is.numeric(shift))
		stop("shift must be a numeric.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold < 0 || threshold > 1)
		stop("threshold must be between zero and one.")
	if (!is.numeric(pseudoknots))
		stop("pseudoknots must be a numeric.")
	if (pseudoknots < 0 || pseudoknots > 3)
		stop("pseudoknots must be between zero and three.")
	if (!is.numeric(weight))
		stop("weight must be a numeric.")
	if (length(weight)!=1 && length(weight)!=length(myXStringSet))
		stop("Length of weight must equal one or the length of the myXStringSet.")
	if (length(weight)==1) {
		weight <- rep(1, length(myXStringSet))
	} else {
		if (!isTRUE(all.equal(1, mean(weight))))
			stop("The mean of weight must be 1.")
	}
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
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
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
	}
	
	ans <- .Call("predictDBN",
		myXStringSet,
		type,
		minOccupancy,
		impact,
		avgProdCorr,
		slope,
		shift,
		weight,
		pseudoknots,
		threshold,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (type==2L) {
		s <- strsplit(ans, "")[[1]]
		n <- sum(s != ".")
		c1 <- c2 <- c3 <- integer(n/2)
		i <- 0L
		parens <- square <- curly <- straight <- integer()
		for (j in seq_along(s)) {
			if (s[j] != ".") {
				if (s[j]=="(") {
					parens <- c(parens, j)
				} else if (s[j]==")") {
					i <- i + 1L
					c1[i] <- parens[length(parens)]
					c2[i] <- j
					length(parens) <- length(parens) - 1L
				} else if (s[j]=="[") {
					square <- c(square, j)
				} else if (s[j]=="]") {
					i <- i + 1L
					c1[i] <- square[length(square)]
					c2[i] <- j
					c3[i] <- 1L
					length(square) <- length(square) - 1L
				} else if (s[j]=="{") {
					curly <- c(curly, j)
				} else if (s[j]=="}") {
					i <- i + 1L
					c1[i] <- curly[length(curly)]
					c2[i] <- j
					c3[i] <- 2L
					length(curly) <- length(curly) - 1L
				} else if (s[j]=="<") {
					straight <- c(straight, j)
				} else if (s[j]==">") {
					i <- i + 1L
					c1[i] <- straight[length(straight)]
					c2[i] <- j
					c3[i] <- 3L
					length(straight) <- length(straight) - 1L
				}
			}
		}
		
		o <- order(c3, c1, c2)
		ans <- matrix(c(c1[o], c2[o], c3[o]),
			ncol=3,
			dimnames=list(NULL, c("(", ")", "order")))
	}
	
	if (type==3L) {
		rownames(ans) <- c(".", "(", ")")
	} else if (type==4L) {
		ans <- lapply(ans,
			function(x) {
				rownames(x) <- c(".", "(", ")")
				return(x)
			})
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(ans)
}
