StaggerAlignment <- function(myXStringSet,
	tree=NULL,
	threshold=3,
	fullLength=FALSE,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	type <- switch(class(myXStringSet),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet."))
	if (length(myXStringSet) < 3)
		return(myXStringSet)
	u <- unique(width(myXStringSet))
	if (length(u)!=1)
		stop("Sequences in myXStringSet must be the same width (aligned).")
	if (u < 1) # no changes can be made
		return(myXStringSet)
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold <= 0)
		stop("threshold must be greater than zero.")
	full <- list()
	if (is.list(fullLength)) {
		if (length(fullLength)==1L) {
			full$right <- full$left <- fullLength[[1]]
		} else if (length(fullLength)==2L) {
			if (is.null(names(fullLength))) {
				full$left <- fullLength[[1]]
				full$right <- fullLength[[2]]
			} else if (all(names(fullLength) %in% c("left", "right"))) {
				full <- fullLength
			} else {
				stop("The names of fullLength are not 'left' and 'right'.")
			}
		} else {
			stop("fullLength must be a list with 1 or 2 elements.")
		}
		if (length(full$left)==1L) {
			full$left <- rep(full$left, length(myXStringSet))
		} else if (length(full$left)!=length(myXStringSet)) {
			stop("'left' component of fullLength is not length 1 or the same length as myXStringSet.")
		}
		if (length(full$right)==1L) {
			full$right <- rep(full$right, length(myXStringSet))
		} else if (length(full$right)!=length(myXStringSet)) {
			stop("'right' component of fullLength is not length 1 or the same length as myXStringSet.")
		}
	} else {
		if (!is.logical(fullLength)) {
			stop("fullLength must be a logical.")
		}
		if (length(fullLength)==1L) {
			full$right <- full$left <- rep(fullLength, length(myXStringSet))
		} else if (length(fullLength)==length(myXStringSet)) {
			full$right <- full$left <- fullLength
		} else {
			stop("fullLength is not length 1 or the same length as myXStringSet.")
		}
	}
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
	
	if (is.null(tree)) {
		if (verbose) {
			cat("Calculating distance matrix:\n")
			flush.console()
		}
		
		d <- DistanceMatrix(myXStringSet,
			correction="JC",
			processors=processors,
			verbose=verbose)
		
		if (verbose) {
			cat("Constructing neighbor-joining tree:\n")
			flush.console()
		}
		
		suppressWarnings(tree <- IdClusters(d,
			method="NJ",
			asDendrogram=TRUE,
			processors=processors,
			verbose=verbose))
	} else {
		if (class(tree)!="dendrogram")
			stop("tree must be a dendrogram.")
		if (!all(unlist(tree) %in% seq_along(myXStringSet)))
			stop("tree is incompatible with myXStringSet.")
	}
	
	.assignIndels <- function(x) {
		if (is.leaf(x)) {
			attr(x, "state") <- s[x, pos]
			return(x)
		}
		
		x[[1]] <- .assignIndels(x[[1]])
		x[[2]] <- .assignIndels(x[[2]])
		a1 <- attributes(x[[1]])
		a2 <- attributes(x[[2]])
		
		a <- attributes(x)
		a[["state"]] <- unique(c(a1$state, a2$state))
		
		indels <- c(0L, 0L, 0L) # insertions, deletions, mixed
		# record insertions
		gap1 <- .Call("any", a1$state, PACKAGE="DECIPHER")
		gap2 <- .Call("any", a2$state, PACKAGE="DECIPHER")
		if (!is.na(gap1) &&
			!is.na(gap2) &&
			xor(gap1, gap2)) {
			if (gap2) {
				#attr(x[[1]], "edgePar") <- list(col = "plum")
				attr(x[[1]], "ins") <- TRUE
				indels[1] <- indels[1] + 1L
			} else {
				#attr(x[[2]], "edgePar") <- list(col = "plum")
				attr(x[[2]], "ins") <- TRUE
				indels[1] <- indels[1] + 1L
			}
		}
		
		# record deletions
		gap1 <- .Call("all", a1$state, PACKAGE="DECIPHER")
		gap2 <- .Call("all", a2$state, PACKAGE="DECIPHER")
		if (!is.na(gap1) &&
			!is.na(gap2) &&
			xor(gap1, gap2)) {
			if (gap1) {
				#attr(x[[1]], "edgePar") <- list(col = "green")
				indels[2] <- indels[2] + 1L
			} else {
				#attr(x[[2]], "edgePar") <- list(col = "green")
				indels[2] <- indels[2] + 1L
			}
		}
		
		if (!is.null(a1$indels))
			indels <- indels + a1$indels
		if (!is.null(a2$indels))
			indels <- indels + a2$indels
		
		if (indels[1] != 0L) { # insertions in subtree
			if ((indels[1] + indels[3]) >= (1 + indels[2])) {
				# subtree can be better explained by one insertion
				indels <- c(1L, indels[2], indels[2])
				x <- .Call("clearIns", x, PACKAGE="DECIPHER")
				# mark branch as an insertion
				#a[["edgePar"]] <- list(col = "plum")
				a[["ins"]] <- TRUE
			}
		}
		
		a[["indels"]] <- indels
		attributes(x) <- a
		
		return(x)
	}
	
	.groupIns <- function(x) {
		if (!is.null(attr(x, "ins"))) {
			groups[[length(groups) + 1L]] <<- as.integer(unlist(x))
		} else if (!is.leaf(x) &&
			attr(x, "indels")[1] > 0) {
			# keep descending tree
			.groupIns(x[[1]])
			.groupIns(x[[2]])
		}
		
		return(NULL)
	}
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Staggering insertions and deletions:\n")
		flush.console()
		pBar <- txtProgressBar(style=3, max=100)
		percentComplete <- before <- 0L
	}
	
	s <- matrix(nrow=length(myXStringSet),
		ncol=u)
	t <- TerminalChar(myXStringSet)
	for (i in seq_along(myXStringSet)) {
		x <- .Call("subsetXStringSet",
			myXStringSet,
			i,
			type,
			processors)
		x <- strsplit(as.character(x), "", fixed=T)[[1]]
		s[i,] <- x=="-" | x=="."
		
		# exclude terminal gaps
		if (!full$left[i] && t[i, 1] > 0)
			s[i, 1:t[i, 1]] <- NA
		if (!full$right[i] && t[i, 2] > 0)
			s[i, (u - t[i, 2] + 1):u] <- NA
	}
	
	v <- seq_along(myXStringSet)
	ns <- names(myXStringSet)
	pos <- u # start at end
	changeMade <- FALSE
	while (pos > 0L) {
		y <- .assignIndels(tree)
		
		runLength <- 1L
		while ((pos - runLength) > 0L) {
			s1 <- s[, pos]
			s2 <- s[, pos - runLength]
			w1 <- which(is.na(s1))
			w2 <- which(is.na(s2))
			if (length(w1)==length(w2) &&
				length(w1) < length(s1) &&
				all(w1==w2)) {
				if (length(w1) > 0) {
					if (all(s1[-w1]==s2[-w2])) {
						runLength <- runLength + 1L
					} else {
						break
					}
				} else {
					if (all(s1==s2)) {
						runLength <- runLength + 1L
					} else {
						break
					}
				}
			} else {
				break
			}
		}
		
		if (is.null(attr(y, "ins"))) {
			# position does not exist in root sequence
			indels <- attr(y, "indels")
			if (indels[1] > 1L &&
				(indels[1] + indels[3])/indels[2] < threshold) {
				# stagger insertions
				changeMade <- TRUE
				groups <- list()
				.groupIns(y)
				
				for (i in 2:length(groups)) {
					g <- groups[[i]]
					seqs <- .Call("subsetXStringSet",
						myXStringSet,
						v[-g],
						type,
						processors)
					myXStringSet[-g] <- .Call("insertGaps",
						seqs,
						pos + 1L,
						runLength,
						type,
						processors,
						PACKAGE="DECIPHER")
					seqs <- .Call("subsetXStringSet",
						myXStringSet,
						g,
						type,
						processors)
					myXStringSet[g] <- .Call("insertGaps",
						seqs,
						pos - runLength + 1L,
						runLength,
						type,
						processors,
						PACKAGE="DECIPHER")
				}
			}
		}
		
		if (verbose) {
			percentComplete <- as.integer(100*(1 - pos/u))
			if (percentComplete > before) {
				setTxtProgressBar(pBar, percentComplete)
				before <- percentComplete
			}
		}
		
		pos <- pos - runLength
	}
	
	if (changeMade)
		zzz <- .Call("consolidateGaps",
			myXStringSet,
			type,
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
		cat("\n")
	}
	
	names(myXStringSet) <- ns
	
	return(myXStringSet)
}
