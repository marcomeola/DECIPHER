DigestDNA <- function(sites,
	myDNAStringSet,
	type="fragments",
	strand="both",
	processors=1) {
	
	# error checking
	if (!is.character(sites))
		stop("sites must be a character vector.")
	if (any(is.na(sites)))
		stop("sites cannot be NA.")
	if (is.character(myDNAStringSet))
		myDNAStringSet <- DNAStringSet(myDNAStringSet)
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (length(myDNAStringSet)==0)
		stop("myDNAStringSet is empty.")
	TYPES <- c("fragments", "positions")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	STRANDS <- c("both", "top", "bottom")
	strand <- pmatch(strand[1], STRANDS)
	if (is.na(strand))
		stop("Invalid strand.")
	if (strand == -1)
		stop("Ambiguous strand.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	# parse sites
	DNA_LOOKUP <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N")
	names(DNA_LOOKUP) <- c("T", "G", "C", "A", "K", "Y", "W", "S", "R", "M", "B", "D", "H", "V", "N")
	cut1 <- cut2 <- integer(length(sites)) # position cut is before
	s <- strsplit(sites, "", fixed=TRUE)
	rc_sites <- character(length(sites))
	for (i in seq_along(s)) {
		site <- s[[i]]
		w <- which(site=="/")
		if (length(w)==1) {
			if (site[length(site)]==")") {
				bot <- as.integer(paste(site[(w + 1):(length(site) - 1)],
					collapse=""))
				if (is.na(bot))
					stop("Improperly formatted site:  ", s[[i]])
				top <- which(site=="(")
				if (length(top) != 1L)
					stop("Improperly formatted site:  ", s[[i]])
				if (!all(site[1:(top - 1)] %in% DNA_LOOKUP))
					stop("Unexpected characters found in site:  ", s[[i]])
				sites[i] <- paste(site[1:(top - 1)],
					collapse="")
				rc_sites[i] <- paste(DNA_LOOKUP[site[(top - 1):1]],
					collapse="")
				top <- as.integer(paste(site[(top + 1):(w - 1)],
					collapse=""))
				if (is.na(top))
					stop("Improperly formatted site:  ", s[[i]])
				cut1[i] <- top + nchar(sites[i]) + 1
				cut2[i] <- -bot + 1
			} else {
				if (all(site[-w] %in% DNA_LOOKUP)) {
					cut1[i] <- cut2[i] <- w
					sites[i] <- paste(site[-w], collapse="")
				} else {
					stop("Unexpected characters found in site:  ", s[[i]])
				}
				rc_sites[i] <- paste(rev(DNA_LOOKUP[site[-w]]),
					collapse="")
			}
		} else if (length(w) > 1) { # multiple cut sites
			stop("Multiple cut sites are not supported:  ", s[[i]])
		} else {
			stop("No cut location(s) found in site:  ", s[[i]])
		}
	}
	
	# search for sites
	p <- sites != rc_sites
	cuts_top <- cuts_bot <- integer()
	ns <- names(myDNAStringSet)
	names(myDNAStringSet) <- seq_along(myDNAStringSet)
	for (i in seq_along(sites)) {
		v_top <- vmatchPattern(sites[i],
			myDNAStringSet,
			fixed="subject")
		v_top <- unlist(v_top)
		cuts_top <- c(cuts_top,
			setNames(start(v_top) + cut1[i] - 1,
				names(v_top)))
		cuts_bot <- c(cuts_bot,
			setNames(end(v_top) - cut2[i] + 2,
				names(v_top)))
		if (p[i]) { # not a palindromic site
			# search for the site's reverse complement
			v_bot <- vmatchPattern(rc_sites[i],
				myDNAStringSet,
				fixed="subject")
			v_bot <- unlist(v_bot)
			cuts_top <- c(cuts_top,
				setNames(start(v_bot) + cut2[i] - 1,
					names(v_bot)))
			cuts_bot <- c(cuts_bot,
				setNames(end(v_bot) - cut1[i] + 2,
					names(v_bot)))
		}
	}
	
	# remove out-of-bounds sites
	ws <- width(myDNAStringSet)
	if (length(cuts_top) > 0) {
		w <- c(which(cuts_top > ws[as.integer(names(cuts_top))]),
			which(cuts_top < 2))
		if (length(w) > 0)
			cuts_top <- cuts_top[-unique(w)]
	}
	if (length(cuts_bot) > 0) {
		w <- c(which(cuts_bot > ws[as.integer(names(cuts_bot))]),
			which(cuts_bot < 2))
		if (length(w) > 0)
			cuts_bot <- cuts_bot[-unique(w)]
	}
	
	# initialize a list of cut positions
	if (strand==1L) {
		value <- list(top=integer(0),
			bottom=integer(0))
	} else if (strand==2L) {
		value <- list(top=integer(0))
	} else { # strand==3L
		value <- list(bottom=integer(0))
	}
	cuts <- lapply(seq_along(myDNAStringSet),
		function(x) {
			value
		})
	names(cuts) <- seq_along(myDNAStringSet)
	
	# record cut positions in top strand
	if (strand==1 || strand==2) {
		m <- match(names(cuts_top),
			seq_along(myDNAStringSet))
		names(cuts_top) <- NULL
		for (i in unique(m)) {
			w <- which(m==i)
			cuts[as.character(i)][[1]][[1]] <- sort(unique(cuts_top[w]))
		}
	}
	
	# record cut positions in bottom strand
	if (strand==1 || strand==3) {
		m <- match(names(cuts_bot),
			seq_along(myDNAStringSet))
		names(cuts_bot) <- NULL
		for (i in unique(m)) {
			w <- which(m==i)
			cuts[as.character(i)][[1]][[ifelse(strand==1, 2, 1)]] <- sort(unique(ws[i] - cuts_bot[w] + 2))
		}
	}
	
	# return cut positions
	if (type==2) {
		if (!is.null(ns))
			names(cuts) <- ns
		return(cuts)
	}
	
	# cut myDNAStringSet into fragments
	if (strand==1 || strand==3)
		rc <- reverseComplement(myDNAStringSet)
	fragments <- list()
	for (i in seq_along(myDNAStringSet)) {
		if (strand==1 || strand==2) {
			cut <- cuts[[i]]$top
			if (length(cut) > 0) {
				top <- extractAt(myDNAStringSet[[i]],
					IRanges(start=c(1, cut),
						end=c(cut - 1, ws[i])))
			} else {
				top <- .Call("subsetXStringSet",
					myDNAStringSet,
					i,
					1L,
					processors,
					PACKAGE="DECIPHER")
			}
		} else {
			top <- .Call("subsetXStringSet",
				myDNAStringSet,
				integer(),
				1L,
				processors,
				PACKAGE="DECIPHER")
		}
		if (strand==1 || strand==3) {
			cut <- cuts[[i]]$bottom
			if (length(cut) > 0) {
				bot <- extractAt(rc[[i]],
					IRanges(start=c(1, cut),
						end=c(cut - 1, ws[i])))
			} else {
				bot <- .Call("subsetXStringSet",
					rc,
					i,
					1L,
					processors,
					PACKAGE="DECIPHER")
			}
		} else {
			bot <- .Call("subsetXStringSet",
				myDNAStringSet,
				integer(),
				1L,
				processors,
				PACKAGE="DECIPHER")
		}
		fragments[[i]] <- setNames(.Call("appendXStringSets",
				top,
				bot,
				1L,
				processors,
				PACKAGE="DECIPHER"),
			c(rep("top", length(top)),
				rep("bottom", length(bot))))
	}
	
	fragments <- relist(do.call(base::c,
			fragments),
		fragments)
	
	if (!is.null(ns)) {
		names(fragments) <- ns
	} else {
		names(fragments) <- seq_along(myDNAStringSet)
	}
	
	return(fragments)
}
