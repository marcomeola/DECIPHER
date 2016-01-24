DesignArray <- function(myDNAStringSet,
	maxProbeLength=24,
	minProbeLength=20,
	maxPermutations=4,
	numRecordedMismatches=500,
	numProbes=10,
	start=1,
	end=NULL,
	maxOverlap=5,
	hybridizationFormamide=10,
	minMeltingFormamide=15,
	maxMeltingFormamide=20,
	minScore=-1e12,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(hybridizationFormamide))
		stop("hybridizationFormamide must be a positive numeric.")
	if (hybridizationFormamide <= 0)
		stop("hybridizationFormamide must be a positive numeric.")
	if (!is.numeric(minMeltingFormamide))
		stop("minMeltingFormamide must be a positive numeric.")
	if (minMeltingFormamide <= 0)
		stop("minMeltingFormamide must be a positive numeric.")
	if (!is.numeric(maxMeltingFormamide))
		stop("maxMeltingFormamide must be a positive numeric.")
	if (maxMeltingFormamide <= minMeltingFormamide)
		stop("maxMeltingFormamide must be a greater than minMeltingFormamide.")
	if (!is.numeric(minScore))
		stop("minScore must be a numeric.")
	if (minScore >= 100)
		stop("minScore must be less than the maximum possible score (100).")
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	numF <- length(myDNAStringSet)
	if (!is.numeric(maxProbeLength))
		stop("maxProbeLength must be a positive integer.")
	if (maxProbeLength < 1)
		stop("maxProbeLength must be a positive integer.")
	if ((maxProbeLength < 19) || (maxProbeLength > 26))
		warning("maxProbeLength is not in ideal range from 19 to 26.")
	if (!is.numeric(minProbeLength))
		stop("minProbeLength must be a positive integer.")
	if (minProbeLength < 1)
		stop("minProbeLength must be a positive integer.")
	if ((minProbeLength < 19) || (minProbeLength > 26))
		warning("minProbeLength is not in ideal range from 19 to 26.")
	if (!is.numeric(maxPermutations))
		stop("maxPermutations must be a positive integer.")
	if (maxPermutations < 1)
		stop("maxPermutations must be a positive integer")
	if (!is.numeric(numRecordedMismatches))
		stop("numRecordedMismatches must be a positive integer.")
	if (numRecordedMismatches < 1)
		stop("numRecordedMismatches must be a positive integer")
	if (!is.numeric(numProbes))
		stop("numProbes must be a positive integer.")
	if (numProbes < 1)
		stop("numProbes must be a positive integer")
	if (!is.numeric(start))
		stop("start must be a positive integer.")
	if (start < 1)
		stop("start must be a positive integer")
	if (!is.numeric(maxOverlap))
		stop("maxOverlap must be a positive integer.")
	if (maxOverlap < 0)
		stop("maxOverlap must be a positive integer")
	if (numRecordedMismatches*length(myDNAStringSet) > (2^31 - 1)/2)
		stop("numRecordedMismatches is too large.")
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
	
	if (is.null(end)) {
		end <- min(width(myDNAStringSet))
	} else {
		if (!is.numeric(end))
			stop("end must be a positive integer.")
		if (end > min(width(myDNAStringSet)))
			stop("end is longer than the shortest sequence.")
	}
	if (!is.numeric(start))
		stop("start must be a positive integer.")
	if (start < 1)
		stop("start must be a positive integer")
	if (end < (start + maxProbeLength))
		stop("end must be greater than start.")
	
	if (numF < 2)
		stop("myDNAStringSet must contain more than one sequence.")
	
	# initialize a progress bar
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
	} else {
		pBar <- NULL
	}
	
	maxW <- unique(width(myDNAStringSet))
	if (length(maxW)!=1)
		stop("\nSequences are not aligned.\n")
	
	probes <- .Call("designProbes",
		myDNAStringSet,
		maxProbeLength,
		minProbeLength,
		maxPermutations,
		numRecordedMismatches,
		numProbes,
		start,
		end,
		maxOverlap,
		hybridizationFormamide,
		minMeltingFormamide,
		maxMeltingFormamide,
		minScore,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	w <- which(probes[[1]][,9] != minScore)
	if (is.null(names(myDNAStringSet)))
		names(myDNAStringSet) <- 1:length(myDNAStringSet)
	f <- function(x, myNames) {
		w1 <- which(x[(numRecordedMismatches + 1):(2*numRecordedMismatches)] != -1)
		w1 <- w1[order(x[w1], decreasing=TRUE)]
		if (length(w1) > 0) {
			return(paste(myNames[x[numRecordedMismatches + w1] + 1],
				" (",
				round(x[w1],
					digits=1),
				"%)",
				sep="",
				collapse=", "))
		} else {
			return("")
		}
	}
	mismatches <- unlist(apply(probes[[4]][w,, drop=FALSE], 1, f, names(myDNAStringSet)))
	
	p <- data.frame(name=I(names(myDNAStringSet)[probes[[1]][w,1] + 1]),
		start=I(probes[[1]][w,2]),
		length=I(probes[[1]][w,3]),
		start_aligned=I(probes[[1]][w,6]),
		end_aligned=I(probes[[1]][w,7]),
		permutations=I(probes[[1]][w,4]),
		score=I(probes[[1]][w,5]),
		formamide=I(probes[[1]][w,8]),
		hyb_eff=I(probes[[1]][w,9]),
		target_site=I(probes[[2]][w]),
		probes=I(probes[[3]][w]),
		mismatches=I(mismatches))
	
	if (verbose) {
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(p)
}
