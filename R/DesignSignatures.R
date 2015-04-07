.CalculateEfficiencyPCR <- function(primer,
	target,
	temp,
	dGini,
	P,
	deltaGrules,
	taqEfficiency=FALSE,
	maxDistance=.4,
	maxGaps=2,
	align=FALSE,
	processors=NULL) {
	
	# error checking
	if (is(primer, "DNAStringSet"))
		primer <- as.character(primer)
	if (is(target, "DNAStringSet"))
		target <- as.character(target)
	if (!is.character(primer))
		stop("primer must be a character vector.")
	if (!is.character(target))
		stop("target must be a character vector.")
	if (!is.numeric(P))
		stop("P must be a numeric.")
	if (!(P > 0))
		stop("P must be greater than zero.")
	if (!is.numeric(temp))
		stop("temp must be a numeric.")
	if (temp < -273)
		stop("temp must be greater than or equal to absolute zero.")
	
	l <- length(primer)
	if (l==0)
		stop("No primer specified.")
	if (l!=length(target))
		stop("primer is not the same length as target.")
	
	if (align) {
		p <- pairwiseAlignment(primer,
			target,
			type="global",
			gapOpen=-5,
			gapExtension=-5)
		primer <- as.character(pattern(p))
		target <- as.character(subject(p))
	} else {
		if (any(nchar(target) != nchar(primer)))
			stop("primer and target must be aligned (equal length).")
	}
	
	if (taqEfficiency) {
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
		
		seqs2 <- as.character(reverseComplement(DNAStringSet(target)))
		
		# calculate elongation efficiency
		eff_taq <- .Call("terminalMismatch",
			primer,
			seqs2,
			maxDistance,
			maxGaps,
			processors,
			PACKAGE="DECIPHER")
	} else {
		eff_taq <- numeric(l)
		eff_taq[] <- 1
	}
	
	dG <- dGini + .Call("calculateDeltaG",
		primer,
		target,
		deltaGrules,
		PACKAGE="DECIPHER")
	
	RT <- .0019871*(273.15 + temp) # [kcal/mol]
	eff <- P*exp(-dG/RT)/(1 + P*exp(-dG/RT))
	return(eff*eff_taq)
}

DesignSignatures <- function(dbFile,
	tblName="DNA",
	identifier="",
	focusID=NA,
	type="melt",
	resolution=0.5,
	levels=10,
	enzymes=NULL,
	minLength=17,
	maxLength=26,
	maxPermutations=4,
	annealingTemp=64,
	P=4e-7,
	monovalent=70e-3,
	divalent=3e-3,
	dNTPs=8e-04,
	minEfficiency=0.8,
	ampEfficiency=0.5,
	numPrimerSets=100,
	minProductSize=70,
	maxProductSize=400,
	kmerSize=8,
	searchPrimers=500,
	maxDictionary=20000,
	primerDimer=1e-7,
	taqEfficiency=TRUE,
	processors=NULL,
	verbose=TRUE) {
	
	# error checking:
	TYPES <- c("melt", "length", "sequence")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.numeric(monovalent))
		stop("monovalent cations (Na+, K+) concentration must be a numeric.")
	if (!is.numeric(divalent))
		stop("divalent cations (Mg++) concentration must be a numeric.")
	if (!is.numeric(dNTPs))
		stop("dNTPs concentration must be a numeric.")
	ions <- monovalent + 3.3*sqrt(divalent - dNTPs)
	if (ions < .01 || is.nan(ions))
		stop("Sodium equivilent concentration must be at least 0.01M.")
	if (!is.numeric(minLength))
		stop("minLength must be a numeric.")
	if (floor(minLength)!=minLength)
		stop("minLength must be a whole number.")
	if (!is.numeric(levels))
		stop("levels must be a numeric.")
	if (floor(levels)!=levels)
		stop("levels must be a whole number.")
	if (levels < 2)
		stop("levels must be at least two.")
	if (!is.numeric(maxLength))
		stop("maxLength must be a numeric.")
	if (floor(maxLength)!=maxLength)
		stop("maxLength must be a whole number.")
	if (minLength > maxLength)
		stop("minLength must be less or equal to maxLength.")
	if (minLength < 8)
		stop("minLength must be at least 8 nucleotides.")
	if (!is.numeric(annealingTemp))
		stop("annealingTemp must be a numeric.")
	if (annealingTemp < 10)
		stop("annealingTemp must be at least 10 degrees Celsius.")
	if (annealingTemp >= 90)
		stop("annealingTemp must be less than 90 degrees Celsius.")
	if (!is.numeric(P))
		stop("P must be a numeric.")
	if (!(P > 0))
		stop("P must be greater than zero.")
	if (!is.numeric(minEfficiency))
		stop("minEfficiency must be a numeric.")
	if (minEfficiency < 0 || minEfficiency > 1)
		stop("minEfficiency must be between zero and one.")
	if (!is.numeric(ampEfficiency))
		stop("ampEfficiency must be a numeric.")
	if (ampEfficiency < 0 || ampEfficiency > minEfficiency)
		stop("ampEfficiency must be between zero and minEfficiency.")
	if (!is.numeric(numPrimerSets))
		stop("numPrimerSets must be a numeric.")
	if (floor(numPrimerSets)!=numPrimerSets)
		stop("numPrimerSets must be a whole number.")
	if (numPrimerSets <= 0)
		stop("numPrimerSets must be greater than zero.")
	if (!is.numeric(minProductSize))
		stop("minProductSize must be a numeric.")
	if (floor(minProductSize)!=minProductSize)
		stop("minProductSize must be a whole number.")
	if (minProductSize <= 2*maxLength)
		stop("minProductSize must be greater than 2*maxLength.")
	if (!is.numeric(maxProductSize))
		stop("maxProductSize must be a numeric.")
	if (floor(maxProductSize)!=maxProductSize)
		stop("maxProductSize must be a whole number.")
	if (maxProductSize < minProductSize)
		stop("maxProductSize must be greater than or equal to minProductSize.")
	if (!is.numeric(kmerSize))
		stop("kmerSize must be a numeric.")
	if (floor(kmerSize)!=kmerSize)
		stop("kmerSize must be a whole number.")
	if (kmerSize > minLength)
		stop("kmerSize must be less than or equal to minLength.")
	if (kmerSize > 9)
		stop("kmerSize can be at most 9.")
	if (kmerSize < 4)
		stop("kmerSize must be at least 4.")
	if (!is.numeric(maxPermutations))
		stop("maxPermutations must be a numeric.")
	if (floor(maxPermutations)!=maxPermutations)
		stop("maxPermutations must be a whole number.")
	if (maxPermutations <= 0)
		stop("maxPermutations must be greater than zero.")
	if (!is.numeric(searchPrimers))
		stop("searchPrimers must be a numeric.")
	if (floor(searchPrimers)!=searchPrimers)
		stop("searchPrimers must be a whole number.")
	if (searchPrimers <= numPrimerSets)
		stop("searchPrimers must be greater than numPrimerSets.")
	if (!is.numeric(maxDictionary))
		stop("maxDictionary must be a numeric.")
	if (floor(maxDictionary)!= maxDictionary)
		stop("maxDictionary must be a whole number.")
	if (maxDictionary < searchPrimers)
		stop("maxDictionary must be as large as searchPrimers.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(taqEfficiency))
		stop("taqEfficiency must be a logical.")
	if (!is.numeric(primerDimer))
		stop("primerDimer must be a numeric.")
	if (primerDimer <= 0)
		stop("primerDimer must be greater than zero.")
	if (primerDimer > 1)
		stop("primerDimer must be less than or equal to one.")
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
	if (!is.numeric(resolution))
		stop("resolution must be a numeric.")
	if (length(resolution)==1) {
		if (type==1L) { # type is "melt"
			if (resolution < 0.1)
				stop("resolution must be at least 0.1.")
			if (resolution >= 3)
				stop("resolution must be less than 3.")
			# resolution is in degrees Celsius
			maxTemp <- 110 # melt temp is often over-predicted
			resolution <- seq(annealingTemp,
				maxTemp,
				by=resolution)
		} else  {
			if (resolution < 1)
				stop("resolution must be at least 1.")
			if (resolution != floor(resolution))
				stop("resolution must be a whole number.")
			if (type==2L) { # type is "length"
				# resolution specifies the number of bins
				resolution <- seq(minProductSize,
					maxProductSize,
					length.out=resolution + 1)
			} else {
				# resolution specifies the k-mer size
				if (resolution > 6)
					stop("resolution can be at most 6 if type is 'sequence'.")
			}
		}
	} else {
		if (type==3L)
			stop("resolution must be a single number if type is 'sequence'.")
		# resolution specifies boundaries
		if (is.unsorted(resolution))
			stop("resolution must be monotonically increasing.")
		if (type==1L) { # type is "melt"
			maxTemp <- resolution[length(resolution)]
			if (resolution[1] > annealingTemp ||
				maxTemp < annealingTemp)
				stop("The maximum resolution must surpass annealingTemp.")
		} else { # type is "length"
			if (resolution[1] > minProductSize ||
				resolution[length(resolution)] < maxProductSize)
				stop("resolution must span minProductSize to maxProductSize.")
		}
	}
	if (!is.null(enzymes)) {
		if (type==3L)
			stop("enzymes must be NULL if type is 'sequence'.")
		if (!is.character(enzymes))
			stop("enzymes must be a named character vector.")
		if (is.null(names(enzymes)))
			stop("enzymes must be a named character vector.")
		if (any(duplicated(names(enzymes))))
			stop("enzymes must have unique names.")
	}
	
	# initialize database
	driver = dbDriver("SQLite")
	if (is.character(dbFile)) {
		dbConn = dbConnect(driver, dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn = dbFile
		if (!inherits(dbConn,"SQLiteConnection")) 
			stop("'dbFile' must be a character string or SQLiteConnection.")
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	
	searchExpression <- paste("select distinct id from",
		tblName)
	rs <- dbSendQuery(dbConn, searchExpression)
	searchResult <- fetch(rs, n=-1)
	ids <- searchResult$id
	dbClearResult(rs)
	
	if (identifier[1]=="") {
		identifier <- ids
	} else {
		w <- which(!(identifier %in% ids))
		if (length(w) > 0)
			stop("identifier(s) not in dbFile: ",
				paste(identifier[w], collapse=", "))
	}
	
	if (!is.na(focusID[1])) {
		focusID <- match(focusID, identifier)
		if (any(is.na(focusID)))
			stop("focusID not found in identifier.")
	}
	
	# set the seed for repeatable sampling
	last.seed <- .Random.seed
	set.seed(1234)
	on.exit({.Random.seed <- last.seed})
	
	# load dS and dH rules
	data("deltaSrules",
		"deltaHrules",
		envir=environment())
	# apply salt correction
	deltaSrules <- deltaSrules + 0.368*log(ions)/1000
	# calculate free energy parameters
	deltaGrules <- deltaHrules - (273.15 + annealingTemp)*deltaSrules
	# initiation free energy based on GC and AT end pairings
	dSini <- -0.0028*c(0, 1, 2) + 0.0041*c(2, 1, 0) + 0.368*log(ions)/1000
	dHini <- 0.1*c(0, 1, 2) + 2.3*c(2, 1, 0)
	# dGini by terminal base-pairings (left/right)
	dGini <- dHini - (273.15 + annealingTemp)*dSini # AT/AT, GC/AT, GC/GC
	
	# tally the relative fraction of all possible k-mers
	kmers <- numeric(4^kmerSize)
	names(kmers) <- mkAllStrings(c("A", "C", "G", "T"), kmerSize)
	longest <- 0
	if (verbose) {
		time.1 <- Sys.time()
		cat("Tallying ",
			kmerSize,
			"-mers for ",
			length(identifier),
			ifelse(length(identifier)==1L,
				" group",
				" groups"),
			":\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(style=3)
	}
	for (i in seq_along(identifier)) {
		dna <- SearchDB(dbConn,
			tblName=tblName,
			identifier=identifier[i],
			removeGaps="all",
			verbose=FALSE)
		w <- width(dna)
		dna <- unlist(dna)
		if (is.na(focusID[1]) && length(dna) > longest) {
			longest <- length(dna)
			names(longest) <- i
		}
		dna <- replaceAt(dna,
			IRanges(start=cumsum(w) + 1L,
				width=0),
			paste(rep("-", kmerSize),
				collapse=""))
		t <- oligonucleotideFrequency(dna,
			kmerSize,
			as.prob=TRUE)
		m <- match(names(t), names(kmers))
		w <- which(!is.na(m))
		kmers[m[w]] <- kmers[m[w]] + t[w]
		if (verbose)
			setTxtProgressBar(pBar, i/length(identifier))
	}
	
	if (is.na(focusID[1]))
		focusID <- as.integer(names(longest))
	
	# remove k-mers that bind too strongly
	primers <- paste(paste(rep("AA",
				floor(minLength - kmerSize - 1)/2),
			collapse=""),
		ifelse((floor(minLength - kmerSize - 1) %% 2) == 1,
			"T",
			""),
		names(kmers),
		sep="")
	
	.dGinis <- function(primers) {
		left <- substring(primers, 1, 1)
		right <- substring(primers, nchar(primers), nchar(primers))
		GC_term <- (left=="G" | left=="C") + (right=="G" | right=="C")
		
		return(dGini[GC_term + 1])
	}
	
	eff <- .CalculateEfficiencyPCR(primers,
		primers,
		annealingTemp,
		.dGinis(primers),
		P,
		deltaGrules,
		taqEfficiency=FALSE,
		processors=processors)
	w <- which(eff > minEfficiency)
	if (length(w) > 0)
		kmers <- kmers[-w]
	w <- which(kmers >= quantile(kmers, 0.95) & kmers > 0)
	kmers <- names(kmers[w])
	
	# design primers based on the focusID group(s)
	if (verbose) {
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		time.1 <- Sys.time()
		if (is.na(focusID[1])) {
			cat("\nDesigning primer sequences based on the group '",
				identifier[focusID],
				"':\n",
				sep="")
		} else if (length(focusID)==1) {
			cat("\nDesigning primer sequences based on the group '",
				identifier[focusID],
				"':\n",
				sep="")
		} else {
			cat("\nDesigning primer sequences based on ",
				length(focusID),
				" groups:\n",
				sep="")
		}
		flush.console()
		pBar <- txtProgressBar(style=3)
	}
	
	pdict <- PDict(kmers,
		tb.start=1,
		tb.width=kmerSize)
	dna <- DNAStringSet()
	for (i in 1:length(focusID)) {
		dna <- c(dna,
			SearchDB(dbConn,
				tblName=tblName,
				identifier=identifier[focusID[i]],
				removeGaps="all",
				verbose=FALSE))
	}
	w <- width(dna)
	dna <- unlist(dna)
	dna <- replaceAt(dna,
		IRanges(start=cumsum(w) + 1L,
			width=0),
		paste(rep("-", kmerSize),
			collapse=""))
	
	checkBounds <- function(x, start, end) {
		w <- which(x < start | x > end)
		if (length(w) > 0) {
			return(x[-w])
		} else {
			return(x)
		}
	}
	
	m <- matchPDict(pdict, dna)
	starts <- startIndex(m)
	starts <- lapply(starts,
		checkBounds,
		start=1 + maxLength - kmerSize,
		end=length(dna) - maxLength + 1)
	
	hits <- which(unlist(lapply(starts, length)) > 0)
	f.primers <- r.primers <- DNAStringSet()
	for (i in seq_along(hits)) {
		f.tails <- extractAt(dna,
			at=IRanges(start=starts[[hits[i]]] - (maxLength - kmerSize),
				width=maxLength))
		r.tails <- extractAt(dna,
			at=IRanges(start=starts[[hits[i]]],
				width=maxLength))
		
		# determine primer lengths
		a <- alphabetFrequency(f.tails)
		if (dim(a)[1]==1) {
			ambig <- which(sum(a[, 5:18]) > 0)
		} else {
			ambig <- which(rowSums(a[, 5:18]) > 0)
		}
		if (length(ambig) > 0)
			f.tails <- f.tails[-ambig]
		if (length(f.tails) > 0) {
			u <- unique(f.tails)
			m <- match(f.tails, u)
			u <- as.character(u)
			W <- 1:length(u)
			for (j in minLength:maxLength) {
				primers <- substring(u[W], maxLength - j + 1, maxLength)
				if (length(primers)==0)
					next
				eff <- .CalculateEfficiencyPCR(primers,
					primers,
					annealingTemp,
					.dGinis(primers),
					P,
					deltaGrules,
					taqEfficiency=FALSE,
					processors=processors)
				w <- which(eff > minEfficiency)
				if (length(w) > 0) {
					u[W[w]] <- primers[w]
					W <- W[-w]
					if (length(W)==0)
						break
				}
			}
			
			w <- which(nchar(u)==minLength)
			if (length(w) > 0) {
				primers <- substring(u[w], 2, minLength)
				eff <- .CalculateEfficiencyPCR(primers,
					primers,
					annealingTemp,
					.dGinis(primers),
					P,
					deltaGrules,
					taqEfficiency=FALSE,
					processors=processors)
				remove <- which(eff > minEfficiency)
				if (length(remove) > 0)
					W <- c(W, w[remove])
			}
			
			f.tails <- DNAStringSet(u)[m]
			if (length(W) > 0)
				f.tails <- f.tails[-which(m %in% W)]
		}
		
		a <- alphabetFrequency(r.tails)
		if (dim(a)[1]==1) {
			ambig <- which(sum(a[, 5:18]) > 0)
		} else {
			ambig <- which(rowSums(a[, 5:18]) > 0)
		}
		if (length(ambig) > 0)
			r.tails <- r.tails[-ambig]
		if (length(r.tails) > 0) {
			u <- unique(r.tails)
			m <- match(r.tails, u)
			u <- as.character(reverseComplement(u))
			W <- 1:length(u)
			for (j in minLength:maxLength) {
				primers <- substring(u[W], maxLength - j + 1, maxLength)
				if (length(primers)==0)
					next
				eff <- .CalculateEfficiencyPCR(primers,
					primers,
					annealingTemp,
					.dGinis(primers),
					P,
					deltaGrules,
					taqEfficiency=FALSE,
					processors=processors)
				w <- which(eff > minEfficiency)
				if (length(w) > 0) {
					u[W[w]] <- primers[w]
					W <- W[-w]
					if (length(W)==0)
						break
				}
			}
			
			w <- which(nchar(u)==minLength)
			if (length(w) > 0) {
				primers <- substring(u[w], 2, minLength)
				eff <- .CalculateEfficiencyPCR(primers,
					primers,
					annealingTemp,
					.dGinis(primers),
					P,
					deltaGrules,
					taqEfficiency=FALSE,
					processors=processors)
				remove <- which(eff > minEfficiency)
				if (length(remove) > 0)
					W <- c(W, w[remove])
			}
			
			r.tails <- reverseComplement(DNAStringSet(u))[m]
			if (length(W) > 0)
				r.tails <- r.tails[-which(m %in% W)]
		}
		
		# sort primers by 3'-end and tabulate
		t <- table(reverse(f.tails))
		f.tails <- reverse(DNAStringSet(names(t)))
		names(f.tails) <- t
		if (length(f.tails) > 0) {
			o <- order(t, decreasing=TRUE)
			if (length(o) > 10)
				o <- o[1:10]
			f.primers <- c(f.primers, f.tails[o])
			
			# remove repeated primers to
			# prevent mirror primer sets
			w <- which(r.tails %in% f.tails[o])
			if (length(w) > 0)
				r.tails <- r.tails[-w]
		}
		
		t <- table(reverse(r.tails))
		r.tails <- reverse(DNAStringSet(names(t)))
		names(r.tails) <- t
		if (length(r.tails) > 0) {
			o <- order(t, decreasing=TRUE)
			if (length(o) > 10)
				o <- o[1:10]
			r.primers <- c(r.primers, r.tails[o])
		}
		
		if (verbose)	
			setTxtProgressBar(pBar, i/length(hits))
	}
	
	if (length(f.primers) > maxDictionary) {
		o <- order(as.integer(names(f.primers)), decreasing=TRUE)
		f.primers <- f.primers[o[1:maxDictionary]]
	}
	if (length(r.primers) > maxDictionary) {
		o <- order(as.integer(names(r.primers)), decreasing=TRUE)
		r.primers <- r.primers[o[1:maxDictionary]]
	}
	
	# count the number of hits in each group
	if (verbose) {
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		time.1 <- Sys.time()
		cat("\nSelecting the most common primer sequences:\n")
		flush.console()
		pBar <- txtProgressBar(style=3)
	}
	
	f.p <- f.primers
	names(f.p) <- 1:length(f.p)
	f.pdict <- PDict(f.p,
		tb.end=-1,
		tb.width=kmerSize)
	f.w <- width(f.pdict)
	r.p <- r.primers
	names(r.p) <- 1:length(r.p)
	r.pdict <- PDict(r.p,
		tb.start=1,
		tb.width=kmerSize)
	r.w <- width(r.pdict)
	if (maxPermutations > 0) {
		f.a <- f.cm <- array(0,
			dim=c(length(f.pdict),
				5,
				maxLength))
		r.a <- r.cm <- array(0,
			dim=c(length(r.pdict),
				5,
				maxLength))
		for (j in seq_along(f.primers)) {
			f.cm[j,,] <- consensusMatrix(f.primers[j],
				width=maxLength,
				baseOnly=TRUE)
		}
		for (j in seq_along(r.primers)) {
			r.cm[j,,] <- consensusMatrix(r.primers[j],
				width=maxLength,
				baseOnly=TRUE)
		}
	}
	for (i in seq_along(identifier)) {
		dna <- SearchDB(dbConn,
			tblName=tblName,
			identifier=identifier[i],
			removeGaps="all",
			verbose=FALSE)
		w <- width(dna)
		dna <- unlist(dna)
		dna <- replaceAt(dna,
			IRanges(start=cumsum(w) + 1L,
				width=0),
			paste(rep("-", kmerSize),
				collapse=""))
		dna <- xscat(dna, reverseComplement(dna))
		
		p.f <- countPDict(f.pdict, dna)
		m.f <- matchPDict(f.pdict,
			dna,
			min.mismatch=1,
			max.mismatch=1 + ceiling(log(maxPermutations, 2)))
		starts <- startIndex(m.f)
		starts <- lapply(starts,
			checkBounds,
			start=1 + maxLength - kmerSize,
			end=length(dna) - maxLength + 1)
		hits.f <- unlist(lapply(starts, length))
		
		if (maxPermutations > 1) {
			for (j in which(p.f > 0)) {
				f.a[j,,] <- f.a[j,,] + f.cm[j,,]*p.f[j]/(hits.f[j] + p.f[j])
			}
			
			for (j in which(hits.f > 0)) {
				f.hits <- extractAt(dna,
					at=IRanges(start=starts[[j]],
						width=f.w[j]))
				f.a[j,,] <- f.a[j,,] + consensusMatrix(f.hits,
					width=maxLength,
					baseOnly=TRUE)/(hits.f[j] + p.f[j])
			}
		}
		
		p.r <- countPDict(r.pdict, dna)
		m.r <- matchPDict(r.pdict,
			dna,
			min.mismatch=1,
			max.mismatch=1 + ceiling(log(maxPermutations, 2)))
		starts <- startIndex(m.r)
		starts <- lapply(starts,
			checkBounds,
			start=1 + maxLength - kmerSize,
			end=length(dna) - maxLength + 1)
		hits.r <- unlist(lapply(starts, length))
		
		if (maxPermutations > 1) {
			for (j in which(p.r > 0)) {
				r.a[j,,] <- r.a[j,,] + r.cm[j,,]*p.r[j]/(hits.r[j] + p.r[j])
			}
			
			for (j in which(hits.r > 0)) {
				r.hits <- extractAt(dna,
					at=IRanges(start=starts[[j]],
						width=r.w[j]))
				r.a[j,,] <- r.a[j,,] + consensusMatrix(r.hits,
					width=maxLength,
					baseOnly=TRUE)/(hits.r[j] + p.r[j])
			}
		}
		
		# add hits normalized by max achieved
		max_hits <- max(p.f + hits.f, p.r + hits.r)
		names(f.primers) <- as.numeric(names(f.primers)) + (p.f + hits.f)/max_hits
		names(r.primers) <- as.numeric(names(r.primers)) + (p.r + hits.r)/max_hits
		
		if (verbose)
			setTxtProgressBar(pBar, i/(length(identifier) - 1))
	}
	
	# eliminate primers with fewest hits
	if (length(f.primers) > searchPrimers) {
		t <- as.numeric(names(f.primers))
		q <- 1 - searchPrimers/length(f.primers)
		q <- ifelse(q < 0, 0, q)
		f.w <- which(t > quantile(t, q))
		f.w <- c(f.w,
			sample(which(t==quantile(t, q)),
				searchPrimers - length(f.w)))
		f.primers <- f.primers[f.w]
	} else {
		f.w <- seq_along(f.primers)
	}
	if (length(r.primers) > searchPrimers) {
		t <- as.numeric(names(r.primers))
		q <- 1 - searchPrimers/length(r.primers)
		q <- ifelse(q < 0, 0, q)
		r.w <- which(t > quantile(t, q))
		r.w <- c(r.w,
			sample(which(t==quantile(t, q)),
				searchPrimers - length(r.w)))
		r.primers <- r.primers[r.w]
	} else {
		r.w <- seq_along(r.primers)
	}
	
	# add ambiguity to primers to capture more targets
	w.f <- width(f.primers)
	w.r <- width(r.primers)
	if (maxPermutations > 1) {
		# normalize relative frequencies to a max of 1
		f.a <- f.a/length(identifier)
		r.a <- r.a/length(identifier)
		
		for (i in seq_along(f.w)) {
			a <- seq <- f.a[f.w[i], 1:4, 1:w.f[i]]
			o <- apply(a,
				2,
				rank,
				ties.method="first")
			o <- order(o, a, decreasing=TRUE)
			seq[] <- 0
			seq[o[1:w.f[i]]] <- 1
			prev <- seq
			count <- 1L
			a <- a/rep(colSums(a), each=4)
			# incorporate ambiguity up to maxPermutations
			repeat {
				if (a[o[w.f[i] + count]] < 0.01)
					break # not worth incorporating
				
				# incorporate this base
				seq[o[w.f[i] + count]] <- 1
				perms <- prod(colSums(seq))
				
				if (perms > maxPermutations) {
					seq[o[w.f[i] + count]] <- 0
					break # too many permutations
				}
				
				# iterate to next base
				count <- count + 1L
			}
			
			if (any(seq != prev)) {
				# replace primer with consensus sequence
				rownames(seq) <- c("A", "C", "G", "T")
				seq <- seq/rep(colSums(seq), each=4)
				f.primers[[i]] <- DNAString(consensusString(seq,
					threshold=1e-5,
					ambiguityMap=IUPAC_CODE_MAP))
			}
		}
		
		for (i in seq_along(r.w)) {
			a <- seq <- r.a[r.w[i], 1:4, 1:w.r[i]]
			o <- apply(a,
				2,
				rank,
				ties.method="first")
			o <- order(o, a, decreasing=TRUE)
			seq[] <- 0
			seq[o[1:w.r[i]]] <- 1
			prev <- seq
			count <- 1L
			a <- a/rep(colSums(a), each=4)
			repeat {
				if (a[o[w.r[i] + count]] < 0.01)
					break # not worth incorporating
				
				# incorporate this base
				seq[o[w.r[i] + count]] <- 1
				perms <- prod(colSums(seq))
				
				if (perms > maxPermutations) {
					seq[o[w.r[i] + count]] <- 0
					break # too many permutations
				}
				
				# iterate to next base
				count <- count + 1L
			}
			
			if (any(seq != prev)) {
				# replace primer with consensus sequence
				rownames(seq) <- c("A", "C", "G", "T")
				seq <- seq/rep(colSums(seq), each=4)
				r.primers[[i]] <- DNAString(consensusString(seq,
					threshold=1e-5,
					ambiguityMap=IUPAC_CODE_MAP))
			}
		}
	}
	
	.staggeredPrimerDimer <- function(primer1,
		primer2) {
		primer1 <- .Call("expandAmbiguities",
			primer1,
			"T",
			PACKAGE="DECIPHER")
		primer2 <- .Call("expandAmbiguities",
			primer2,
			"T",
			PACKAGE="DECIPHER")
		
		p <- mapply(function(p1, p2) {
				n1 <- nchar(p1)[1]
				n2 <- nchar(p2)[1]
				n <- min(n1, n2)
				
				e <- expand.grid(1:length(p1),
					1:length(p2),
					4:n) # start at 4 base pairs
				
				p1 <- p1[e[, 1]]
				p2 <- p2[e[, 2]]
				
				p1 <- substring(p1,
					n1 - e[, 3] + 1,
					n1)
				p2 <- substring(p2,
					1,
					e[, 3])
				
				list(p1, p2)
			},
			primer1,
			primer2,
			SIMPLIFY=FALSE)
		
		effs <- lapply(p,
			function(x) {
				max(.CalculateEfficiencyPCR(x[[1]],
					x[[2]],
					annealingTemp,
					dGini[1],
					P,
					deltaGrules,
					taqEfficiency=taqEfficiency,
					processors=processors,
					maxDistance=1,
					maxGaps=maxLength))
			})
		
		return(unlist(effs))
	}
	
	effs.f <- .staggeredPrimerDimer(f.primers,
		reverseComplement(f.primers))
	w <- which(effs.f < primerDimer)
	if (length(w) > 0) {
		f.primers <- f.primers[w]
	} else {
		stop("Not enough primers met the specified constraints.")
	}
	effs.r <- .staggeredPrimerDimer(r.primers,
		reverseComplement(r.primers))
	w <- which(effs.r < primerDimer)
	if (length(w) > 0) {
		r.primers <- r.primers[w]
	} else {
		stop("Not enough primers met the specified constraints.")
	}
	
	# find the PCR products for each group by pair of primers
	if (verbose) {
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		time.1 <- Sys.time()
		cat("\nDetermining PCR products from each group:\n")
		flush.console()
		pBar <- txtProgressBar(style=3)
	}
	
	threePrimeWidth <- floor(kmerSize/2)
	f.p <- f.primers
	names(f.p) <- seq_along(f.p)
	f.pdict <- PDict(f.p,
		tb.end=-1,
		tb.width=threePrimeWidth)
	f.pdict2 <- PDict(f.p,
		tb.end=-threePrimeWidth - 1,
		tb.width=threePrimeWidth)
	r.p <- r.primers
	names(r.p) <- seq_along(r.p)
	r.pdict <- PDict(r.p,
		tb.start=1,
		tb.width=threePrimeWidth)
	r.pdict2 <- PDict(r.p,
		tb.start=threePrimeWidth + 1,
		tb.width=threePrimeWidth)
	r.primers.rc <- as.character(reverseComplement(r.p))
	
	if (type==3L) { # type is 'sequence'
		bins <- 4^resolution # possible k-mers
	} else { # type is 'melt' or 'length'
		bins <- length(resolution) - 1 # number of bins
	}
	maxBins <- floor(log(.Machine$integer.max,
		levels))
	ints <- ceiling(bins/maxBins)
	initialize <- function() {
		matrix(as.integer(NA),
			nrow=50000,
			ncol=4 + ints,
			dimnames=list(NULL,
				c("Identifier",
					"Forward",
					"Reverse",
					"Peaks",
					paste("Signature.",
						seq_len(ints),
						sep=""))))
	}
	
	if (type==3L) { # type is 'sequence'
		s <- 1:(4^resolution)
	} else { # type is 'melt' or 'length'
		s <- resolution
		if (type==1L) # define midpoint melt temps for bins
			midpoints <- (s[2:length(s)] + s[1:(length(s) - 1)])/2
	}
	binMat <- matrix(seq_len(maxBins*ints),
		nrow=maxBins)
	.bin <- function(x) {
		l <- .bincode(x,
			breaks=s,
			include.lowest=TRUE)
		t <- tabulate(l, nbins=maxBins*ints)
		b <- round(t*(levels - 1)/max(t))
		l <- levels^(0:(maxBins - 1))
		binMat[] <- b[binMat]
		y <- colSums(binMat*l)
		return(as.integer(y))
	}
	.relist <- function(flesh, skeleton) {
		ind <- 1L
		result <- skeleton
		for (i in which(unlist(lapply(skeleton, length)) > 0)) {
			r <- result[[i]]
			l <- length(r)
			r <- flesh[seq.int(ind, length.out=l)]
			w <- which(r == -1)
			if (length(w) > 0)
				r <- r[-w]
			result[[i]] <- r
			ind <- ind + l
		}
		result
	}
	
	amplicons <- initialize()
	count <- 0
	for (i in seq_along(identifier)) {
		dna <- SearchDB(dbConn,
			tblName=tblName,
			identifier=identifier[i],
			removeGaps="all",
			verbose=FALSE)
		w <- width(dna)
		if (max(width(dna)) < minProductSize)
			next
		
		dna <- unlist(dna)
		dna <- replaceAt(dna,
			IRanges(start=cumsum(w) + 1L,
				width=0),
			paste(rep("-", maxProductSize),
				collapse=""))
		dna <- xscat(dna, reverseComplement(dna)) # both strands
		
		m.f <- matchPDict(f.pdict,
			dna,
			max.mismatch=3 + ceiling(log(maxPermutations, 2)),
			fixed="subject")
		m.f2 <- matchPDict(f.pdict2,
			dna,
			max.mismatch=3 + ceiling(log(maxPermutations, 2)),
			fixed="subject")
		starts <- startIndex(m.f)
		starts2 <- startIndex(m.f2)
		index <- unlist(lapply(starts, length))
		index2 <- unlist(lapply(starts2, length))
		if (all(index==0) && all(index2==0))
			next
		keep <- integer(sum(index, index2))
		pos <- 0L # position in keep
		counter <- 0L # counter along m.f
		counter2 <- sum(index) # counter along m.f2
		for (j in which(index > 0 | index2 > 0)) {
			l <- index[j]
			if (l > 0)
				keep[(pos + 1L):(pos <- pos + l)] <- (counter + 1L):(counter <- counter + l)
			
			# add hits in m.f2 missing from m.f
			w <- which(!(starts2[[j]] %in% starts[[j]]))
			l <- length(w)
			if (l > 0) {
				keep[(pos + 1L):(pos <- pos + l)] <- counter2 + w
				starts[[j]] <- c(starts[[j]],
					starts2[[j]][w])
			}
			counter2 <- counter2 + index2[j]
		}
		w <- which(keep==0)
		if (length(w) > 0)
			keep <- keep[-w]
		index <- unlist(lapply(starts, length)) # reevaluate index
			targets <- extractAllMatches(dna, m.f)
		targets2 <- extractAllMatches(dna, m.f2)
		targets <- suppressWarnings(as.character(targets,
			check.limits=FALSE))
		targets2 <- suppressWarnings(as.character(targets2,
			check.limits=FALSE))
		targets <- c(targets, targets2)[keep]
		
		primers <- rep(as.character(f.p), index)
		w <- which(nchar(targets) > nchar(primers))
		if (length(w) > 0)
			targets[w] <- substring(targets[w],
				nchar(targets[w]) - nchar(primers[w]) + 1,
				nchar(targets[w]))
		w <- which(nchar(targets) < nchar(primers))
		if (length(w) > 0)
			primers[w] <- substring(primers[w],
				nchar(primers[w]) - nchar(targets[w]) + 1,
				nchar(primers[w]))
		eff.f <- .CalculateEfficiencyPCR(primers,
			targets,
			annealingTemp,
			.dGinis(primers),
			P,
			deltaGrules,
			taqEfficiency=taqEfficiency,
			processors=processors)
		w <- which(eff.f < ampEfficiency^2)
		if (length(w) > 0) {
			temp <- unlist(starts)
			temp[w] <- eff.f[w] <- -1
			eff.f <- .relist(eff.f, starts)
			starts <- .relist(temp, starts)
		} else {
			eff.f <- .relist(eff.f, starts)
		}
		
		index <- unlist(lapply(starts, length))
		if (all(index==0))
			next
		
		m.r <- matchPDict(r.pdict,
			dna,
			max.mismatch=3 + ceiling(log(maxPermutations, 2)),
			fixed="subject")
		m.r2 <- matchPDict(r.pdict2,
			dna,
			max.mismatch=3 + ceiling(log(maxPermutations, 2)),
			fixed="subject")
		ends <- endIndex(m.r)
		ends2 <- endIndex(m.r2)
		index <- unlist(lapply(ends, length))
		index2 <- unlist(lapply(ends2, length))
		if (all(index==0) && all(index2==0))
			next
		keep <- integer(sum(index, index2))
		pos <- 0L # position in keep
		counter <- 0L # counter along m.r
		counter2 <- sum(index) # counter along m.r2
		for (j in which(index > 0 | index2 > 0)) {
			l <- index[j]
			if (l > 0)
				keep[(pos + 1L):(pos <- pos + l)] <- (counter + 1L):(counter <- counter + l)
			
			# add hits in m.r2 missing from m.r
			w <- which(!(ends2[[j]] %in% ends[[j]]))
			l <- length(w)
			if (l > 0) {
				keep[(pos + 1L):(pos <- pos + l)] <- counter2 + w
				ends[[j]] <- c(ends[[j]],
					ends2[[j]][w])
			}
			counter2 <- counter2 + index2[j]
		}
		w <- which(keep==0)
		if (length(w) > 0)
			keep <- keep[-w]
		index <- unlist(lapply(ends, length)) # reevaluate index
		targets <- extractAllMatches(dna, m.r)
		targets <- reverseComplement(targets)
		targets2 <- extractAllMatches(dna, m.r2)
		targets2 <- reverseComplement(targets2)
		targets <- suppressWarnings(as.character(targets,
			check.limits=FALSE))
		targets2 <- suppressWarnings(as.character(targets2,
			check.limits=FALSE))
		targets <- c(targets, targets2)[keep]
		
		primers <- rep(r.primers.rc, index)
		w <- which(nchar(targets) > nchar(primers))
		if (length(w) > 0)
			targets[w] <- substring(targets[w],
				1,
				nchar(primers[w]))
		w <- which(nchar(targets) < nchar(primers))
		if (length(w) > 0)
			primers[w] <- substring(primers[w],
				1,
				nchar(targets[w]))
		eff.r <- .CalculateEfficiencyPCR(primers,
			targets,
			annealingTemp,
			.dGinis(primers),
			P,
			deltaGrules,
			taqEfficiency=taqEfficiency,
			processors=processors)
		w <- which(eff.r < ampEfficiency^2)
		if (length(w) > 0) {
			temp <- unlist(ends)
			temp[w] <- -1
			ends <- .relist(temp, ends)
			eff.r <- eff.r[-w]
		}
		
		index <- unlist(lapply(ends, length))
		if (all(index==0))
			next
		index <- rep(1:length(index), index)
		ends <- as.integer(unlist(ends))
		o <- order(ends)
		ends <- ends[o]
		index <- index[o]
		eff.r <- eff.r[o]
		
		for (j in which(unlist(lapply(starts, length)) > 0)) {
			start <- starts[[j]]
			eff <- eff.f[[j]]
			W <- b <- e <- effs <- list()
			for (k in seq_along(start)) {
				w <- .Call("boundedMatches",
					as.integer(ends),
					as.integer(start[k] + minProductSize),
					as.integer(start[k] + maxProductSize - 1L),
					PACKAGE="DECIPHER")
				if (length(w) > 0) {
					# geometric mean > minimum efficiency
					eff.pairs <- sqrt(eff[k]*eff.r[w])
					z <- which(eff.pairs > ampEfficiency)
					if (length(z)==0)
						next
					
					w <- w[z]
					b[[k]] <- rep(start[k], length(w))
					e[[k]] <- ends[w]
					W[[k]] <- w
					effs[[k]] <- eff.pairs[z]
				}
			}
			W <- unlist(W)
			
			if (length(W) > 0) {
				b <- unlist(b)
				e <- unlist(e)
				effs <- unlist(effs)
				if (type==1L) { # melt
					ind <- index[W]
					u <- as.integer(names(table(ind)))
					t <- matrix(nrow=length(u), ncol=ints)
					products <- extractAt(dna,
						at=IRanges(start=b + w.f[j],
							end=e - w.r[ind]))
					products <- xscat(f.p[j],
						products,
						r.primers.rc[ind])
					out <- MeltDNA(products,
						type="melt",
						temps=midpoints,
						ions=ions)
					for (p in seq_along(u)) {
						w <- which(ind==u[p])
						for (k in 1:length(w)) {
							mP <- out[, w[k]]
							if (k==1) {
								mPs <- mP*effs[w[k]]
							} else {
								mPs <- mPs + mP*effs[w[k]]
							}
						}
						mPs <- mPs/sum(effs[w]) # weighted average melt profile
						
						t[p,] <- .bin(rep(midpoints,
								round((levels - 1)*mPs)))
					}
				} else if (type==2L) { # length
					reps <- ceiling(effs*(levels - 1L))
					peaks <- rep(e - b + 1L, reps)
					ind <- rep(index[W], reps)
					t <- tapply(peaks, ind, .bin)
					u <- as.integer(names(t))
					t <- matrix(unlist(t),
						ncol=ints,
						byrow=TRUE)
				} else { # type is 'sequence'
					ind <- index[W]
					u <- names(table(ind))
					t <- matrix(nrow=length(u), ncol=ints)
					products <- extractAt(dna,
						at=IRanges(start=b + w.f[j],
							end=e - w.r[ind]))
					oligos <- oligonucleotideFrequency(products,
						width=resolution,
						with.labels=FALSE)
					oligos <- rowsum(oligos, ind, reorder=FALSE)
					for (p in seq_along(u)) {
						counts <- oligos[u[p],]
						counts <- counts/max(counts)*levels
						w <- which(counts > 0)
						t[p,] <- .bin(rep(w,
							ceiling(counts[w])))
					}
					u <- as.integer(u)
				}
				
				l <- dim(t)[1]
				if ((count + l) > dim(amplicons)[1])
					amplicons <- rbind(amplicons, initialize())
				amplicons[(count + 1):(count + l), "Identifier"] <- i
				amplicons[(count + 1):(count + l), "Forward"] <- j
				amplicons[(count + 1):(count + l), "Reverse"] <- u
				amplicons[(count + 1):(count + l), "Peaks"] <- table(ind)
				amplicons[(count + 1):(count + l), 4 + 1:ints] <- t
				count <- count + l
			}
			if (verbose)
				setTxtProgressBar(pBar, (i + j/length(starts) - 1)/length(identifier))
		}
		if (verbose)	
			setTxtProgressBar(pBar, i/length(identifier))
	}
	
	w <- which(is.na(amplicons[, 1]))
	if (length(w) > 0)
		amplicons <- amplicons[-w, ]
	
	# find the best combinations of forward and reverse primers
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		time.1 <- Sys.time()
		cat("\nScoring signature primer pair combinations:\n")
		flush.console()
		pBar <- txtProgressBar(style=3)
		i <- 0
		previous <- 0
	}
	
	.revBin <- function(y) {
		ints <- integer()
		for (i in 1:length(y)) {
			int <- integer()
			while (y[i] > 0) {
				int <- c(y[i] %% levels, int)
				y[i] <- y[i] %/% levels
			}
			ints <- c(rep(0L,
				ifelse(i < length(y),
					maxBins - length(int),
					bins - (length(int) + length(ints)))),
				int,
				ints)
		}
		return(rev(ints))
	}
	
	pairs <- paste(amplicons[, "Forward"], amplicons[, "Reverse"])
	u <- unique(pairs)
	if (verbose)
		tot <- length(u)
	m <- match(pairs, u)
	o <- order(m)
	amplicons <- amplicons[o,]
	m <- m[o]
	sigs <- prods <- numeric(length(u))
	begin <- 1L
	for (i in 1:length(u)) {
		w <- .Call("multiMatch", m, i, begin, PACKAGE="DECIPHER")
		begin <- w[length(w)]
		
		if (length(w)==1)
			next
		
		sigs[i] <- .Call("intDist",
			amplicons[w, 4 + 1:ints],
			levels,
			bins,
			maxBins,
			length(w),
			length(identifier),
			PACKAGE="DECIPHER")
		prods[i] <- sum(amplicons[w, "Peaks"])
		
		if (verbose && round(i/tot, 2) > previous) {
			previous <- i/tot
			setTxtProgressBar(pBar, previous)
		}
	}
	
	o <- order(sigs, prods, decreasing=TRUE)
	
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		time.1 <- Sys.time()
		cat("\nChoosing optimal forward and reverse pairs:\n")
		flush.console()
		pBar <- txtProgressBar(style=3)
	}
	
	primers <- data.frame(forward_primer=I(character(numPrimerSets)),
		reverse_primer=I(character(numPrimerSets)),
		score=I(numeric(numPrimerSets)),
		coverage=I(numeric(numPrimerSets)),
		products=I(integer(numPrimerSets)))
	
	count <- 0L
	for (i in seq_along(o)) {
		w <- which(m==o[i])
		forward <- f.primers[amplicons[w[1], "Forward"]]
		reverse <- r.primers[amplicons[w[1], "Reverse"]]
		
		# check for probability of dimer artifacts above primerDimer
		pD <- .staggeredPrimerDimer(forward,
			reverse)
		if (pD >= primerDimer)
			next
		
		reverse <- reverseComplement(reverse)
		
		count <- count + 1L
		primers[count, "forward_primer"] <- as.character(forward)
		primers[count, "reverse_primer"] <- as.character(reverse)
		primers[count, "score"] <- sigs[o[i]]
		primers[count, "coverage"] <- length(w)/length(identifier)
		primers[count, "products"] <- prods[o[i]]
		
#		signatures <- list()
#		for (j in 1:length(w)) {
#			signatures[[j]] <- .revBin(amplicons[w[j], 4 + 1:ints])
#		}
#		if (type==3L) { # sequence
#			leadingZeros <- 1L # display all oligos
#		} else { # type is melt or length
#			leadingZeros <- min(unlist(lapply(signatures,
#				function(x) which(x > 0)[1])))
#		}
#		signatures <- unlist(lapply(signatures,
#			function(x) paste(x[leadingZeros:length(x)],
#				collapse=ifelse(levels <= 10, "", " "))))
#		u <- unique(signatures)
#		temp <- ""
#		for (j in seq_along(u)) {
#			W <- which(signatures==u[j])
#			identifiers <- identifier[amplicons[w[W], "Identifier"]]
#			temp <- paste(ifelse(j==1,
#					"",
#					temp),
#				paste(signatures[W[1]],
#					" (",
#					paste(identifiers, collapse=", "),
#					")",
#					sep=""),
#				sep=ifelse(j==1, "", "; "))
#		}
#		w <- which(!(1:length(identifier) %in% amplicons[w, "Identifier"]))
#		if (length(w) > 0)
#			temp <- paste(temp,
#				paste(paste(rep("0",
#						nchar(signatures[1])),
#						collapse=ifelse(levels <= 10, "", " ")),
#					" (",
#					paste(identifier[w], collapse=", "),
#					")",
#					sep=""),
#				sep=ifelse(length(u) > 0, " ;", ""))
#		primers[count, "signatures"] <- temp
		
		setTxtProgressBar(pBar, count/numPrimerSets)
		
		if (count==numPrimerSets)
			break
	}
	
	if (count < numPrimerSets) {
		warning("Not enough primers meet the specified constraints.")
		if (verbose)
			setTxtProgressBar(pBar, 1)
	}
	
	# find the best restriction enzyme to digest the amplicons
	if (length(enzymes) > 0) {
		if (verbose) {
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			time.1 <- Sys.time()
			cat("\nFinding the best restriction enzyme:\n")
			flush.console()
			pBar <- txtProgressBar(style=3)
		}
		
		f.p <- DNAStringSet(primers$forward_primer)
		names(f.p) <- seq_along(f.p)
		f.pdict <- PDict(f.p,
			tb.end=-1,
			tb.width=threePrimeWidth)
		f.pdict2 <- PDict(f.p,
			tb.end=-threePrimeWidth - 1,
			tb.width=threePrimeWidth)
		r.p <- reverseComplement(DNAStringSet(primers$reverse_primer))
		names(r.p) <- seq_along(r.p)
		r.pdict <- PDict(r.p,
			tb.start=1,
			tb.width=threePrimeWidth)
		r.pdict2 <- PDict(r.p,
			tb.start=threePrimeWidth + 1,
			tb.width=threePrimeWidth)
		r.primers.rc <- as.character(reverseComplement(r.p))
		
		initialize <- function() {
			matrix(as.integer(NA),
				nrow=50000,
				ncol=4 + ints,
				dimnames=list(NULL,
					c("Identifier",
						"Primers",
						"Fragments",
						"Enzyme",
						paste("Signature.",
							seq_len(ints),
							sep=""))))
		}
		
		fragments <- initialize()
		count <- 0
		w.f <- width(f.p)
		w.r <- width(r.p)
		for (i in seq_along(identifier)) {
			dna <- SearchDB(dbConn,
				tblName=tblName,
				identifier=identifier[i],
				removeGaps="all",
				verbose=FALSE)
			w <- width(dna)
			if (max(width(dna)) < minProductSize)
				next
			
			dna <- unlist(dna)
			dna <- replaceAt(dna,
				IRanges(start=cumsum(w) + 1L,
					width=0),
				paste(rep("-", maxProductSize),
					collapse=""))
			dna <- xscat(dna, reverseComplement(dna)) # both strands
			
			m.f <- matchPDict(f.pdict,
				dna,
				max.mismatch=3 + ceiling(log(maxPermutations, 2)),
				fixed="subject")
			m.f2 <- matchPDict(f.pdict2,
				dna,
				max.mismatch=3 + ceiling(log(maxPermutations, 2)),
				fixed="subject")
			starts <- startIndex(m.f)
			starts2 <- startIndex(m.f2)
			index <- unlist(lapply(starts, length))
			index2 <- unlist(lapply(starts2, length))
			if (all(index==0) && all(index2==0))
				next
			keep <- integer(sum(index, index2))
			pos <- 0L # position in keep
			counter <- 0L # counter along m.f
			counter2 <- sum(index) # counter along m.f2
			for (j in which(index > 0 | index2 > 0)) {
				l <- index[j]
				if (l > 0)
					keep[(pos + 1L):(pos <- pos + l)] <- (counter + 1L):(counter <- counter + l)
				
				# add hits in m.f2 missing from m.f
				w <- which(!(starts2[[j]] %in% starts[[j]]))
				l <- length(w)
				if (l > 0) {
					keep[(pos + 1L):(pos <- pos + l)] <- counter2 + w
					starts[[j]] <- c(starts[[j]],
						starts2[[j]][w])
				}
				counter2 <- counter2 + index2[j]
			}
			w <- which(keep==0)
			if (length(w) > 0)
				keep <- keep[-w]
			index <- unlist(lapply(starts, length)) # reevaluate index
				targets <- extractAllMatches(dna, m.f)
			targets2 <- extractAllMatches(dna, m.f2)
			targets <- suppressWarnings(as.character(targets,
				check.limits=FALSE))
			targets2 <- suppressWarnings(as.character(targets2,
				check.limits=FALSE))
			targets <- c(targets, targets2)[keep]
			
			p <- rep(as.character(f.p), index)
			w <- which(nchar(targets) > nchar(p))
			if (length(w) > 0)
				targets[w] <- substring(targets[w],
					nchar(targets[w]) - nchar(p[w]) + 1,
					nchar(targets[w]))
			w <- which(nchar(targets) < nchar(p))
			if (length(w) > 0)
				p[w] <- substring(p[w],
					nchar(p[w]) - nchar(targets[w]) + 1,
					nchar(p[w]))
			eff.f <- .CalculateEfficiencyPCR(p,
				targets,
				annealingTemp,
				.dGinis(p),
				P,
				deltaGrules,
				taqEfficiency=taqEfficiency,
				processors=processors)
			w <- which(eff.f < ampEfficiency^2)
			if (length(w) > 0) {
				temp <- unlist(starts)
				temp[w] <- eff.f[w] <- -1
				eff.f <- .relist(eff.f, starts)
				starts <- .relist(temp, starts)
			} else {
				eff.f <- .relist(eff.f, starts)
			}
			
			index <- unlist(lapply(starts, length))
			if (all(index==0))
				next
			
			m.r <- matchPDict(r.pdict,
				dna,
				max.mismatch=3 + ceiling(log(maxPermutations, 2)),
				fixed="subject")
			m.r2 <- matchPDict(r.pdict2,
				dna,
				max.mismatch=3 + ceiling(log(maxPermutations, 2)),
				fixed="subject")
			ends <- endIndex(m.r)
			ends2 <- endIndex(m.r2)
			index <- unlist(lapply(ends, length))
			index2 <- unlist(lapply(ends2, length))
			if (all(index==0) && all(index2==0))
				next
			keep <- integer(sum(index, index2))
			pos <- 0L # position in keep
			counter <- 0L # counter along m.r
			counter2 <- sum(index) # counter along m.r2
			for (j in which(index > 0 | index2 > 0)) {
				l <- index[j]
				if (l > 0)
					keep[(pos + 1L):(pos <- pos + l)] <- (counter + 1L):(counter <- counter + l)
				
				# add hits in m.r2 missing from m.r
				w <- which(!(ends2[[j]] %in% ends[[j]]))
				l <- length(w)
				if (l > 0) {
					keep[(pos + 1L):(pos <- pos + l)] <- counter2 + w
					ends[[j]] <- c(ends[[j]],
						ends2[[j]][w])
				}
				counter2 <- counter2 + index2[j]
			}
			w <- which(keep==0)
			if (length(w) > 0)
				keep <- keep[-w]
			index <- unlist(lapply(ends, length)) # reevaluate index
			targets <- extractAllMatches(dna, m.r)
			targets <- reverseComplement(targets)
			targets2 <- extractAllMatches(dna, m.r2)
			targets2 <- reverseComplement(targets2)
			targets <- suppressWarnings(as.character(targets,
				check.limits=FALSE))
			targets2 <- suppressWarnings(as.character(targets2,
				check.limits=FALSE))
			targets <- c(targets, targets2)[keep]
			
			p <- rep(r.primers.rc, index)
			w <- which(nchar(targets) > nchar(p))
			if (length(w) > 0)
				targets[w] <- substring(targets[w],
					1,
					nchar(p[w]))
			w <- which(nchar(targets) < nchar(p))
			if (length(w) > 0)
				p[w] <- substring(p[w],
					1,
					nchar(targets[w]))
			eff.r <- .CalculateEfficiencyPCR(p,
				targets,
				annealingTemp,
				.dGinis(p),
				P,
				deltaGrules,
				taqEfficiency=taqEfficiency,
				processors=processors)
			w <- which(eff.r < ampEfficiency^2)
			if (length(w) > 0) {
				temp <- unlist(ends)
				temp[w] <- eff.r[w] <- -1
				eff.r <- .relist(eff.r, ends)
				ends <- .relist(temp, ends)
			} else {
				eff.r <- .relist(eff.r, ends)
			}
			
			index <- unlist(lapply(ends, length))
			if (all(index==0))
				next
			
			# loop through pairs of primers
			for (j in seq_along(starts)) {
				start <- starts[[j]]
				end <- ends[[j]]
				eff <- eff.f[[j]]
				effr <- eff.r[[j]]
				W <- b <- e <- effs <- list()
				for (k in seq_along(start)) {
					w <- .Call("boundedMatches",
						as.integer(end),
						as.integer(start[k] + minProductSize),
						as.integer(start[k] + maxProductSize - 1L),
						PACKAGE="DECIPHER")
					if (length(w) > 0) {
						# geometric mean > minimum efficiency
						eff.pairs <- sqrt(eff[k]*effr[w])
						z <- which(eff.pairs > ampEfficiency)
						if (length(z)==0)
							next
						
						w <- w[z]
						b[[k]] <- rep(start[k], length(w))
						e[[k]] <- end[w]
						W[[k]] <- w
						effs[[k]] <- eff.pairs[z]
					}
				}
				W <- unlist(W)
				
				if (length(W) > 0) {
					b <- unlist(b)
					e <- unlist(e)
					effs <- unlist(effs)
					products <- extractAt(dna,
						at=IRanges(start=b + w.f[j],
							end=e - w.r[j]))
					products <- xscat(f.p[j],
						products,
						r.primers.rc[j])
					ws <- width(products)
					for (en in seq_along(enzymes)) {
						cuts <- DigestDNA(enzymes[en],
							products,
							type="positions",
							strand="top",
							processors=processors)
						pieces <- unlist(lapply(cuts,
							function(x)
								length(x[[1]])))
						digested <- which(pieces > 0)
						if (length(digested)==0)
							next # no cut sites
						
						prods <- list()
						for (k in seq_along(digested)) {
							prods[[k]] <- extractAt(products[[k]],
								IRanges(start=c(1, cuts[[k]][[1]]),
									end=c(cuts[[k]][[1]] - 1, ws[k])))
						}
						
						pieces <- pieces[digested] + 1L
						prods <- do.call(base::c,
							prods)
						
						if (type==1L) { # melt
							EFFS <- rep(effs[digested], pieces)
							
							remove <- which(width(prods) < 3)
							if (length(remove) > 0) {
								if (length(remove)==length(prods))
									next
								prods <- prods[-remove]
								EFFS <- EFFS[-remove]
							}
							
							t <- matrix(nrow=1, ncol=ints)
							
							out <- MeltDNA(prods,
								type="melt",
								temps=midpoints,
								ions=ions)
							for (p in seq_len(dim(out)[2])) {
								mP <- out[, p]
								if (p==1) {
									mPs <- mP*EFFS[p]
								} else {
									mPs <- mPs + mP*EFFS[p]
								}
							}
							mPs <- mPs/sum(EFFS) # weighted average melt profile
							
							t[1,] <- .bin(rep(midpoints,
									round((levels - 1)*mPs)))
						} else { # length
							t <- matrix(.bin(width(prods)),
								ncol=ints,
								byrow=TRUE)
						}
						
						l <- dim(t)[1]
						if ((count + l) > dim(fragments)[1])
							fragments <- rbind(fragments, initialize())
						fragments[(count + 1):(count + l), "Identifier"] <- i
						fragments[(count + 1):(count + l), "Primers"] <- j
						fragments[(count + 1):(count + l), "Fragments"] <- sum(pieces)
						fragments[(count + 1):(count + l), "Enzyme"] <- en
						fragments[(count + 1):(count + l), 4 + 1:ints] <- t
						count <- count + l
					}
				}
				
				if (verbose)
					setTxtProgressBar(pBar, (i + j/length(starts) - 1)/length(identifier))
			}
			if (verbose)	
				setTxtProgressBar(pBar, i/length(identifier))
		}
		
		w <- which(is.na(fragments[, 1]))
		if (length(w) > 0)
			fragments <- fragments[-w, ]
		
		primers <- cbind(primers,
			data.frame(enzyme=I(character(dim(primers)[1])),
				digest_score=I(numeric(dim(primers)[1])),
				fragments=I(integer(dim(primers)[1]))))
		
		if (dim(fragments)[1] > 0) {
			pairs <- paste(fragments[, "Primers"],
				fragments[, "Enzyme"])
			u <- unique(pairs)
			m <- match(pairs, u)
			o <- order(m)
			fragments <- fragments[o,]
			m <- m[o]
			sigs <- prods <- numeric(length(u))
			pSets <- e <- integer(length(u))
			begin <- 1L
			for (i in 1:length(u)) {
				w <- .Call("multiMatch", m, i, begin, PACKAGE="DECIPHER")
				begin <- w[length(w)]
				
				pSets[i] <- fragments[w[1], "Primers"]
				e[i] <- fragments[w[1], "Enzyme"]
				prods[i] <- sum(fragments[w, "Fragments"])
				
				if (length(w)==1)
					next # sigs[i] is 0
				
				sigs[i] <- .Call("intDist",
					fragments[w, 4 + 1:ints],
					levels,
					bins,
					maxBins,
					length(w),
					length(identifier),
					PACKAGE="DECIPHER")
			}
			
			primers2 <- primers[pSets,]
			primers2[, "digest_score"] <- sigs
			primers2[, "fragments"] <- prods
			primers2[, "enzyme"] <- names(enzymes)[e]
			
			# combine with unrepresented primer sets
			w <- which(!(seq_len(dim(primers)[1]) %in% pSets))
			primers <- rbind(primers2,
				primers[w,])
		}
		
		o <- order(primers[, "score"] + primers[, "digest_score"],
			primers[, "products"] + primers[, "fragments"],
			decreasing=TRUE)
		primers <- primers[o,]
		
		rownames(primers) <- seq_len(dim(primers)[1])
	}
	
	if (verbose) {
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
	}
	
	return(primers)
}
