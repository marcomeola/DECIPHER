.primerDimer <- function(primer1,
	primer2,
	temp,
	P,
	ions,
	minOverlap=6,
	processors=NULL) {
	
	RT <- .0019871*(273.15 + temp)
	
	# primer-dimer formation
	dG1 <- as.numeric(system(paste("hybrid-min -n DNA -t",
		temp,
		"-T",
		temp,
		"-N",
		ions,
		"-E -q",
		paste(primer1,
			primer2,
			sep=" ",
			collapse=" ")),
		intern=TRUE))
	K1 <- exp(-dG1/RT)
	
	# primer1 folding
	dG2 <- as.numeric(system(paste("hybrid-ss-min -n DNA -t",
		temp,
		"-T",
		temp,
		"-N",
		ions,
		"-E -q",
			paste(primer1, collapse=" ")),
		intern=TRUE))
	K2 <- exp(-dG2/RT)
	
	# primer2 folding
	dG3 <- as.numeric(system(paste("hybrid-ss-min -n DNA -t",
		temp,
		"-T",
		temp,
		"-N",
		ions,
		"-E -q",
			paste(primer2, collapse=" ")),
		intern=TRUE))
	K3 <- exp(-dG3/RT)
	
	K <- K1/((1 + K2)*(1 + K3))
	
	# Solve P + T <-> PT for PT when P not >> T
	# K*PT^2 - (1 + K*P*T)*PT + K*P*T = 0
	# with simplification P = T (primers in same conc.)
	coeffs <- matrix(c(K*P^2, -1-2*K*P, K), ncol=3)
	roots <- apply(coeffs, 1, polyroot)
	roots <- lapply(roots, function(x) return(Re(x)/P))
	eff_hyb <- unlist(lapply(roots, function (x) {
		w <- which(x <= 1)
		return(x[w[length(w)]])
	}))
	
	primer2 <- reverseComplement(DNAStringSet(primer2))
	primer2 <- sapply(strsplit(toString(primer2), ", ", fixed=TRUE), `[`)
	primer2 <- paste("------------------------------", primer2, "------------------------------", sep="")
	
	mapping <- matrix(0.05, nrow=16, ncol=16)
	dimnames(mapping) <- list(c(names(IUPAC_CODE_MAP), "-"), c(names(IUPAC_CODE_MAP), "-"))
	mapping["A", c("A","M","R","W","V","H","D","N")] <- .7
	mapping["G", c("G","R","S","K","V","D","D","N")] <- 1
	mapping["C", c("C","M","S","Y","V","H","B","N")] <- 1
	mapping["T", c("T","W","Y","K","H","D","B","N")] <- .7
	mapping["G", "T"] <- .2
	mapping["T", "G"] <- .2
	w <- which(mapping==0)
	mapping[w] <- 0
	mapping[,"-"] <- 0.2
	mapping["-",] <- 0.2
	
	p <- pairwiseAlignment(primer1,
		primer2,
		type="global-local",
		gapOpen=-5, gapExt=-5, fuzzyMatrix=mapping)
	
	primer1 <- sapply(strsplit(toString(pattern(p)), ", ", fixed=TRUE), `[`)
	primer2 <- sapply(strsplit(toString(subject(p)), ", ", fixed=TRUE), `[`)
	primer2 <- reverseComplement(DNAStringSet(as.vector(primer2)))
	
	t <- TerminalChar(primer2)
	primer2 <- sapply(strsplit(toString(primer2), ", ", fixed=TRUE), `[`)
	primer1 <- substr(primer1, t[,2] + 1, t[,2] + t[,3])
	primer2 <- substr(primer2, t[,1] + 1, t[,1] + t[,3])
	
	eff1 <- .Call("terminalMismatch", primer1, primer2, .2, 0, processors, PACKAGE="DECIPHER")
	eff2 <- .Call("terminalMismatch", primer2, primer1, .2, 0, processors, PACKAGE="DECIPHER")
	eff_taq <- sqrt(eff1*eff2)
	
	n <- nchar(primer1)
	w <- which(n < minOverlap)
	if (length(w) > 0)
		eff_taq[w] <- 0
	
	return(eff_hyb*eff_taq)
}

DesignPrimers <- function(tiles,
	identifier="",
	start=1,
	end=NULL,
	minLength=17,
	maxLength=26,
	maxPermutations=4,
	minCoverage=0.9,
	minGroupCoverage=0.2,
	annealingTemp=64,
	P=4e-7,
	monovalent=70e-3,
	divalent=3e-3,
	dNTPs=.8e-3,
	minEfficiency=.8,
	worstScore=-Inf,
	numPrimerSets=0,
	minProductSize=75,
	maxProductSize=1200,
	maxSearchSize=1500,
	batchSize=1000,
	maxDistance=.4,
	primerDimer=1e-7,
	ragged5Prime=TRUE,
	taqEfficiency=TRUE,
	induceMismatch=FALSE,
	processors=NULL,
	verbose=TRUE) {
	
	# error checking
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
	if (!is.numeric(maxLength))
		stop("maxLength must be a numeric.")
	if (floor(maxLength)!=maxLength)
		stop("maxLength must be a whole number.")
	if (minLength > maxLength)
		stop("minLength must be less or equal to maxLength.")
	if (minLength < 7)
		stop("minLength must be at least 7 nucleotides.")
	ids <- unique(tiles$id)
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
	if (!is.numeric(worstScore))
		stop("worstScore must be a numeric.")
	if (worstScore > 0)
		stop("worstScore must be less than or equal to zero.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize)!=batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
	if (!is.numeric(numPrimerSets))
		stop("numPrimerSets must be a numeric.")
	if (floor(numPrimerSets)!=numPrimerSets)
		stop("numPrimerSets must be a whole number.")
	if (numPrimerSets < 0)
		stop("numPrimerSets must be greater than or equal to zero.")
	if (!is.numeric(minProductSize))
		stop("minProductSize must be a numeric.")
	if (floor(minProductSize)!=minProductSize)
		stop("minProductSize must be a whole number.")
	if (minProductSize <= 0)
		stop("minProductSize must be greater than zero.")
	if (!is.numeric(maxProductSize))
		stop("maxProductSize must be a numeric.")
	if (floor(maxProductSize)!=maxProductSize)
		stop("maxProductSize must be a whole number.")
	if (maxProductSize <= 0)
		stop("maxProductSize must be greater than zero.")
	if (!is.numeric(maxSearchSize))
		stop("maxSearchSize must be a numeric.")
	if (floor(maxSearchSize)!=maxSearchSize)
		stop("maxSearchSize must be a whole number.")
	if (maxSearchSize <= 0)
		stop("maxSearchSize must be greater than zero.")
	if (!is.numeric(maxPermutations))
		stop("maxPermutations must be a numeric.")
	if (floor(maxPermutations)!=maxPermutations)
		stop("maxPermutations must be a whole number.")
	if (maxPermutations <= 0)
		stop("maxPermutations must be greater than zero.")
	if (!is.numeric(maxDistance))
		stop("maxDistance must be a numeric.")
	if (maxDistance < 0)
		stop("maxDistance must be greater or equal to zero.")
	if (maxDistance > 1)
		stop("maxDistance must be less than or equal to one.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(induceMismatch) && !is.logical(induceMismatch))
		stop("induceMismatch must be a logical or integer between 2 and 6.")
	if (is.numeric(induceMismatch)) {
		if (floor(induceMismatch)!=induceMismatch)
			stop("induceMismatch must be a whole number.")
		if (induceMismatch < 2 || induceMismatch > 6)
			stop("induceMismatch must be between 2 and 6.")
		pos <- induceMismatch
		induceMismatch <- TRUE
	} else {
		pos <- 6
	}
	if (!is.logical(ragged5Prime))
		stop("ragged5Prime must be a logical.")
	if (!is.logical(taqEfficiency))
		stop("taqEfficiency must be a logical.")
	if (!is.numeric(start))
		stop("start must be a numeric.")
	if (!is.numeric(end) && !is.null(end))
		stop("end must be a numeric or NULL.")
	if (start < 1)
		stop("start must be greater than zero.")
	if (floor(start)!=start)
		stop("start must be a whole number.")
	if (!is.numeric(minCoverage))
		stop("minCoverage must be a numeric.")
	if (minCoverage > 1 || minCoverage < 0)
		stop("minCoverage must be between zero and one.")
	if (!is.numeric(minGroupCoverage))
		stop("minGroupCoverage must be a numeric.")
	if (minGroupCoverage > 1 || minGroupCoverage < 0)
		stop("minGroupCoverage must be between zero and one.")
	if (!is.numeric(primerDimer))
		stop("primerDimer must be a numeric.")
	if (primerDimer <= 0)
		stop("primerDimer must be greater than zero.")
	if (primerDimer > 1)
		stop("primerDimer must be less than or equal to one.")
	if (!is.null(end)) {
		if (end < start)
			stop("end must be greater than start.")
		if (floor(end)!=end)
			stop("end must be a whole number.")
	}
	if (identifier[1]=="") {
		identifier <- ids
	} else {
		w <- which(!(identifier %in% ids))
		if (length(w) > 0)
			stop("identifier not in tiles: ",
				paste(identifier[w], collapse=", "))
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
	
	primers_all <- data.frame()
	
	if (numPrimerSets > 0) {
		l <- length(identifier)
		pSets <- data.frame(identifier=I(rep(identifier, each=numPrimerSets)),
			start_forward=I(integer(l*numPrimerSets)),
			start_reverse=I(integer(l*numPrimerSets)),
			product_size=I(integer(l*numPrimerSets)),
			start_aligned_forward=I(integer(l*numPrimerSets)),
			start_aligned_reverse=I(integer(l*numPrimerSets)),
			permutations_forward=I(integer(l*numPrimerSets)),
			permutations_reverse=I(integer(l*numPrimerSets)),
			score_forward=I(numeric(l*numPrimerSets)),
			score_reverse=I(numeric(l*numPrimerSets)),
			score_set=I(numeric(l*numPrimerSets)),
			forward_primer=I(matrix(character(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations)),
			reverse_primer=I(matrix(character(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations)),
			forward_efficiency=I(matrix(numeric(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations)),
			reverse_efficiency=I(matrix(numeric(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations)),
			forward_coverage=I(matrix(numeric(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations)),
			reverse_coverage=I(matrix(numeric(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations)),
			mismatches_forward=I(character(l*numPrimerSets)),
			mismatches_reverse=I(character(l*numPrimerSets)),
			mismatches_set=I(character(l*numPrimerSets)))
		if (induceMismatch) {
			pSets$forward_MM <- I(matrix(character(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations))
			pSets$forward_efficiency_MM <- I(matrix(numeric(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations))
			pSets$reverse_MM <- I(matrix(character(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations))
			pSets$reverse_efficiency_MM <- I(matrix(numeric(),
				nrow=l*numPrimerSets,
				ncol=maxPermutations))
		}
	}
	l_ids <- length(ids)
	for (id in identifier) {
		if (verbose) {
			time.1 <- Sys.time()
			cat("\n", id, sep="")
			flush.console()
		}
		
		if (is.null(end)) {
			w <- which(tiles$id==id &
				tiles$misprime==0 &
				tiles$start_aligned >= start)
		} else {
			w <- which(tiles$id==id &
				tiles$misprime==0 &
				tiles$start_aligned >= start &
				tiles$end_aligned <= end)
		}
		if (length(w)==0) {
			warning("Skipped because no target sites met the specified constraints: ", id)
			next
		}
		
		target <- tiles[w,]
		
		uw <- unique(target$width)
		if (length(uw) > 1) {
			warning("Skipped due to multiple sequence widths: ", id)
			next
		}
		
		max_w <- max(nchar(target$target_site))
		if (max_w < maxLength) {
			warning("Skipped because maxLength is greater than width of target sites: ", id)
			next
		}
		
		u <- unique(target$start_aligned)
		l <- length(u)
		primers <- data.frame(identifier=I(rep(id,l)),
			start_forward=I(integer(l)),
			start_reverse=I(integer(l)),
			start_aligned_forward=I(integer(l)),
			start_aligned_reverse=I(integer(l)),
			permutations_forward=I(integer(l)),
			permutations_reverse=I(integer(l)),
			score_forward=I(numeric(l)),
			score_reverse=I(numeric(l)),
			forward_primer=I(matrix(character(),
				nrow=l,
				ncol=maxPermutations)),
			reverse_primer=I(matrix(character(),
				nrow=l,
				ncol=maxPermutations)),
			forward_efficiency=I(matrix(numeric(),
				nrow=l,
				ncol=maxPermutations)),
			reverse_efficiency=I(matrix(numeric(),
				nrow=l,
				ncol=maxPermutations)),
			forward_coverage=I(matrix(numeric(),
				nrow=l,
				ncol=maxPermutations)),
			reverse_coverage=I(matrix(numeric(),
				nrow=l,
				ncol=maxPermutations)),
			mismatches_forward=I(character(l)),
			mismatches_reverse=I(character(l)))
		
		targets_F <- array("",
			c(l, maxPermutations, maxLength - minLength + 1),
			dimnames=list(position=NULL,permutation=NULL,length=NULL))
		targets_R <- array("",
			c(l, maxPermutations, maxLength - minLength + 1),
			dimnames=list(position=NULL,permutation=NULL,length=NULL))
		
		begin <- 1L
		for (i in 1:l) {
			# find all target_sites
			#w <- which(target$start==u[i])
			w <- .Call("multiMatch", target$start_aligned, u[i], begin, PACKAGE="DECIPHER")
			begin <- w[1]
			target_site <- target$target_site[w]
			perms <- length(target_site)
			
			if (perms > maxPermutations)
				next
			
			if (minCoverage > sum(target$coverage[w]))
				next
			
			if (minGroupCoverage > sum(target$groupCoverage[w]))
				next
			
			primers$permutations_forward[i] <-
				primers$permutations_reverse[i] <- perms
			primers$forward_coverage[i, 1:perms] <-
				primers$reverse_coverage[i, 1:perms] <- target$coverage[w]
			
			# determine optimal length for priming
			for (j in 1:perms) {
				n <- nchar(target_site[j])
				for (p in minLength:ifelse(n > maxLength, maxLength, n)) {
					if (ragged5Prime) {
						targets_F[i, j, p - minLength + 1] <- substr(target_site[j], n - p + 1, n)
						targets_R[i, j, p - minLength + 1] <- substr(target_site[j], 1, p)
					} else {
						targets_F[i, j, p - minLength + 1] <- substr(target_site[j], 1, p)
						targets_R[i, j, p - minLength + 1] <- substr(target_site[j], n - p + 1, n)
					}
				}
			}
			primers$start_aligned_forward[i] <- u[i]
			primers$start_aligned_reverse[i] <- target$end_aligned[w[1]]
			
			primers$start_forward[i] <- target$start[w[1]]
			primers$start_reverse[i] <- target$end[w[1]]
		}
		
		w <- which(primers$permutations_forward==0)
		if (length(w)==dim(primers)[1]) {
			warning("Skipped because all target sites have too many permutations: ", id)
			next
		}
		if (length(w) > 0) {
			primers <- primers[-w,]
			targets_F <- targets_F[-w,,]
			dim(targets_F) <- c(l - length(w), maxPermutations, maxLength - minLength + 1)
			targets_R <- targets_R[-w,,]
			dim(targets_R) <- c(l - length(w), maxPermutations, maxLength - minLength + 1)
			l <- dim(primers)[1]
		}
		
		# Choose optimal length for each forward primer
		w <- integer()
		empty <- integer()
		for (j in 1:(maxLength - minLength + 1)) {
			if (length(empty) > 0) {
				empty <- c((1:(l*maxPermutations))[-empty][w],
					which(targets_F[,,j]==""),
					empty)
			} else {
				empty <- c((1:(l*maxPermutations))[w],
					which(targets_F[,,j]==""))
			}
			t <- targets_F[,,j]
			if (length(empty) > 0) {
				t <- t[-empty]
			} else {
				t <- as.vector(t)
			}
			if (length(t)==0)
				break
			eff <- CalculateEfficiencyPCR(t,
				strsplit(toString(reverseComplement(DNAStringSet(t))),
					", ",
					fixed=TRUE)[[1]],
				annealingTemp,
				P,
				ions,
				batchSize,
				taqEfficiency=taqEfficiency,
				maxDistance,
				processors=processors)
			w <- which(eff > minEfficiency)
			if (length(empty) > 0) {
				primers$forward_primer[-empty][w] <- targets_F[,,j][-empty][w]
				primers$forward_efficiency[-empty][w] <- eff[w]
			} else {
				primers$forward_primer[w] <- targets_F[,,j][w]
				primers$forward_efficiency[w] <- eff[w]
			}
		}
		
		# Choose optimal length for each reverse primer
		w <- integer()
		empty <- integer()
		for (j in 1:(maxLength - minLength + 1)) {
			if (length(empty) > 0) {
				empty <- c((1:(l*maxPermutations))[-empty][w],
					which(targets_R[,,j]==""),
					empty)
			} else {
				empty <- c((1:(l*maxPermutations))[w],
					which(targets_R[,,j]==""))
			}
			t <- targets_R[,,j]
			if (length(empty) > 0) {
				t <- t[-empty]
			} else {
				t <- as.vector(t)
			}
			if (length(t)==0)
				break
			eff <- CalculateEfficiencyPCR(strsplit(toString(reverseComplement(DNAStringSet(t))),
					", ",
					fixed=TRUE)[[1]],
				t,
				annealingTemp,
				P,
				ions,
				batchSize,
				taqEfficiency=taqEfficiency,
				maxDistance,
				processors=processors)
			w <- which(eff > minEfficiency)
			if (length(empty) > 0) {
				primers$reverse_primer[-empty][w] <- strsplit(toString(reverseComplement(DNAStringSet(targets_R[,,j][-empty][w]))),
					", ",
					fixed=TRUE)[[1]]
				primers$reverse_efficiency[-empty][w] <- eff[w]
			} else {
				primers$reverse_primer[w] <- strsplit(toString(reverseComplement(DNAStringSet(targets_R[,,j][w]))),
					", ",
					fixed=TRUE)[[1]]
				primers$reverse_efficiency[w] <- eff[w]
			}
		}
		
		w <- integer()
		for (i in 1:l) {
			if (primers$start_aligned_forward[i]==0) {
				w <- c(w,i) # empty primer
				next
			}
			p <- primers$permutations_forward[i]
			pf <- length(which(primers$forward_primer[i,]!=""))
			pr <- length(which(primers$reverse_primer[i,]!=""))
			if ((p != pr) && (p != pf)) {
				w <- c(w,i)
			} else {
				if (p != pf) {
					primers$permutations_forward[i] <- 0
					primers$forward_primer[i,] <- NA
					primers$forward_efficiency[i,] <- NA
					primers$score_forward[i] <- -Inf
				} else {
					df <- duplicated(primers$forward_primer[i,], incomparables=NA)
					l <- length(which(df))
					if (l > 0) {
						for (j in which(df)) {
							w1 <- which(primers$forward_primer[i,1:(j - 1)]==primers$forward_primer[i,j])[1]
							primers$forward_coverage[i,w1] <- primers$forward_coverage[i,w1] + primers$forward_coverage[i,j]
							primers$forward_coverage[i,j] <- NA
						}
						primers$permutations_forward[i] <- primers$permutations_forward[i] - l
						primers$forward_primer[i,] <- c(primers$forward_primer[i,!df],
							rep(NA, l))
						primers$forward_efficiency[i,] <- c(primers$forward_efficiency[i,!df],
							rep(NA, l))
						primers$forward_coverage[i,] <- c(primers$forward_coverage[i,!df],
							rep(NA, l))
					}
					if (ragged5Prime) # correct starting position
						primers$start_forward[i] <- primers$start_forward[i] + max_w - max(nchar(primers$forward_primer[i,]))
				}
				if (p != pr) {
					primers$permutations_reverse[i] <- 0
					primers$reverse_primer[i,] <- NA
					primers$reverse_efficiency[i,] <- NA
					primers$score_reverse[i] <- -Inf
				} else {
					dr <- duplicated(primers$reverse_primer[i,], incomparables=NA)
					l <- length(which(dr))
					if (l > 0) {
						for (j in which(dr)) {
							w1 <- which(primers$reverse_primer[i,1:(j - 1)]==primers$reverse_primer[i,j])[1]
							primers$reverse_coverage[i,w1] <- primers$reverse_coverage[i,w1] + primers$reverse_coverage[i,j]
							primers$reverse_coverage[i,j] <- NA
						}
						primers$permutations_reverse[i] <- primers$permutations_reverse[i] - l
						primers$reverse_primer[i,] <- c(primers$reverse_primer[i,!dr],
							rep(NA, l))
						primers$reverse_efficiency[i,] <- c(primers$reverse_efficiency[i,!dr],
							rep(NA, l))
						primers$reverse_coverage[i,] <- c(primers$reverse_coverage[i,!dr],
							rep(NA, l))
					}
					if (ragged5Prime) # correct starting position
						primers$start_reverse[i] <- primers$start_reverse[i] - max_w + max(nchar(primers$reverse_primer[i,]))
				}
			}
		}
		
		if (length(w) > 0)
			primers <- primers[-w,]
		d <- dim(primers)[1]
		if (d==0) {
			warning("No primers met the specified constraints: ",id)
			next
		}
		if (verbose) {
			cat(" (", d, "):\n", sep="")
			flush.console()
			pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
		}
		
		if (induceMismatch) {
			empty <- which(is.na(primers$forward_primer))
			if (length(empty)==0) {
				n <- nchar(primers$forward_primer)
				nt <- strsplit(toString(reverseComplement(DNAStringSet(primers$forward_primer))),
					", ",
					fixed=TRUE)[[1]]
				X <- substr(primers$forward_primer, n - pos + 1, n - pos + 1)
			} else {
				n <- nchar(primers$forward_primer[-empty])
				nt <- strsplit(toString(reverseComplement(DNAStringSet(primers$forward_primer[-empty]))),
					", ",
					fixed=TRUE)[[1]]
				X <- substr(primers$forward_primer[-empty], n - pos + 1, n - pos + 1)
			}
			primers$forward_MM <- I(matrix(nrow=dim(primers$forward_primer)[1],
				ncol=dim(primers$forward_primer)[2]))
			primers$forward_efficiency_MM <- I(matrix(nrow=dim(primers$forward_primer)[1],
				ncol=dim(primers$forward_primer)[2]))
			for (i in c("A", "C", "T", "G")) {
				w <- which(X != i)
				if (length(w)==0)
					next
				if (length(empty)==0) {
					t <- primers$forward_primer[w]
				} else {
					t <- primers$forward_primer[-empty][w]
				}
				x <- X[w]
				substr(t, n[w] - pos + 1, n[w] - pos + 1) <- i
				eff_PM <- CalculateEfficiencyPCR(t,
					strsplit(toString(reverseComplement(DNAStringSet(t))),
						", ",
						fixed=TRUE)[[1]],
					annealingTemp,
					P,
					ions,
					batchSize,
					taqEfficiency=taqEfficiency,
					maxDistance,
					processors=processors)
				W <- which(eff_PM >= minEfficiency)
				w <- w[W]
				if (length(w)==0)
					next
				t <- t[W]
				x <- x[W]
				eff_PM <- eff_PM[W]
				eff_MM <- CalculateEfficiencyPCR(t,
					nt[w],
					annealingTemp,
					P,
					ions,
					batchSize,
					taqEfficiency=taqEfficiency,
					maxDistance,
					processors=processors)
				W <- which(is.na(primers$forward_MM[w]))
				nnas <- which(!is.na(primers$forward_MM[w]))
				if (length(nnas) > 0)
					W <- sort(c(W, nnas[which(eff_MM[nnas] > primers$forward_efficiency_MM[w[nnas]])]))
				W <- W[which(eff_MM[W] >= 0.1)] # must have 10% MM efficiency
				if (length(W) > 0) {
					if (length(empty)==0) {
						primers$forward_efficiency_MM[w[W]] <- eff_MM[W]
						primers$forward_efficiency[w[W]] <- eff_PM[W]
						primers$forward_MM[w[W]] <- paste(i,
							"/",
							strsplit(toString(reverseComplement(DNAStringSet(x[W]))),
								", ",
								fixed=TRUE)[[1]],
							sep="")
						primers$forward_primer[w[W]] <- t[W]
					} else {
						primers$forward_efficiency_MM[-empty][w[W]] <- eff_MM[W]
						primers$forward_efficiency[-empty][w[W]] <- eff_PM[W]
						primers$forward_MM[-empty][w[W]] <- paste(i,
							"/",
							strsplit(toString(reverseComplement(DNAStringSet(x[W]))),
								", ",
								fixed=TRUE)[[1]],
							sep="")
						primers$forward_primer[-empty][w[W]] <- t[W]
					}
				}
			}
			
			empty <- which(is.na(primers$reverse_primer))
			if (length(empty)==0) {
				n <- nchar(primers$reverse_primer)
				nt <- strsplit(toString(reverseComplement(DNAStringSet(primers$reverse_primer))),
					", ",
					fixed=TRUE)[[1]]
				X <- substr(primers$reverse_primer, n - pos + 1, n - pos + 1)
			} else {
				n <- nchar(primers$reverse_primer[-empty])
				nt <- strsplit(toString(reverseComplement(DNAStringSet(primers$reverse_primer[-empty]))),
					", ",
					fixed=TRUE)[[1]]
				X <- substr(primers$reverse_primer[-empty], n - pos + 1, n - pos + 1)
			}
			primers$reverse_MM <- I(matrix(nrow=dim(primers$reverse_primer)[1],
				ncol=dim(primers$reverse_primer)[2]))
			primers$reverse_efficiency_MM <- I(matrix(nrow=dim(primers$reverse_primer)[1],
				ncol=dim(primers$reverse_primer)[2]))
			for (i in c("A", "C", "T", "G")) {
				w <- which(X != i)
				if (length(w)==0)
					next
				if (length(empty)==0) {
					t <- primers$reverse_primer[w]
				} else {
					t <- primers$reverse_primer[-empty][w]
				}
				x <- X[w]
				substr(t, n[w] - pos + 1, n[w] - pos + 1) <- i
				eff_PM <- CalculateEfficiencyPCR(t,
					strsplit(toString(reverseComplement(DNAStringSet(t))),
						", ",
						fixed=TRUE)[[1]],
					annealingTemp,
					P,
					ions,
					batchSize,
					taqEfficiency=taqEfficiency,
					maxDistance,
					processors=processors)
				W <- which(eff_PM >= minEfficiency)
				w <- w[W]
				if (length(w)==0)
					next
				t <- t[W]
				x <- x[W]
				eff_PM <- eff_PM[W]
				eff_MM <- CalculateEfficiencyPCR(t,
					nt[w],
					annealingTemp,
					P,
					ions,
					batchSize,
					taqEfficiency=taqEfficiency,
					maxDistance,
					processors=processors)
				W <- which(is.na(primers$reverse_MM[w]))
				nnas <- which(!is.na(primers$reverse_MM[w]))
				if (length(nnas) > 0)
					W <- sort(c(W, nnas[which(eff_MM[nnas] > primers$reverse_efficiency_MM[w[nnas]])]))
				W <- W[which(eff_MM[W] >= 0.1)] # must have 10% MM efficiency
				if (length(W) > 0) {
					if (length(empty)==0) {
						primers$reverse_efficiency_MM[w[W]] <- eff_MM[W]
						primers$reverse_efficiency[w[W]] <- eff_PM[W]
						primers$reverse_MM[w[W]] <- paste(i,
							"/",
							strsplit(toString(reverseComplement(DNAStringSet(x[W]))),
								", ",
								fixed=TRUE)[[1]],
							sep="")
						primers$reverse_primer[w[W]] <- t[W]
					} else {
						primers$reverse_efficiency_MM[-empty][w[W]] <- eff_MM[W]
						primers$reverse_efficiency[-empty][w[W]] <- eff_PM[W]
						primers$reverse_MM[-empty][w[W]] <- paste(i,
							"/",
							strsplit(toString(reverseComplement(DNAStringSet(x[W]))),
								", ",
								fixed=TRUE)[[1]],
							sep="")
						primers$reverse_primer[-empty][w[W]] <- t[W]
					}
				}
			}
		}
		
		for (k in 1:length(ids)) {
			if (ids[k] == id)
				next
			
			w <- which(tiles$id==ids[k])
			target <- tiles[w,]
			
			if (any(uw != unique(target$width)))
				next
			
			combos_F <- numeric(d)
			count_F <- 1
			combos_R <- numeric(d)
			count_R <- 1
			nontargets_F <- targets_F <- character(d*sum(primers$permutations_forward)*max(table(target$start)))
			nontargets_R <- targets_R <- character(d*sum(primers$permutations_reverse)*max(table(target$start)))
			
			begin <- 1L
			for (i in 1:d) {
				# find all non-target_sites
				w <- .Call("multiMatchUpper", target$start_aligned, primers$start_aligned_forward[i], begin, PACKAGE="DECIPHER")
				
				if (length(w)==0) # target does not overlap nontarget
					break
				
				begin <- w[1]
				target_site <- target$target_site[w]
				
				j <- length(target_site)
				if (j==0) # no nontarget
					next
				
				w <- .Call("multiMatchCharNotNA", primers$forward_primer[i,], PACKAGE="DECIPHER")
				l <- length(w)
				
				if (l > 0) {
					combos_F[i] <- j*l
					targets_F[count_F:(count_F + combos_F[i] - 1)] <- rep(primers$forward_primer[i,w],
						each=j)
					nontargets_F[count_F:(count_F + combos_F[i] - 1)] <- rep(target_site,
						l)
					count_F <- count_F + combos_F[i]
				}
			}
			
			begin <- 1L#dim(target)[1]
			for (i in 1:d) {
				# find all non-target_sites
				w <- .Call("multiMatchUpper", target$end_aligned, primers$start_aligned_reverse[i], begin, PACKAGE="DECIPHER")
				
				if (length(w)==0) # target does not overlap nontarget
					break
				
				begin <- w[1]
				target_site <- target$target_site[w]
				
				j <- length(target_site)
				if (j==0) # no nontarget
					next
				
				w <- .Call("multiMatchCharNotNA", primers$reverse_primer[i,], PACKAGE="DECIPHER")
				l <- length(w)
				
				if (l > 0) {
					combos_R[i] <- j*l
					targets_R[count_R:(count_R + combos_R[i] - 1)] <- rep(primers$reverse_primer[i,w],
						each=j)
					nontargets_R[count_R:(count_R + combos_R[i] - 1)] <- rep(target_site,
						l)
					count_R <- count_R + combos_R[i]
				}
			}
			
			if ((count_F==1) && (count_R==1)) # sites do not overlap
				next
			
			if (count_F != 1) {
				eff_F <- CalculateEfficiencyPCR(targets_F[1:(count_F - 1)],
					strsplit(toString(reverseComplement(DNAStringSet(nontargets_F[1:(count_F - 1)]))),
						", ",
						fixed=TRUE)[[1]],
					annealingTemp,
					P,
					ions,
					batchSize,
					taqEfficiency=taqEfficiency,
					maxDistance,
					processors=processors)
				
				eff_max_F <- numeric(d)
				count_F <- 1
				for (i in 1:d) {
					if (combos_F[i] != 0) {
						eff_max_F[i] <- max(eff_F[count_F:(count_F + combos_F[i] - 1)])
						count_F <- count_F + combos_F[i]
					}
				}
				
				primers$score_forward <- primers$score_forward - eff_max_F
				
				w <- which(eff_max_F > .0001)
				if (length(w) > 0) {
					primers$mismatches_forward[w] <- paste(primers$mismatches_forward[w],
						ids[k],
						" (",
						signif(100*eff_max_F[w],
							digits=3),
						"%) ",
						sep="")
				}
			}
			
			if (count_R != 1) {
				eff_R <- CalculateEfficiencyPCR(targets_R[1:(count_R - 1)],
					nontargets_R[1:(count_R - 1)],
					annealingTemp,
					P,
					ions,
					batchSize,
					taqEfficiency=taqEfficiency,
					maxDistance,
					processors=processors)
				
				eff_max_R <- numeric(d)
				count_R <- 1
				for (i in 1:d) {
					if (combos_R[i] != 0) {
						eff_max_R[i] <- max(eff_R[count_R:(count_R + combos_R[i] - 1)])
						count_R <- count_R + combos_R[i]
					}
				}
				
				primers$score_reverse <- primers$score_reverse - eff_max_R
				
				w <- which(eff_max_R > .0001)
				if (length(w) > 0) {
					primers$mismatches_reverse[w] <- paste(primers$mismatches_reverse[w],
						ids[k],
						" (",
						signif(100*eff_max_R[w],
							digits=3),
						"%) ",
						sep="")
				}
			}
			
			if (worstScore > -Inf) {
				w_F <- primers$score_forward < worstScore
				if (any(w_F)) {
					primers$forward_primer[w_F,] <- NA
					primers$forward_efficiency[w_F,] <- NA
					primers$mismatches_forward[w_F] <- ""
					primers$score_forward[w_F] <- -Inf
					primers$permutations_forward[w_F] <- 0
					primers$forward_coverage[w_F,] <- NA
				}
				w_R <- primers$score_reverse < worstScore
				if (any(w_R)) {
					primers$reverse_primer[w_R,] <- NA
					primers$reverse_efficiency[w_R,] <- NA
					primers$mismatches_reverse[w_R] <- ""
					primers$score_reverse[w_R] <- -Inf
					primers$permutations_reverse[w_R] <- 0
					primers$reverse_coverage[w_R,] <- NA
				}
				w <- which(w_F & w_R)
				if (length(w) > 0)
					primers <- primers[-w,]
				d <- dim(primers)[1]
				if (d==0) {
					warning("All primers were below the worstScore: ",id)
					next
				}
			}
			
			if (verbose)
				setTxtProgressBar(pBar,
					floor(100*k/l_ids))
		}
		
		if (verbose)
			setTxtProgressBar(pBar, 100)
		
		if (numPrimerSets > 0) {
			if (verbose) {
				close(pBar)
				if (numPrimerSets > 1) {
					cat("Determining Best Primer Pairs:\n")
				} else {
					cat("Determining Best Primer Pair:\n")
				}
				flush.console()
				pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
			}
			
			# order primers by score
			searchSpace <- ifelse(numPrimerSets > 100, numPrimerSets, 100)
			f_F <- which(is.finite(primers$score_forward))
			o_F <- order(primers$score_forward[f_F],
				-1*primers$permutations_forward[f_F],
				rowSums(as.matrix(primers$forward_coverage[f_F,])),
				decreasing=TRUE)
			s_F <- ifelse(length(f_F) > searchSpace, searchSpace, length(f_F))
			f_R <- which(is.finite(primers$score_reverse))
			o_R <- order(primers$score_reverse[f_R],
				-1*primers$permutations_reverse[f_R],
				rowSums(as.matrix(primers$reverse_coverage[f_R,]), na.rm=TRUE),
				decreasing=TRUE)
			s_R <- ifelse(length(f_R) > searchSpace, searchSpace, length(f_R))
			
			# calculate the set's efficiency
			MMs_F <- strsplit(primers$mismatches_forward[f_F[o_F][1:s_F]], " (", fixed=TRUE)
			ls_F <- unlist(lapply(MMs_F, length))
			ls_F <- ifelse(ls_F > 0, ls_F - 1, 0)
			index_F <- rep(1:length(ls_F), ls_F)
			if (length(index_F) > 0) {
				MMs_F <- unlist(strsplit(unlist(MMs_F), "%) ", fixed=TRUE))
				effs_F <- as.numeric(MMs_F[seq(2, length(MMs_F), 2)])
				MMs_F <- MMs_F[seq(1, length(MMs_F), 2)]
			}
			
			MMs_R <- strsplit(primers$mismatches_reverse[f_R[o_R][1:s_R]], " (", fixed=TRUE)
			ls_R <- unlist(lapply(MMs_R, length))
			ls_R <- ifelse(ls_R > 0, ls_R - 1, 0)
			index_R <- rep(1:length(ls_R), ls_R)
			if (length(index_R) > 0) {
				MMs_R <- unlist(strsplit(unlist(MMs_R), "%) ", fixed=TRUE))
				effs_R <- as.numeric(MMs_R[seq(2, length(MMs_R), 2)])
				MMs_R <- MMs_R[seq(1, length(MMs_R), 2)]
			}
			
			m <- matrix(0,
				nrow=s_F,
				ncol=s_R,
				dimnames=list(f_F[o_F][1:s_F],
					f_R[o_R][1:s_R]))
			for (i in 1:s_F) {
				w_F <- which(index_F==i)
				if (length(w_F)==0)
					next
				for (j in 1:s_R) {
					w_R <- which(index_R==j)
					if (length(w_R)==0)
						next
					MMs <- match(MMs_F[w_F], MMs_R[w_R])
					w <- which(!is.na(MMs))
					if (length(w) > 0)
						m[i, j] <- sum(tapply(c(effs_F[w_F[w]], effs_R[w_R[MMs[w]]]),
							rep(1:length(w), 2),
							function(x) {return(prod(x)/100)}))
				}
			}
			
			# determine which primers will meet the specified constraints
			dims <- dim(m)
			p <- outer(primers$permutations_forward[f_F[o_F][1:s_F]],
				primers$permutations_reverse[f_R[o_R][1:s_R]],
				FUN="+")
			c <- -1*outer(ifelse(rep(s_F, s_F)==1,
					sum(primers$forward_coverage[f_F[o_F][1:s_F],], na.rm=TRUE),
					rowSums(as.matrix(primers$forward_coverage[f_F[o_F][1:s_F],]), na.rm=TRUE)),
				ifelse(rep(s_R, s_R)==1,
					sum(primers$reverse_coverage[f_R[o_R][1:s_R],], na.rm=TRUE),
					rowSums(as.matrix(primers$reverse_coverage[f_R[o_R][1:s_R],]), na.rm=TRUE)),
				FUN="*")
			s <- -1*outer(primers$score_forward[f_F[o_F][1:s_F]],
				primers$score_reverse[f_R[o_R][1:s_R]],
				FUN="+")
			z <- -1*outer(primers$start_forward[f_F[o_F][1:s_F]],
				primers$start_reverse[f_R[o_R][1:s_R]],
				FUN="-") - 1
			z <- ifelse(z > minProductSize & z < maxProductSize, 0, 1)
			o <- order(z, m, p, c, s)
			g_F <- g_R <- integer()
			space <- 0
			starts_F <- starts_R <- integer()
			for (i in 1:length(o)) {
				w_o <- c((o[i] - 1)%%dims[1] + 1,
					(o[i] - 1)%/%dims[1] + 1)
				
				# check for primer-dimer formation
				seqs1 <- c(rep(primers$forward_primer[f_F[o_F][w_o[1]],],
						maxPermutations),
					rep(primers$forward_primer[f_F[o_F][w_o[1]],],
						maxPermutations),
					rep(primers$reverse_primer[f_R[o_R][w_o[2]],],
						maxPermutations))
				seqs2 <- c(rep(primers$reverse_primer[f_R[o_R][w_o[2]],],
						each=maxPermutations),
					rep(primers$forward_primer[f_F[o_F][w_o[1]],],
						each=maxPermutations),
					rep(primers$reverse_primer[f_R[o_R][w_o[2]],],
						each=maxPermutations))
				remove <- which(is.na(seqs1) | is.na(seqs2))
				if (length(remove) > 0) {
					seqs1 <- seqs1[-remove]
					seqs2 <- seqs2[-remove]
				}
				eff_pd <- .primerDimer(seqs1, seqs2, annealingTemp, P, ions, processors=processors)
				if (any(eff_pd >= primerDimer)) # primer-dimer
					next
				
				start_F <- primers$start_forward[f_F[o_F][w_o[1]]]
				start_R <- primers$start_reverse[f_R[o_R][w_o[2]]]
				productSize <- start_R - start_F + 1
				if ((productSize > minProductSize) &&
					(productSize < maxProductSize)) {
					space <- space + 1
					if (length(which(g_F==w_o[1])) == 0) {
						g_F <- c(g_F, w_o[1])
						starts_F <- min(starts_F,
							primers$start_aligned_forward[f_F[o_F][w_o[1]]])
					}
					if (length(which(g_R==w_o[2])) == 0) {
						g_R <- c(g_R, w_o[2])
						starts_R <- max(starts_R,
							primers$start_aligned_reverse[f_R[o_R][w_o[2]]])
					}
				}
				if (space >= numPrimerSets)
					break
			}
			s_F <- length(g_F)
			s_R <- length(g_R)
			
			if ((s_F == 0) || (s_R == 0)) {
				warning("Not enough primers to meet the specified constraints: ", id)
				next
			}
			
			# exhaustively rescore top forward primers
			for (i in 1:s_F) {
				# initialize mismatch scores
				primers$mismatches_forward[f_F[o_F][g_F[i]]] <- ""
				primers$score_forward[f_F[o_F][g_F[i]]] <- 0
				for (k in 1:length(ids)) {
					if (ids[k] == id)
						next
					
					w <- which(tiles$id==ids[k])
					if (length(w) > 0)
						w <- w[which(tiles$end_aligned[w] < starts_R)]
					#if (length(w) > 0)
					#	w <- w[which(tiles$end[w] > (max(tiles$end[w]) - maxSearchSize))]
					if (length(w) > 0)
						w <- w[which((tiles$start[w] >= primers$start_forward[f_F[o_F][g_F[i]]] - maxSearchSize) &
							(tiles$start[w] <= primers$start_forward[f_F[o_F][g_F[i]]] + maxSearchSize))]
					if (length(w)==0)
						next
					target <- tiles[w,]
					
					if (any(uw != unique(target$width)))
						next
					
					w <- .Call("multiMatchCharNotNA", primers$forward_primer[f_F[o_F][g_F[i]],], PACKAGE="DECIPHER")
					l <- length(target$target_site)
					targets <- rep(primers$forward_primer[f_F[o_F][g_F[i]],w], each=l)
					nontargets <- strsplit(toString(reverseComplement(DNAStringSet(rep(target$target_site, length(w))))),
						", ",
						fixed=TRUE)[[1]]
					
					eff <- CalculateEfficiencyPCR(targets,
						nontargets,
						annealingTemp,
						P,
						ions,
						batchSize,
						taqEfficiency=taqEfficiency,
						maxDistance,
						processors=processors)
					
					if (ragged5Prime) {
						w_max <- which(eff==max(eff))
						w_max <- w_max[length(w_max)]
					} else {
						w_max <- which.max(eff)
					}
					eff_max <- eff[w_max]
					
					primers$score_forward[f_F[o_F][g_F[i]]] <- primers$score_forward[f_F[o_F][g_F[i]]] - eff_max
					
					if (eff_max > .0001) {
						s1 <- targets[w_max]
						s2 <- reverseComplement(DNAStringSet(nontargets[w_max]))
						p <- pairwiseAlignment(s1,
							paste("----", s2, "----", sep=""),
							type="global-local",
							gapOpen=-5,
							gapExtension=-5)
						s1 <- toString(pattern(p))
						s2 <- toString(subject(p))
						s2 <- toString(reverseComplement(DNAString(s2)))
						primers$mismatches_forward[f_F[o_F][g_F[i]]] <- paste(primers$mismatches_forward[f_F[o_F][g_F[i]]],
							ids[k],
							" (",
							signif(100*eff_max,
								digits=3),
							"%,",
							s1,
							"/",
							s2,
							") ",
							sep="")
					}
					
					if (verbose)
						setTxtProgressBar(pBar, floor(100*((i - 1)*l_ids + k)/(s_F + s_R)/l_ids))
				}
			}
			
			# exhaustively rescore top reverse primers
			for (i in 1:s_R) {
				# initialize mismatch scores
				primers$mismatches_reverse[f_R[o_R][g_R[i]]] <- ""
				primers$score_reverse[f_R[o_R][g_R[i]]] <- 0
				for (k in 1:length(ids)) {
					if (ids[k] == id)
						next
					
					w <- which(tiles$id==ids[k])
					if (length(w) > 0)
						w <- w[which(tiles$start_aligned[w] > starts_F)]
					#if (length(w) > 0)
					#	w <- w[which(tiles$start[w] < (min(tiles$start[w]) + maxSearchSize))]
					if (length(w) > 0)
						w <- w[which((tiles$end[w] >= primers$start_reverse[f_R[o_R][g_R[i]]] - maxSearchSize) &
							(tiles$end[w] <= primers$start_reverse[f_R[o_R][g_R[i]]] + maxSearchSize))]
					if (length(w)==0)
						next
					target <- tiles[w,]
					
					if (any(uw != unique(target$width)))
						next
					
					w <- .Call("multiMatchCharNotNA", primers$reverse_primer[f_R[o_R][g_R[i]],], PACKAGE="DECIPHER")
					l <- length(target$target_site)
					targets <- rep(primers$reverse_primer[f_R[o_R][g_R[i]],w], each=l)
					nontargets <- rep(target$target_site, length(w))
					
					eff <- CalculateEfficiencyPCR(targets,
						nontargets,
						annealingTemp,
						P,
						ions,
						batchSize,
						taqEfficiency=taqEfficiency,
						maxDistance,
						processors=processors)
					
					if (ragged5Prime) {
						w_max <- which.max(eff)
					} else {
						w_max <- which(eff==max(eff))
						w_max <- w_max[length(w_max)]
					}
					
					eff_max <- eff[w_max]
					
					primers$score_reverse[f_R[o_R][g_R[i]]] <- primers$score_reverse[f_R[o_R][g_R[i]]] - eff_max
					
					if (eff_max > .0001) {
						s1 <- targets[w_max]
						s2 <- reverseComplement(DNAStringSet(nontargets[w_max]))
						p <- pairwiseAlignment(s1,
							paste("----", s2, "----", sep=""),
							type="global-local",
							gapOpen=-5,
							gapExtension=-5)
						s1 <- toString(pattern(p))
						s2 <- toString(subject(p))
						s2 <- toString(reverseComplement(DNAString(s2)))
						primers$mismatches_reverse[f_R[o_R][g_R[i]]] <- paste(primers$mismatches_reverse[f_R[o_R][g_R[i]]],
							ids[k],
							" (",
							signif(100*eff_max,
								digits=3),
							"%,",
							s1,
							"/",
							s2,
							") ",
							sep="")
					}
					
					if (verbose)
						setTxtProgressBar(pBar, floor(100*((i + s_F - 1)*l_ids + k)/(s_F + s_R)/l_ids))
				}
			}
			
			# calculate the set's efficiency
			MMs_F <- strsplit(primers$mismatches_forward[f_F[o_F][g_F[1:s_F]]], " (", fixed=TRUE)
			ls_F <- unlist(lapply(MMs_F, length))
			ls_F <- ifelse(ls_F > 0, ls_F - 1, 0)
			index_F <- rep(1:length(ls_F), ls_F)
			if (length(index_F) > 0) {
				MMs_F <- strsplit(unlist(MMs_F), "%,", fixed=TRUE)
				MMs_F <- unlist(strsplit(unlist(MMs_F), ") ", fixed=TRUE))
				effs_F <- as.numeric(MMs_F[seq(2, length(MMs_F), 3)])
				MMs_F <- MMs_F[seq(1, length(MMs_F), 3)]
			}
			
			MMs_R <- strsplit(primers$mismatches_reverse[f_R[o_R][g_R[1:s_R]]], " (", fixed=TRUE)
			ls_R <- unlist(lapply(MMs_R, length))
			ls_R <- ifelse(ls_R > 0, ls_R - 1, 0)
			index_R <- rep(1:length(ls_R), ls_R)
			if (length(index_R) > 0) {
				MMs_R <- strsplit(unlist(MMs_R), "%,", fixed=TRUE)
				MMs_R <- unlist(strsplit(unlist(MMs_R), ") ", fixed=TRUE))
				effs_R <- as.numeric(MMs_R[seq(2, length(MMs_R), 3)])
				MMs_R <- MMs_R[seq(1, length(MMs_R), 3)]
			}
			
			m <- matrix(0,
				nrow=s_F,
				ncol=s_R,
				dimnames=list(f_F[o_F][g_F[1:s_F]],
					f_R[o_R][g_R[1:s_R]]))
			for (i in 1:s_F) {
				w_F <- which(index_F==i)
				if (length(w_F)==0)
					next
				for (j in 1:s_R) {
					w_R <- which(index_R==j)
					if (length(w_R)==0)
						next
					MMs <- match(MMs_F[w_F], MMs_R[w_R])
					w <- which(!is.na(MMs))
					if (length(w) > 0)
						m[i, j] <- sum(tapply(c(effs_F[w_F[w]], effs_R[w_R[MMs[w]]]),
							rep(1:length(w), 2),
							function(x) {return(prod(x)/100)}))
				}
			}
			
			# choose the best primer sets
			d <- dimnames(m)
			w <- which(pSets$identifier==id)
			p <- outer(primers$permutations_forward[f_F[o_F][g_F[1:s_F]]],
				primers$permutations_reverse[f_R[o_R][g_R[1:s_R]]],
				FUN="+")
			c <- -1*outer(ifelse(rep(s_F, s_F)==1,
					sum(primers$forward_coverage[f_F[o_F][g_F[1:s_F]],], na.rm=TRUE),
					rowSums(as.matrix(primers$forward_coverage[f_F[o_F][g_F[1:s_F]],]), na.rm=TRUE)),
				ifelse(rep(s_R, s_R)==1,
					sum(primers$reverse_coverage[f_R[o_R][g_R[1:s_R]],], na.rm=TRUE),
					rowSums(as.matrix(primers$reverse_coverage[f_R[o_R][g_R[1:s_R]],]), na.rm=TRUE)),
				FUN="*")
			s <- -1*outer(primers$score_forward[f_F[o_F][g_F[1:s_F]]],
				primers$score_reverse[f_R[o_R][g_R[1:s_R]]],
				FUN="+")
			o <- order(m, p, c, s)
			j <- 0
			dims <- dim(m)
			for (k in 1:numPrimerSets) {
				if ((j + 1) > length(o)) {
					warnings("Not enough primer sets meet the specified constaints:", id)
					break
				}
				for (j in (j + 1):length(o)) {
					w_o <- c((o[j] - 1)%%dims[1] + 1,
						(o[j] - 1)%/%dims[1] + 1)
					f <- as.numeric(d[[1]][w_o[1]])
					r <- as.numeric(d[[2]][w_o[2]])
					
					# check for primer-dimer formation
					seqs1 <- c(rep(primers$forward_primer[f,],
							maxPermutations),
						rep(primers$forward_primer[f,],
							maxPermutations),
						rep(primers$reverse_primer[r,],
							maxPermutations))
					seqs2 <- c(rep(primers$reverse_primer[r,],
							each=maxPermutations),
						rep(primers$forward_primer[f,],
							each=maxPermutations),
						rep(primers$reverse_primer[r,],
							each=maxPermutations))
					remove <- which(is.na(seqs1) | is.na(seqs2))
					if (length(remove) > 0) {
						seqs1 <- seqs1[-remove]
						seqs2 <- seqs2[-remove]
					}
					eff_pd <- .primerDimer(seqs1, seqs2, annealingTemp, P, ions, processors=processors)
					if (any(eff_pd >= primerDimer)) # primer-dimer
						next
					
					start_F <- primers$start_forward[f]
					start_R <- primers$start_reverse[r]
					productSize <- start_R - start_F + 1
					if ((productSize > minProductSize) &&
						(productSize < maxProductSize) &&
						!(any(eff_pd >= primerDimer)))
						break
				}
				
				if ((productSize > minProductSize) &&
					(productSize < maxProductSize) &&
					!(any(eff_pd >= primerDimer))) {
					pSets$start_forward[w[k]] <- primers$start_forward[f]
					pSets$start_reverse[w[k]] <- primers$start_reverse[r]
					pSets$product_size[w[k]] <- productSize
					pSets$start_aligned_forward[w[k]] <- primers$start_aligned_forward[f]
					pSets$start_aligned_reverse[w[k]] <- primers$start_aligned_reverse[r]
					pSets$permutations_forward[w[k]] <- primers$permutations_forward[f]
					pSets$permutations_reverse[w[k]] <- primers$permutations_reverse[r]
					pSets$score_forward[w[k]] <- primers$score_forward[f]
					pSets$score_reverse[w[k]] <- primers$score_reverse[r]
					pSets$score_set[w[k]] <- -m[o[j]]/100
					pSets$forward_primer[w[k],] <- primers$forward_primer[f,]
					pSets$reverse_primer[w[k],] <- primers$reverse_primer[r,]
					pSets$forward_efficiency[w[k],] <- primers$forward_efficiency[f,]
					pSets$reverse_efficiency[w[k],] <- primers$reverse_efficiency[r,]
					pSets$forward_coverage[w[k],] <- primers$forward_coverage[f,]
					pSets$reverse_coverage[w[k],] <- primers$reverse_coverage[r,]
					pSets$mismatches_forward[w[k]] <- primers$mismatches_forward[f]
					pSets$mismatches_reverse[w[k]] <- primers$mismatches_reverse[r]
					
					if (induceMismatch) {
						pSets$forward_MM[w[k],] <- primers$forward_MM[f,]
						pSets$forward_efficiency_MM[w[k],] <- primers$forward_efficiency_MM[f,]
						pSets$reverse_MM[w[k],] <- primers$reverse_MM[r,]
						pSets$reverse_efficiency_MM[w[k],] <- primers$reverse_efficiency_MM[r,]
					}
					
					w_F <- which(index_F==w_o[1])
					w_R <- which(index_R==w_o[2])
					MMs <- match(MMs_F[w_F], MMs_R[w_R])
					w_S <- which(!is.na(MMs))
					if (length(w_S)==0)
						next
					effs <- tapply(c(effs_F[w_F[w_S]], effs_R[w_R[MMs[w_S]]]),
							rep(1:length(w_S), 2),
							function(x) {return(prod(x)/100)})
					# record set mismatches
					w1 <- which(effs >= .1)
					if (length(w1) > 0)
						pSets$mismatches_set[w[k]] <- paste(MMs_F[w_F][w_S][w1],
							" (",
							signif(effs[w1],
								digits=1),
							"%)",
							sep="",
							collapse=", ")
				} else {
					warning("Not enough primers to meet the specified constraints: ", id)
					break
				}
			}
			
			if (verbose)
				setTxtProgressBar(pBar, 100)
		}
		
		if (verbose) {
			close(pBar)
			time.2 <- Sys.time()
			cat("\n")
			print(signif(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
		}
		
		if (numPrimerSets == 0)
			primers_all <- rbind(primers_all, primers)
	}
	
	if (numPrimerSets > 0) {
		w <- which(pSets$start_forward==0)
		if (length(w) > 0)
			pSets <- pSets[-w,]
		return(pSets)
	} else {
		w <- which(primers_all$start_forward==0)
		if (length(w) > 0)
			primers_all <- primers_all[-w,]
		return(primers_all)
	}
}
