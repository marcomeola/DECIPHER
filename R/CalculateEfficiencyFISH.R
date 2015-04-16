CalculateEfficiencyFISH <- function(probe,
	target,
	temp,
	P,
	ions,
	FA,
	batchSize=1000) {
	
	# error checking
	if (is.character(probe))
		probe <- toupper(probe)
	if (is(probe, "DNAStringSet"))
		probe <- strsplit(toString(probe), ", ", fixed=TRUE)[[1]]
	if (is.character(target))
		target <- .Call("replaceChar",
			target,
			"U",
			"T",
			PACKAGE="DECIPHER")
	if (is(target, "RNAStringSet"))
		target <- DNAStringSet(target)
	if (is(target, "DNAStringSet"))
		target <- strsplit(toString(target), ", ", fixed=TRUE)[[1]]
	if (!is.character(probe))
		stop("probe must be a DNAStringSet or character vector.")
	if (!is.character(target))
		stop("target must be a DNAStringSet or character vector.")
	if (!is.numeric(ions))
		step("ions must be a numeric.")
	if (ions < .01 || is.nan(ions))
		stop("Sodium equivilent concentration must be at least 0.01M.")
	if (!is.numeric(P))
		stop("P must be a numeric.")
	if (!(P > 0))
		stop("P must be greater than zero.")
	if (!is.numeric(FA))
		stop("FA must be a numeric.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize)!=batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
	if (!is.numeric(temp))
		stop("temp must be a numeric.")
	
	RT <- .0019871*(273.15 + temp) # [kcal/mol]
	l <- length(probe)
	n <- nchar(probe)
	
	if (l==0)
		stop("No probe specified.")
	if (l!=length(target))
		stop("probe is not the same length as target.")
	
	# align probe and target
	seqs2 <- reverseComplement(DNAStringSet(target))
	seqs2 <- unlist(strsplit(toString(seqs2), ", ", fixed=TRUE))
	seqs2 <- paste("----", seqs2, "----", sep="")
	p <- pairwiseAlignment(probe,
		seqs2,
		type="global-local",
		gapOpen=-10,
		gapExtension=-10)
	seqs1 <- unlist(strsplit(toString(pattern(p)), ", ", fixed=TRUE))
	seqs2 <- unlist(strsplit(toString(subject(p)), ", ", fixed=TRUE))
	
	deltas <- .Call("calculateFISH", seqs1, seqs2, PACKAGE="DECIPHER")
	dG1_PM_DNARNA <- deltas[,1] - (273.15 + temp)/1000*(deltas[,2] + 0.368*n*log(ions))
	
	# determine ddG1 for mismatched probes
	ddG1_MM_DNARNA <- numeric(l)
	ddG1_MM_DNADNA <- numeric(l)
	ddG1_MM_RNARNA <- numeric(l)
	dG1_PM_UNADNA <- numeric(l)
	dG1_MM_UNADNA <- numeric(l)
	dG1_PM_UNARNA <- numeric(l)
	dG1_MM_UNARNA <- numeric(l)
	MM <- which(deltas[,3] != 0) # mismatched probes/targets
	if (length(MM) > 0) {
		ddG1_MM_DNARNA[MM] <- deltas[MM,3] - (273.15 + temp)/1000*(deltas[MM,4] + 0.368*n[MM]*log(ions))
		ddG1_MM_DNADNA[MM] <- deltas[MM,5] - (273.15 + temp)/1000*(deltas[MM,6] + 0.368*n[MM]*log(ions))
		ddG1_MM_RNARNA[MM] <- deltas[MM,7] - (273.15 + temp)/1000*(deltas[MM,8] + 0.368*n[MM]*log(ions))
		
		target <- .Call("replaceChar", seqs2[MM], "-", "", PACKAGE="DECIPHER")
		target <- reverseComplement(DNAStringSet(target))
		target <- unlist(strsplit(toString(target), ", ", fixed=TRUE))
		target_PM <- unlist(strsplit(toString(reverseComplement(DNAStringSet(probe[MM]))), ", ", fixed=TRUE))
		
		seqs <- paste(probe[MM],
			target_PM,
			sep=" ")
		seq <- unique(seqs)
		ls <- length(seq)
		dG <- numeric(ls)
		for (start in seq(1, ls, batchSize)) {
			end <- ifelse(start + batchSize > ls, ls, start + batchSize)
			dG[start:end] <- as.numeric(system(paste("hybrid-min -n DNA -t",
					temp,
					"-T",
					temp,
					"-N",
					ions,
					"-E -q",
					paste(seq[start:end], collapse=" ")),
				intern=TRUE))
		}
		dG1_PM_UNADNA[MM] <- dG[match(seqs, seq)]
		
		seqs <- paste(probe[MM],
			target,
			sep=" ")
		seq <- unique(seqs)
		ls <- length(seq)
		dG <- numeric(ls)
		for (start in seq(1, ls, batchSize)) {
			end <- ifelse(start + batchSize > ls, ls, start + batchSize)
			dG[start:end] <- as.numeric(system(paste("hybrid-min -n DNA -t",
					temp,
					"-T",
					temp,
					"-N",
					ions,
					"-E -q",
					paste(seq[start:end], collapse=" ")),
				intern=TRUE))
		}
		dG1_MM_UNADNA[MM] <- dG[match(seqs, seq)]
		
		seqs <- paste(probe[MM],
			target_PM,
			sep=" ")
		seq <- unique(seqs)
		ls <- length(seq)
		dG <- numeric(ls)
		for (start in seq(1, ls, batchSize)) {
			end <- ifelse(start + batchSize > ls, ls, start + batchSize)
			dG[start:end] <- as.numeric(system(paste("hybrid-min -n RNA -t",
					temp,
					"-T",
					temp,
					"-E -q",
					paste(seq[start:end], collapse=" ")),
				intern=TRUE))
		}
		dG1_PM_UNARNA[MM] <- dG[match(seqs, seq)]
		
		seqs <- paste(probe[MM],
			target,
			sep=" ")
		seq <- unique(seqs)
		ls <- length(seq)
		dG <- numeric(ls)
		for (start in seq(1, ls, batchSize)) {
			end <- ifelse(start + batchSize > ls, ls, start + batchSize)
			dG[start:end] <- as.numeric(system(paste("hybrid-min -n RNA -t",
					temp,
					"-T",
					temp,
					"-E -q",
					paste(seq[start:end], collapse=" ")),
				intern=TRUE))
		}
		dG1_MM_UNARNA[MM] <- dG[match(seqs, seq)]
	}
	
	ddG1_DNA <- dG1_MM_UNADNA - dG1_PM_UNADNA
	ddG1_RNA <- dG1_MM_UNARNA - dG1_PM_UNARNA
	ddG1_loop_DNA <- ddG1_DNA + ddG1_MM_DNADNA
	ddG1_loop_RNA <- ddG1_RNA + ddG1_MM_RNARNA
	ddG1 <- (ddG1_loop_DNA + ddG1_loop_RNA)/2 - ddG1_MM_DNARNA
	dG1 <- dG1_PM_DNARNA + ddG1
	dG1 <- 0.2558*dG1 - 6.4867
	K1 <- exp(-(dG1 + FA*(0.0175 + 0.0028*n))/RT)
	eff <- P*K1/(1 + P*K1)
	
	FAm <- numeric(l)
	FAm[] <- NA
	for (i in 1:l) {
		f <- function(FA) {
			K1 <- exp(-(dG1[i] + FA*(0.0175 + 0.0028*n[i]))/RT)
			return(0.5 - P*K1/(1 + P*K1))
		}
		try(FAm[i] <- uniroot(f,
				c(-1000,
					1000))$root,
			silent=TRUE)
	}
	
	ans <- matrix(nrow=l,
		ncol=4,
		dimnames=list(1:l,
			c("HybEff",
				"FAm",
				"ddG1",
				"dG1")))
	
	ans[, "HybEff"] <- eff
	ans[, "FAm"] <- FAm
	ans[, "ddG1"] <- ddG1
	ans[, "dG1"] <- dG1
	
	return(ans)
}
