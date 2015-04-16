CalculateEfficiencyPCR <- function(primer,
	target,
	temp,
	P,
	ions,
	batchSize=1000,
	taqEfficiency=TRUE,
	maxDistance=0.4,
	maxGaps=2,
	processors=NULL) {
	
	# error checking
	if (is.character(primer))
		primer <- toupper(primer)
	if (is(primer, "DNAStringSet"))
		primer <- strsplit(toString(primer), ", ", fixed=TRUE)[[1]]
	if (is(target, "DNAStringSet"))
		target <- strsplit(toString(target), ", ", fixed=TRUE)[[1]]
	if (!is.character(primer))
		stop("primer must be a character vector or DNAStringSet.")
	if (!is.character(target))
		stop("target must be a character vector or DNAStringSet.")
	if (!is.numeric(ions))
		stop("ions must be a numeric.")
	if (ions < .01 || is.nan(ions))
		stop("Sodium equivilent concentration must be at least 0.01M.")
	if (!is.numeric(P))
		stop("P must be a numeric.")
	if (!(P > 0))
		stop("P must be greater than zero.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize)!=batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
	if (!is.numeric(temp))
		stop("temp must be a numeric.")
	if (!is.logical(taqEfficiency))
		stop("taqEfficiency must be a logical.")
	
	RT <- .0019871*(273.15 + temp) # [kcal/mol]
	l <- length(primer)
	
	if (l==0)
		stop("No primer specified.")
	if (l!=length(target))
		stop("primer is not the same length as target.")
	
	# align primer and target
	seqs2 <- reverseComplement(DNAStringSet(target))
	seqs2 <- unlist(strsplit(toString(seqs2), ", ", fixed=TRUE))
	seqs2 <- paste("----", seqs2, "----", sep="")
	p <- pairwiseAlignment(primer,
		seqs2,
		type="global-local",
		gapOpen=-5,
		gapExtension=-5)
	t_START <- nchar(seqs2) - p@subject@range@start - p@subject@range@width - 2
	t_END <- nchar(seqs2) - p@subject@range@start - 3
	
	if (taqEfficiency) {
		if (!is.numeric(maxDistance))
			stop("maxDistance must be a numeric.")
		if (maxDistance < 0)
			stop("maxDistance must be greater or equal to zero.")
		if (maxDistance > 1)
			stop("maxDistance must be less than or equal to one.")
		if (!is.numeric(maxGaps))
			stop("maxGaps must be a numeric.")
		if (maxGaps < 0)
			stop("maxGaps must be at least zero.")
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
		
		seqs1 <- unlist(strsplit(toString(pattern(p)), ", ", fixed=TRUE))
		seqs2 <- unlist(strsplit(toString(subject(p)), ", ", fixed=TRUE))
		seqs2 <- reverseComplement(DNAStringSet(seqs2))
		seqs2 <- unlist(strsplit(toString(seqs2), ", ", fixed=TRUE))
		
		# calculate elongation efficiency
		eff_taq <- .Call("terminalMismatch", seqs1, seqs2, maxDistance, maxGaps, processors, PACKAGE="DECIPHER")
	} else {
		eff_taq <- numeric(l)
		eff_taq[] <- 1
	}
	t <- which(eff_taq >= .001)
	tl <- length(t)
	if (tl == 0)
		return(eff_taq)
	
	# duplex formation
	seqs <- paste(primer[t],
		target[t],
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
	dG1 <- numeric(l)
	dG1[t] <- dG[match(seqs, seq)]
	K1 <- exp(-dG1/RT)
	
	eff <- P*K1/(1 + P*K1)*eff_taq
	if (tl != l)
		eff[-t] <- eff_taq[-t]
	z <- which(eff >= .001)
	l <- length(z)
	if (l == 0)
		return(eff)
	primer <- primer[z]
	target <- target[z]
	K1 <- K1[z]
	dG1 <- dG1[z]
	eff_taq <- eff_taq[z]
	t_START <- t_START[z]
	t_END <- t_END[z]
	
	# primer folding
	seqs <- primer
	seq <- unique(seqs)
	ls <- length(seq)
	dG <- numeric(ls)
	for (start in seq(1, ls, batchSize)) {
		end <- ifelse(start + batchSize > ls, ls, start + batchSize)
		dG[start:end] <- as.numeric(system(paste("hybrid-ss-min -n DNA -t",
				temp,
				"-T",
				temp,
				"-N",
				ions,
				"-E -q",
				paste(seq[start:end], collapse=" ")),
			intern=TRUE))
	}
	dG2a <- numeric(l)
	dG2a <- dG[match(seqs, seq)]
	K2a <- exp(-dG2a/RT)
	
	# primer-primer duplex
	seqs <- paste(primer,
		primer,
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
	dG2b <- numeric(l)
	dG2b <- dG[match(seqs, seq)]
	K2b <- exp(-dG2b/RT)
	
	# target folding
	seqs <- substr(target, t_START, t_END)
	seq <- unique(seqs)
	ls <- length(seq)
	dG <- numeric(ls)
	for (start in seq(1, ls, batchSize)) {
		end <- ifelse(start + batchSize > ls, ls, start + batchSize)
		dG[start:end] <- as.numeric(system(paste("hybrid-ss-min -n DNA -t",
				temp,
				"-T",
				temp,
				"-N",
				ions,
				"-E -q",
				paste(seq[start:end], collapse=" ")),
			intern=TRUE))
	}
	dG3 <- numeric(l)
	dG3 <- dG[match(seqs, seq)]
	K3 <- exp(-dG3/RT)
	
	K2eff <- numeric(l)
	K2eff[] <- -1
	w <- which(K2b==0)
	K2eff[w] <- 0
	w <- which(K2b*P/K2a < .01)
	K2eff[w] <- K2a[w]
	w <- which(8*K2b*P > (1 + K2a)^2)
	K2eff[w] <- K2b[w]
	w <- which(K2eff==-1)
	K2eff[w] <- 4*K2b[w]*P/(-1 - K2a[w] + sqrt((1 + K2a[w])^2 - 8*K2b[w]*P)) - 1
	w <- which(K2eff < 0)
	K2eff[w] <- 0
	
	Kov <- K1/((1 + K2eff)*(1 + K3))
	
	eff[z] <- P*Kov/(1 + P*Kov)*eff_taq
	return(eff)
}
