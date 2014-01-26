AlignTranslation <- function(myXStringSet,
	sense="+",
	direction="5' to 3'",
	readingFrame=NA,
	asAAStringSet=FALSE,
	...) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	if (min(width(myXStringSet)) < 2)
		stop("All sequences in myXStringSet must be at least two nucleotides long.")
	if (length(myXStringSet) < 2)
		stop("At least two sequences are required.")
	if (sense != "+" && sense != "-")
		stop('sense must be either "+" or "-".')
	if (direction != "5' to 3'" && direction != "3' to 5'")
		stop('direction must be either "5\' to 3\'" or "3\' to 5\'".')
	if (!is.na(readingFrame) && !is.numeric(readingFrame))
		stop('readingFrame must be a numeric.')
	if (length(readingFrame) != 1 && length(readingFrame) != length(myXStringSet))
		stop('readingFrame must be a single numeric or the length of myXStringSet.')
	if (any(!is.na(readingFrame) & (readingFrame > 3 | readingFrame < 1)))
		stop('readingFrame must be either NA, 1, 2, or 3.')
	if (any(!is.na(readingFrame) & floor(readingFrame) != readingFrame))
		stop('readingFrame must be either NA, 1, 2, or 3.')
	if (!is.logical(asAAStringSet))
		stop('asAAStringSet must be a logical.')
	a <- vcountPattern("-", myXStringSet)
	if (any(a) > 0)
		stop("Gaps must be removed before alignment.")
	
	lkupTable <- c(rep("A", 4), rep("R", 6), rep("N", 2), rep("D", 2), rep("C", 2), rep("Q", 2), rep("E", 2), rep("G", 4), rep("H", 2), rep("I", 3), rep("L", 6), rep("K", 2), "M", rep("F", 2), rep("P", 4), rep("S", 6), rep("T", 4), "W", rep("Y", 2), rep("V", 4), rep("*", 3))
	if (is(myXStringSet, "DNAStringSet")) {
		names(lkupTable) <- c("GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "AAT", "AAC", "GAT", "GAC", "TGT", "TGC", "CAA", "CAG", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG", "CAT", "CAC", "ATT", "ATC", "ATA", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "AAA", "AAG", "ATG", "TTT", "TTC", "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "ACT", "ACC", "ACA", "ACG", "TGG", "TAT", "TAC", "GTT", "GTC", "GTA", "GTG", "TAA", "TGA", "TAG")
	} else { # RNAStringSet
		names(lkupTable) <- c("GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC", "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA", "GUG", "UAA", "UGA", "UAG")
	}
	
	if (sense=="-")
		myXStringSet <- reverseComplement(myXStringSet)
	if (direction=="3' to 5'")
		myXStringSet <- reverse(myXStringSet)
	if (length(readingFrame)==1)
		readingFrame <- rep(readingFrame, length(myXStringSet))
	
	AA <- list()
	rFs <- numeric(length(myXStringSet))
	for (i in 1:length(myXStringSet)) {
		if (is.na(readingFrame[i])) {
			index <- c(0:2, 0) # need to return to zero by default
		} else {
			index <- readingFrame[i] - 1
		}
		starts <- seq(1, width(myXStringSet)[i], 3)
		for (rF in index) {
			rFs[i] <- rF
			codons <- substring(myXStringSet[i], starts + rF, starts + rF + 2)
			if (nchar(codons[length(codons)]) < 3)
				length(codons) <- length(codons) - 1L
			Ns <- grep("[^A|C|G|T]", codons)
			if (length(Ns) > 0) {
				AA[[i]] <- character(length(codons))
				AA[[i]][Ns] <- "X"
				AA[[i]][-Ns] <- lkupTable[codons[-Ns]]
			} else {
				AA[[i]] <- lkupTable[codons]
			}
			if (length(which(AA[[i]]=="*")) < 2) {
				break
			}
		}
	}
	
	AA <- AAStringSet(unlist(lapply(AA, paste, collapse="")))
	AA <- AlignSeqs(myXStringSet=AA, ...)
	if (asAAStringSet) {
		names(AA) <- names(myXStringSet)
		return(AA)
	}
	
	gaps <- vmatchPattern("-", AA)
	starts <- list()
	maxRF <- max(rFs)
	Ls <- numeric(length(myXStringSet))
	for (i in 1:length(myXStringSet)) {
		start <- start(gaps[[i]])
		if (length(start) > 0) {
			w <- which(start + (length(start) - 1):0==width(AA)[i])
			if (length(w) > 0)
				length(start) <- length(start) - length(w)
			start <- start*3 - 2 + rFs[i]
			start <- sort(c(start, start + 1, start + 2))
			start <- start - 0:(length(start) - 1)
			w <- which(start==(rFs[i] + 1))
			if (length(w) > 0)
				start[w] <- 1
		}
		start <- c(rep(1, maxRF - rFs[i]), start)
		starts[[i]] <- start
		Ls[i] <- length(start) + width(myXStringSet)[i]
	}
	
	# add trailing gaps
	maxWidth <- max(Ls)
	for (i in 1:length(myXStringSet)) {
		starts[[i]] <- c(starts[[i]],
			rep(width(myXStringSet)[i] + 1, maxWidth - Ls[i]))
	}
	
	# remove leading common gaps
	start <- min(unlist(lapply(starts, function(x) return(length(which(x==1))))))
	if (start > 0)
		starts <- lapply(starts, function(x) return(x[-(1:start)]))
	
	myXStringSet <- replaceAt(myXStringSet,
		starts,
		"-")
	if (sense=="-")
		myXStringSet <- reverseComplement(myXStringSet)
	if (direction=="3' to 5'")
		myXStringSet <- reverse(myXStringSet)
	
	return(myXStringSet)
}
