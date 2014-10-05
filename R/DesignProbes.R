`.substr<-` <- function(x,
	start,
	stop,
	value) {
	w <- which(start > stop)
	if (length(w) > 0)
		stop(paste("start", start[w], "> stop", stop[w]))
	w <- which(stop < start)
	if (length(w) > 0)
		stop(paste("stop", stop[w], "< start", start[w]))
	w <- which(start < 1)
	if (length(w) > 0)
		stop(paste("start", start[w], "< 1"))
	w <- which(stop > nchar(x))
	if (length(w) > 0)
		stop(paste("stop", stop[w], "> nchar(x)", nchar(x)[w]))
	
	x <- paste(ifelse(start > 1,
			substr(x,
				1,
				start - 1),
			""),
		value,
		ifelse(stop < nchar(x),
			substr(x,
				stop + 1,
				nchar(x)),
			""),
		sep="")
}

.tiledTemplate <- function(tiles) {
	template <- ""
	begin <- 0
	repeat {
		w <- which(tiles$start > begin)[1]
		
		if (is.na(w)) {
			# complete template with last tile
			last <- which(tiles$start==tiles$start[dim(tiles)[1]])[1]
			if (tiles$end[last]==begin)
				break
			template <- paste(template,
				substr(tiles$target_site[last],
					begin - tiles$start[last] + 2,
					nchar(tiles$target_site[last])),
				sep="")
			dist <- tiles$end[last] - nchar(template)
			if (dist > 0)
				template <- paste(template,
					paste(rep("A", dist), collapse=""),
					sep="")
			break
		} else {
			# fill missing areas of template
			dist <- tiles$start[w] - begin
			if (dist > 1)
				template <- paste(template,
					paste(rep("A", dist - 1), collapse=""),
					sep="")
		}
		
		template <- paste(template,
			tiles$target_site[w], sep="")
		begin <- tiles$start[w] + nchar(tiles$target_site[w]) - 1
	}
	
	return(template)
}

.CalculateEfficiencyFISH <- function(probe,
	template,
	fivePrimeEnd,
	threePrimeEnd,
	temp,
	P,
	ions,
	FA,
	batchSize=1000,
	minEfficiency=0.1,
	type="SSU") {
	
	# error checking
	if (is.character(probe))
		probe <- toupper(probe)
	if (is(probe, "DNAStringSet"))
		probe <- strsplit(toString(probe), ", ", fixed=TRUE)[[1]]
	if (is(template, "DNAStringSet"))
		template <- strsplit(toString(template), ", ", fixed=TRUE)[[1]]
	if (!is.character(probe))
		stop("probe must be a DNAStringSet or character vector.")
	if (!is.character(template))
		stop("template must be a DNAStringSet or character vector.")
	if (!is.numeric(ions))
		stop("ions must be a numeric.")
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
	if (any(!is.numeric(fivePrimeEnd)))
		stop("fivePrimeEnd must be a numeric.")
	if (any(fivePrimeEnd < 1))
		stop("fivePrimeEnd must be greater than zero.")
	if (any(fivePrimeEnd != floor(fivePrimeEnd)))
		stop("fivePrimeEnd must be an integer.")
	if (any(!is.numeric(threePrimeEnd)))
		stop("threePrimeEnd must be a numeric.")
	if (any(threePrimeEnd < 1))
		stop("threePrimeEnd must be greater than zero.")
	if (any(threePrimeEnd != floor(threePrimeEnd)))
		stop("threePrimeEnd must be an integer.")
	if (any(fivePrimeEnd > threePrimeEnd))
		stop("fivePrimeEnd must be less than threePrimeEnd.")
	if (!is.numeric(minEfficiency))
		stop("minEfficiency must be a numeric.")
	if (any(minEfficiency < 0))
		stop("minEfficiency must be greater than or equal to zero.")
	if (any(minEfficiency > 1))
		stop("minEfficiency must be less than or equal to one.")
	
	RT <- .0019871*(273.15 + temp) # [kcal/mol]
	l <- length(probe)
	n <- nchar(probe)
	
	if (l==0)
		stop("No probe specified.")
	if (l!=length(template))
		stop("probe is not the same length as template.")
	if (l!=length(fivePrimeEnd))
		stop("probe is not the same length as fivePrimeEnd.")
	if (l!=length(threePrimeEnd))
		stop("threePrimeEnd is not the same length as fivePrimeEnd.")
	target <- substr(template, fivePrimeEnd, threePrimeEnd)
	
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
	#ws <- p@subject@range@width # width of templates
	
	deltas <- .Call("calculateFISH", seqs1, seqs2, PACKAGE="DECIPHER")
	dG1_PM_DNARNA <- deltas[,1] - (273.15 + temp)/1000*(deltas[,2] - 0.368*n*log(ions))
	
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
		ddG1_MM_DNARNA[MM] <- deltas[MM,3] - (273.15 + temp)/1000*(deltas[MM,4] - 0.368*n[MM]*log(ions))
		ddG1_MM_DNADNA[MM] <- deltas[MM,5] - (273.15 + temp)/1000*(deltas[MM,6] - 0.368*n[MM]*log(ions))
		ddG1_MM_RNARNA[MM] <- deltas[MM,7] - (273.15 + temp)/1000*(deltas[MM,8] - 0.368*n[MM]*log(ions))
		
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
	#cat(dG1, sep=",")
	#cat("\n")
	K1 <- exp(-(dG1 + FA*(0.0696 + 0.0079*n))/RT)
	
	eff <- P*K1/(1 + P*K1)
	
	# disregard probes with < minEfficiency (from dG1 alone)
	z <- which(eff >= minEfficiency)
	l <- length(z)
	if (l != length(eff))
		eff[-z] <- 0
	if (l == 0) {
		ans <- matrix(nrow=length(dG1),
		ncol=6,
		dimnames=list(1:length(dG1),
			c("HybEff",
				"FAm",
				"ddG1",
				"dG1",
				"dG2",
				"dG3")))
		ans[, "HybEff"] <- eff
		ans[, "ddG1"] <- ddG1
		ans[, "dG1"] <- dG1
		return(ans)
	}
	
	probe <- probe[z]
	template <- template[z]
	fivePrimeEnd <- fivePrimeEnd[z]
	threePrimeEnd <- threePrimeEnd[z]
	K1 <- K1[z]
	#ws <- ws[z]
	n <- n[z]
	
	# probe folding
	seqs <- probe
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
	#cat(dG2a - 0.1*FA, sep=",")
	#cat("\n")
	K2a <- exp(-(dG2a + 0.1*FA)/RT)
	
	# target folding
	if (type=="SSU") {
		mapping <- matrix(0, nrow=16, ncol=16)
		dimnames(mapping) <- list(c(names(IUPAC_CODE_MAP), "-"), c(names(IUPAC_CODE_MAP), "-"))
		mapping["A",c("A","M","R","W","V","H","D","N")] <- 1
		mapping["G",c("G","R","S","K","V","D","D","N")] <- 1
		mapping["C",c("C","M","S","Y","V","H","B","N")] <- 1
		mapping["T",c("T","W","Y","K","H","D","B","N")] <- 1
		w <- which(mapping==0)
		mapping[w] <- .1
		mapping[,"-"] <- 0
		ecoli <- "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA"
		
		p <- pairwiseAlignment(unique(template), ecoli, type="global-local", fuzzyMatrix=mapping)
		p <- p[match(template, unique(template))]
		
		# determine domain positioning
		ins <- insertion(p)
		starts <- lapply(ins, slot, "start")
		widths <- lapply(ins, slot, "width")
		
		index <- p@subject@range@start > 567 | (p@subject@range@start + p@subject@range@width) <= 567
		domain2start <- ifelse(index,
			1,
			568 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 568)
			if (length(w) > 0)
				domain2start[j] <- domain2start[j] + cum[length(w)]
		}
		
		index <- p@subject@range@start > 913 | (p@subject@range@start + p@subject@range@width) <= 913
		domain3start <- ifelse(index,
			domain2start,
			914 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 913)
			if (length(w) > 0)
				domain3start[j] <- domain3start[j] + cum[length(w)]
		}
		
		index <- p@subject@range@start > 1397 | (p@subject@range@start + p@subject@range@width) <= 1397
		domain4start <- ifelse(index,
			domain3start,
			1398 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 1398)
			if (length(w) > 0)
				domain4start[j] <- domain4start[j] + cum[length(w)]
		}
		
		# find target site position in the alignment
		startsPattern <- p@pattern@range@start - 1
		del <- deletion(p) # gaps added to pattern
		starts <- fivePrimeEnd
		ends <- threePrimeEnd
		ls <- length(p)
		for (j in 1:ls) {
			d <- del[[j]]
			if (length(d)==0)
				next
			w <- which((d@start + startsPattern[j]) <= starts[j])
			if (length(w) > 0)
				starts[j] <- starts[j] + sum(d@width[w])
			w <- which((d@start + startsPattern[j]) <= ends[j])
			if (length(w) > 0)
				ends[j] <- ends[j] + sum(d@width[w])
		}
		
		# determine start and end of each sequence's domain
		domainStarts <- numeric(ls)
		domainEnds <- numeric(ls)
		domains <- matrix(c(rep(1L, ls),
				domain2start,
				domain3start,
				domain4start,
				nchar(template) - p@pattern@range@width + nchar(unlist(strsplit(toString(pattern(p)), ", ", fixed=TRUE))) + 1),
			nrow=ls)
		for (j in 1:ls) {
			w <- which(domains[j,] <= starts[j])
			domainStarts[j] <- domains[j, w[length(w)]]
			w <- which(domains[j,] > ends[j])
			domainEnds[j] <- domains[j, w[1]] - 1
		}
		
		# determine positioning in the original sequence
		for (j in 1:ls) {
			d <- del[[j]]
			if (length(d)==0)
				next
			w <- which((d@start + startsPattern[j] + cumsum(d@width) - d@width) <= domainStarts[j])
			if (length(w) > 0)
				domainStarts[j] <- domainStarts[j] - sum(d@width[w])
			w <- which((d@start + startsPattern[j] + cumsum(d@width) - d@width) <= domainEnds[j])
			if (length(w) > 0)
				domainEnds[j] <- domainEnds[j] - sum(d@width[w])
		}
		
		batchSize <- floor(batchSize*max(n)/max(domainEnds - domainStarts))
		if (batchSize < 1)
			batchSize <- 1
	} else if (type=="LSU") {
		mapping <- matrix(0, nrow=16, ncol=16)
		dimnames(mapping) <- list(c(names(IUPAC_CODE_MAP), "-"), c(names(IUPAC_CODE_MAP), "-"))
		mapping["A",c("A","M","R","W","V","H","D","N")] <- 1
		mapping["G",c("G","R","S","K","V","D","D","N")] <- 1
		mapping["C",c("C","M","S","Y","V","H","B","N")] <- 1
		mapping["T",c("T","W","Y","K","H","D","B","N")] <- 1
		w <- which(mapping==0)
		mapping[w] <- .1
		mapping[,"-"] <- 0
		ecoli <- "GGTTAAGCGACTAAGCGTACACGGTGGATGCCCTGGCAGTCAGAGGCGATGAAGGACGTGCTAATCTGCGATAAGCGTCGGTAAGGTGATATGAACCGTTATAACCGGCGATTTCCGAATGGGGAAACCCAGTGTGATTCGTCACACTATCATTAACTGAATCCATAGGTTAATGAGGCGAACCGGGGGAACTGAAACATCTAAGTACCCCGAGGAAAAGAAATCAACCGAGATTCCCCCAGTAGCGGCGAGCGAACGGGGAGGAGCCCAGAGCCTGAATCAGTGTGTGTGTTAGTGGAAGCGTCTGGAAAGGCGCGCGATACAGGGTGACAGCCCCGTACACAAAAATGCACATACTGTGAGCTCGATGAGTAGGGCGGGACACGTGGTATCCTGTCTGAATATGGGGGGACCATCCTCCAAGGCTAAATACTCCTGACTGACCGATAGTGAACCAGTACCGTGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGTGAAAAAGAACCTGAAACCGTGTACGTACAAGCAGTGGGAGCCTCTTTTATGGGGTGACTGCGTACCTTTTGTATAATGGGTCAGCGACTTATATTCTGTAGCAAGGTTAACCGAATAGGGGAGCCGAAGGGAAACCGAGTCTTAACCGGGCGTTAAGTTGCAGGGTATAGACCCGAAACCCGGTGATCTAGCCATGGGCAGGTTGAAGGTTGGGTAACACTAACTGGAGGACCGAACCGACTAATGTTGAAAAATTAGCGGATGACTTGTGGCTGGGGGTGAAAGGCCAATCAAACCGGGAGATAGCTGGTTCTCCCCGAAAGCTATTTAGGTAGCGCCTCGTGAATTCATCTCCGGGGGTAGAGCACTGTTTCGGCAAGGGGGTCATCCCGACTTACCAACCCGATGCAAACTGCGAATACCGGAGAATGTTATCACGGGAGACATACGGCGGGTGCTAACGTCCGTCGTGAAGAGGGAAACAACCCAGACCGCCAGCTAAGGTCCCAAAGTCATGGTTAAGTGGGAAACGATGTGGGAAGGCCCAGACAGCCAGGATGTTGGCTTAGAAGCAGCCATCATTTAAAGAAAGCGTAATAGCTCACTGGTCGAGTCGGCCTGCGCGGAAGATGTAACGGGGCTAAACCATGCACCGAAGCTGCGGCAGCGACACTGTGTGTTGTTGGGTAGGGGAGCGTTCTGTAAGCCTGTGAAGGTGTACTGTGAGGTATGCTGGAGGTATCAGAAGTGCGAATGCTGACATAAGTAACGATAAAGCGGGTGAAAAGCCCGCTCGCCGGAAGACCAAGGGTTCCTGTCCAACGTTAATCGGGGCAGGGTGAGTCGACCCCTAAGGCGAGGCCGAAAGGCGTAGTCGATGGGAAACAGGTTAATATTCCTGTACTTGGTGTTACTGCGAAGGGGGGACGGAGAAGGCTATGTTGGCCGGGCGACGGTTGTCCCGGTTTAAGCGTGTAGGCTGGTTTTCCAGGCAAATCCGGAAAATCAAGGCTGAGGCGTGATGACGAGGCACTACGGTGCTGAAGCAACAAATGCCCTGCTTCCAGGAAAAGCCTCTAAGCATCAGGTAACATCAAATCGTACCCCAAACCGACACAGGTGGTCAGGTAGAGAATACCAAGGCGCTTGAGAGAACTCGGGTGAAGGAACTAGGCAAAATGGTGCCGTAACTTCGGGAGAAGGCACGCTGATATGTAGGTGAAGTCCCTCGCGGATGGAGCTGAAATCAGTCGAAGATACCAGCTGGCTGCAACTGTTTATTAAAAACACAGCACTGTGCAAACACGAAAGTGGACGTATACGGTGTGACGCCTGCCCGGTGCCGGAAGGTTAATTGATGGGGTCAGCGCAAGCGAAGCTCTTGATCGAAGCCCCGGTAAACGGCGGCCGTAACTATAACGGTCCTAAGGTAGCGAAATTCCTTGTCGGGTAAGTTCCGACCTGCACGAATGGCGTAATGATGGCCAGGCTGTCTCCACCCGAGACTCAGTGAAATTGAACTCGCTGTGAAGATGCAGTGTACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACTGAACATTGAGCCTTGATGTGTAGGATAGGTGGGAGGCTTTGAAGTGTGGACGCCAGTCTGCATGGAGCCGACCTTGAAATACCACCCTTTAATGTTTGATGTTCTAACGTGGACCCGTGATCCGGGTTGCGGACAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCTCCTCCTAAAGAGTAACGGAGGAGCACGAAGGTTGGCTAATCCTGGTCGGACATCAGGAGGTTAGTGCAATGGCATAAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGTGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGGCGCTGGAGAACTGAGGGGGGCTGCTCCTAGTACGAGAGGACCGGAGTGGACGCATCACTGGTGTTCGGGTTGTCATGCCAATGGCACTGCCCGGTAGCTAAATGCGGAAGAGATAAGTGCTGAAAGCATCTAAGCACGAAACTTGCCCCGAGATGAGTTCTCCCTGACCCTTTAAGGGTCCTGAAGGAACGTTGAAGACGACGACGTTGATAGGCCGGGTGTGTAAGCGCAGCGATGCGTTGAGCTAACCGGTACTAATGAACCGTGAGGCTTAACCTT"
		
		p <- pairwiseAlignment(unique(template), ecoli, type="global-local", fuzzyMatrix=mapping)
		p <- p[match(template, unique(template))]
		
		# determine domain positioning
		ins <- insertion(p)
		starts <- lapply(ins, slot, "start")
		widths <- lapply(ins, slot, "width")
		
		index <- p@subject@range@start > 562 | (p@subject@range@start + p@subject@range@width) <= 562
		domain2start <- ifelse(index,
			1,
			563 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 563)
			if (length(w) > 0)
				domain2start[j] <- domain2start[j] + cum[length(w)]
		}
		
		index <- p@subject@range@start > 1270 | (p@subject@range@start + p@subject@range@width) <= 1270
		domain3start <- ifelse(index,
			domain2start,
			1271 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 1271)
			if (length(w) > 0)
				domain3start[j] <- domain3start[j] + cum[length(w)]
		}
		
		index <- p@subject@range@start > 1647 | (p@subject@range@start + p@subject@range@width) <= 1647
		domain4start <- ifelse(index,
			domain3start,
			1648 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 1648)
			if (length(w) > 0)
				domain4start[j] <- domain4start[j] + cum[length(w)]
		}
		
		index <- p@subject@range@start > 2015 | (p@subject@range@start + p@subject@range@width) <= 2015
		domain5start <- ifelse(index,
			domain4start,
			2016 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 2016)
			if (length(w) > 0)
				domain5start[j] <- domain5start[j] + cum[length(w)]
		}
		
		index <- p@subject@range@start > 2626 | (p@subject@range@start + p@subject@range@width) <= 2626
		domain6start <- ifelse(index,
			domain5start,
			2627 - p@subject@range@start)
		index <- !index
		for (j in which(index)) {
			i <- starts[[j]] + p@subject@range@start[j] - 1
			cum <- cumsum(widths[[j]])
			w <- which(i <= 2627)
			if (length(w) > 0)
				domain6start[j] <- domain6start[j] + cum[length(w)]
		}
		
		# find target site position in the alignment
		startsPattern <- p@pattern@range@start - 1
		del <- deletion(p) # gaps added to pattern
		starts <- fivePrimeEnd
		ends <- threePrimeEnd
		ls <- length(p)
		for (j in 1:ls) {
			d <- del[[j]]
			if (length(d)==0)
				next
			w <- which((d@start + startsPattern[j]) <= starts[j])
			if (length(w) > 0)
				starts[j] <- starts[j] + sum(d@width[w])
			w <- which((d@start + startsPattern[j]) <= ends[j])
			if (length(w) > 0)
				ends[j] <- ends[j] + sum(d@width[w])
		}
		
		# determine start and end of each sequence's domain
		domainStarts <- numeric(ls)
		domainEnds <- numeric(ls)
		domains <- matrix(c(rep(1L, ls),
				domain2start,
				domain3start,
				domain4start,
				domain5start,
				domain6start,
				nchar(template) - p@pattern@range@width + nchar(unlist(strsplit(toString(pattern(p)), ", ", fixed=TRUE))) + 1),
			nrow=ls)
		for (j in 1:ls) {
			w <- which(domains[j,] <= starts[j])
			domainStarts[j] <- domains[j, w[length(w)]]
			w <- which(domains[j,] > ends[j])
			domainEnds[j] <- domains[j, w[1]] - 1
		}
		
		# determine positioning in the original sequence
		for (j in 1:ls) {
			d <- del[[j]]
			if (length(d)==0)
				next
			w <- which((d@start + startsPattern[j] + cumsum(d@width) - d@width) <= domainStarts[j])
			if (length(w) > 0)
				domainStarts[j] <- domainStarts[j] - sum(d@width[w])
			w <- which((d@start + startsPattern[j] + cumsum(d@width) - d@width) <= domainEnds[j])
			if (length(w) > 0)
				domainEnds[j] <- domainEnds[j] - sum(d@width[w])
		}
		
		batchSize <- floor(batchSize*max(n)/max(domainEnds - domainStarts))
		if (batchSize < 1)
			batchSize <- 1
	} else { # domain is type nucleotides to either side of the type site
		type <- as.numeric(type)
		batchSize <- floor(batchSize*max(n)/(max(n) + 2*type))
		if (batchSize < 1)
			batchSize <- 1
		
		domainStarts <- ifelse(fivePrimeEnd - type < 1,
			1,
			fivePrimeEnd - type)
		domainEnds <- ifelse(threePrimeEnd + type > nchar(template),
			nchar(template),
			threePrimeEnd + type)
	}
	
	seqs <- substr(template, domainStarts, domainEnds)
	seq <- unique(seqs)
	ls <- length(seq)
	dG3_folded <- numeric(ls)
	for (start in seq(1, ls, batchSize)) {
		end <- ifelse(start + batchSize > ls, ls, start + batchSize)
		dG3_folded[start:end] <- as.numeric(system(paste("hybrid-ss-min -n RNA -t",
				temp,
				"-T",
				temp,
				"-E -q",
				paste(seq[start:end], collapse=" ")),
			intern=TRUE))
	}
	dG3 <- numeric(l)
	dG3 <- dG3_folded[match(seqs, seq)]
	
	#seqs <- substr(template, domainStarts, domainEnds - 1)
	prohibit <- fivePrimeEnd - domainStarts + 1
	groups <- paste(prohibit, "0", threePrimeEnd - fivePrimeEnd + 1, sep=",")
	g <- unique(groups)
	dG3_unfolded <- numeric(l)
	for (i in 1:length(g)) {
		w <- which(groups == g[i])
		seq <- unique(seqs[w])
		ls <- length(seq)
		for (start in seq(1, ls, batchSize)) {
			end <- ifelse(start + batchSize > ls, ls, start + batchSize)
			m <- match(seqs[w], seq[start:end])
			na <- which(!is.na(m))
			dG3_unfolded[w][na] <- as.numeric(system(paste("hybrid-ss-min -n RNA -t",
					temp,
					"-T",
					temp,
					paste("--prohibit=", g[i], sep=""),
					"-E -q",
					paste(seq[start:end], collapse=" ")),
				intern=TRUE))[m][na]
		}
	}
	dG3_unfolded <- ifelse(dG3_unfolded > 0, 0, dG3_unfolded)
	dG3 <- dG3 - dG3_unfolded
	#cat(dG3, sep=",")
	#cat("\n")
	K3 <- exp(-(dG3 + FA*-0.0111*dG3)/RT)
	
	Kov <- K1/((1 + K2a)*(1 + K3))
	eff[z] <- P*Kov/(1 + P*Kov)
	
	FAm <- numeric(length(z))
	FAm[] <- NA
	range <- -1000:1000
	for (i in 1:length(z)) {
		f <- function(FA) {
			K1 <- exp(-(dG1[z][i] + FA*(0.0696 + 0.0079*n[i]))/RT)
			K2a <- exp(-(dG2a[i] + 0.1*FA)/RT)
			K3 <- exp(-(dG3[i] + FA*-0.0111*dG3[i])/RT)
			Kov <- K1/((1 + K2a)*(1 + K3))
			return(0.5 - P*Kov/(1 + P*Kov))
		}
		HybEffs <- f(range)
		try(FAm[i] <- uniroot(f,
				c(range[which.min(HybEffs)],
					range[which.max(HybEffs)] + 1))$root,
			silent=TRUE)
	}
	
	ans <- matrix(nrow=length(dG1),
		ncol=7,
		dimnames=list(1:length(dG1),
			c("HybEff",
				"FAm",
				"ddG1",
				"dG1",
				"dG2",
				"dG3",
				"dGov")))
	ans[, "HybEff"] <- eff
	ans[z, "FAm"] <- FAm
	ans[, "ddG1"] <- ddG1
	ans[, "dG1"] <- dG1
	ans[z, "dG2"] <- dG2a
	ans[z, "dG3"] <- dG3
	ans[z, "dGov"] <- -RT*log(Kov)
	
	return(ans)
}

DesignProbes <- function(tiles,
	identifier="",
	start=1,
	end=NULL,
	minLength=17,
	maxLength=26,
	maxPermutations=4,
	minCoverage=0.9,
	minGroupCoverage=0.2,
	hybTemp=46,
	P=250e-9,
	Na=1,
	FA=35,
	minEfficiency=0.50,
	worstScore=-Inf,
	numProbeSets=0,
	batchSize=1000,
	target="SSU",
	verbose=TRUE) {
	
	# error checking
	TARGETS <- c("SSU","LSU","Other")
	molecule <- pmatch(target, TARGETS)
	if (is.na(molecule))
		stop("Invalid target.  Choose either SSU, LSU, or Other.")
	if (molecule == -1)
		stop("Ambiguous target.  Choose either SSU, LSU, or Other.")
	molecule <- c("SSU", "LSU", 200)[molecule]
	if (!is.numeric(Na))
		stop("[Na+] concentration must be a numeric.")
	if (Na < .01 || is.nan(Na))
		stop("Sodium equivilent concentration must be at least 0.01M.")
	if (Na > 1)
		stop("Sodium equivilent concentration can be at most 1.0M.")
	if (minLength > maxLength)
		stop("minLength must be less or equal to maxLength.")
	ids <- unique(tiles$id)
	if (length(ids)==0)
		stop("No identifiers found in tiles.")
	if (any(grepl("mol,", ids, fixed=TRUE)))
		stop("Identifiers in tiles cannot include 'mol,'.")
	if (any(grepl("%;", ids, fixed=TRUE)))
		stop("Identifiers in tiles cannot include '%;'.")
	if (any(grepl(" (", ids, fixed=TRUE)))
		stop("Identifiers in tiles cannot include ' ('.")
	if (any(grepl(") ", ids, fixed=TRUE)))
		stop("Identifiers in tiles cannot include ') '.")
	if (!is.numeric(hybTemp))
		stop("hybTemp must be a numeric.")
	if (hybTemp < 10)
		stop("hybTemp must be at least 10 degrees Celsius.")
	if (hybTemp >= 90)
		stop("hybTemp must be less than 90 degrees Celsius.")
	if (!is.numeric(FA))
		stop("FA must be a numeric.")
	if (FA < 0)
		stop("FA must be at least 0% (volume/volume).")
	if (FA > 100)
		stop("FA can be at most 100% (volume/volume).")
	if (FA > 50 | FA < 20)
		warning("FA between 20 and 50% (volume/volume) will give better results.")
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
	if (!is.numeric(numProbeSets))
		stop("numProbeSets must be a numeric.")
	if (floor(numProbeSets)!=numProbeSets)
		stop("numProbeSets must be a whole number.")
	if (numProbeSets < 0)
		stop("numProbeSets must be greater than or equal to zero.")
	if (!is.numeric(maxPermutations))
		stop("maxPermutations must be a numeric.")
	if (floor(maxPermutations)!=maxPermutations)
		stop("maxPermutations must be a whole number.")
	if (maxPermutations <= 0)
		stop("maxPermutations must be greater than zero.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
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
	if (class(try(system("hybrid-min -V", intern=TRUE), silent=TRUE))=="try-error")
		stop("OligoArrayAux must be properly installed.  See the Note section in the help page for DesignProbes (enter: ?DesignProbes).")
	
	probes_all <- data.frame()
	
	if (numProbeSets > 0) {
		l <- length(identifier)
		pSets <- data.frame(identifier=I(rep(identifier, each=numProbeSets)),
			start_one=I(integer(l*numProbeSets)),
			start_two=I(integer(l*numProbeSets)),
			start_aligned_one=I(integer(l*numProbeSets)),
			start_aligned_two=I(integer(l*numProbeSets)),
			permutations_one=I(integer(l*numProbeSets)),
			permutations_two=I(integer(l*numProbeSets)),
			score_one=I(numeric(l*numProbeSets)),
			score_two=I(numeric(l*numProbeSets)),
			score_set=I(numeric(l*numProbeSets)),
			one_probe=I(matrix(character(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			two_probe=I(matrix(character(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			one_efficiency=I(matrix(numeric(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			two_efficiency=I(matrix(numeric(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			one_FAm=I(matrix(numeric(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			two_FAm=I(matrix(numeric(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			one_coverage=I(matrix(numeric(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			two_coverage=I(matrix(numeric(),
				nrow=l*numProbeSets,
				ncol=maxPermutations)),
			mismatches_one=I(character(l*numProbeSets)),
			mismatches_two=I(character(l*numProbeSets)),
			mismatches_set=I(character(l*numProbeSets)))
	}
	l_ids <- length(ids)
	for (id in identifier) {
		if (verbose) {
			time.1 <- Sys.time()
			cat("\n", id, sep="")
			flush.console()
		}
		
		w <- which(tiles$id==id)
		
		if (length(w)==0) {
			warning("No target sites met the specified constraints: ", id)
			next
		}
		
		target <- tiles[w,]
		template <- .tiledTemplate(target)
		
		if (is.null(end)) {
			w <- which(tiles$id==id &
				tiles$start_aligned >= start)
		} else {
			w <- which(tiles$id==id &
				tiles$start_aligned >= start &
				tiles$end_aligned <= end)
		}
		if (length(w)==0) {
			warning("No target sites met the specified constraints: ", id)
			next
		}
		
		target <- tiles[w,]
		
		uw <- unique(target$width)
		if (length(uw) > 1) {
			warning("Multiple sequence widths: ", id)
			next
		}
		
		u <- unique(target$start_aligned)
		l <- length(u)
		probes <- data.frame(identifier=I(rep(id,l)),
			start=I(integer(l)),
			start_aligned=I(integer(l)),
			permutations=I(integer(l)),
			score=I(numeric(l)),
			probe=I(matrix(character(),
				nrow=l,
				ncol=maxPermutations)),
			efficiency=I(matrix(numeric(),
				nrow=l,
				ncol=maxPermutations)),
			FAm=I(matrix(numeric(),
				nrow=l,
				ncol=maxPermutations)),
			coverage=I(matrix(numeric(),
				nrow=l,
				ncol=maxPermutations)),
			mismatches=I(character(l)))
		
		targets <- array("",
			c(l, maxPermutations, maxLength - minLength + 1),
			dimnames=list(position=NULL, permutation=NULL, length=NULL))
		
		begin <- 1L
		for (i in 1:l) {
			# find all target_sites
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
			
			probes$permutations[i] <- perms
			probes$coverage[i, 1:perms] <- target$coverage[w]
			
			# determine optimal length
			for (j in 1:perms) {
				n <- nchar(target_site[j])
				for (p in minLength:ifelse(n > maxLength, maxLength, n)) {
					targets[i, j, p - minLength + 1] <- substr(target_site[j], 1, p)
				}
			}
			probes$start_aligned[i] <- u[i]
			probes$start[i] <- target$start[w[1]]
		}
		
		w <- which(probes$permutations==0)
		if (length(w)==dim(probes)[1]) {
			warning("All target sites have too many permutations: ", id)
			next
		}
		if (length(w) > 0) {
			probes <- probes[-w,]
			targets <- targets[-w,,]
			dim(targets) <- c(l - length(w), maxPermutations, maxLength - minLength + 1)
			l <- dim(probes)[1]
		}
		
		# Choose optimal length for each probe
		w <- integer()
		empty <- integer()
		for (j in 1:(maxLength - minLength + 1)) {
			if (length(empty) > 0) {
				empty <- c((1:(l*maxPermutations))[w],
					which(targets[,,j]==""),
					empty)
			} else {
				empty <- c((1:(l*maxPermutations))[w],
					which(targets[,,j]==""))
			}
			t <- targets[,,j]
			if (length(empty) > 0) {
				t <- t[-empty]
			} else {
				t <- as.vector(t)
			}
			if (length(t)==0)
				break
			
			# calculate hybridization efficiency
			eff <- CalculateEfficiencyFISH(unlist(strsplit(toString(reverseComplement(DNAStringSet(t))),
					", ",
					fixed=TRUE)),
				t,
				hybTemp,
				P,
				Na,
				FA,
				batchSize)
			w <- which(eff[, "HybEff"] >= minEfficiency)
			if (length(w)==0)
				next
			
			if (length(empty) > 0) {
				probes$probe[-empty][w] <- strsplit(toString(reverseComplement(DNAStringSet(targets[,,j][-empty][w]))),
					", ",
					fixed=TRUE)[[1]]
				probes$efficiency[-empty][w] <- eff[w, "HybEff"]
				probes$FAm[-empty][w] <- eff[w, "FAm"]
			} else {
				probes$probe[w] <- strsplit(toString(reverseComplement(DNAStringSet(targets[,,j][w]))),
					", ",
					fixed=TRUE)[[1]]
				probes$efficiency[w] <- eff[w, "HybEff"]
				probes$FAm[w] <- eff[w, "FAm"]
			}
			
			a <- apply(probes$probe, 1, function(x) length(which(x != "")))
			w <- which(probes$permutations != a)
			if (length(w) > 0) {
				probes$probe[w,] <- NA
				probes$efficiency[w,] <- NA
				probes$FAm[w,] <- NA
			}
			
			w <- which(probes$probe != "")
		}
		
		w <- integer()
		for (i in 1:l) {
			if (probes$start_aligned[i]==0) {
				w <- c(w,i) # empty probe
				next
			}
			p <- probes$permutations[i]
			
			numPerms <- Inf
			full <- which(probes$probe[i,]!="")
			if (length(full) != 0) {
				c <- suppressWarnings(ConsensusSequence(DNAStringSet(probes$probe[i, full]),
					threshold=0.01,
					verbose=FALSE))
				a <- alphabetFrequency(c)
				numPerms <- sum(2*a[1, 5:10], 3*a[1, 11:14], 4*a[1, 15])
			}
			
			if (numPerms > probes$permutations[i] || p != length(which(probes$probe[i,]!=""))) {
				w <- c(w,i)
				probes$permutations[i] <- 0
				probes$probe[i,] <- NA
				probes$efficiency[i,] <- NA
				probes$score[i] <- -Inf
			} else {
				d <- duplicated(probes$probe[i,], incomparables=NA)
				l <- length(which(d))
				if (l > 0) {
					for (j in which(d)) {
						w1 <- which(probes$probe[i,1:(j - 1)]==probes$probe[i,j])[1]
						probes$coverage[i,w1] <- probes$coverage[i,w1] + probes$coverage[i,j]
						probes$coverage[i,j] <- NA
					}
					probes$permutations[i] <- probes$permutations[i] - l
					probes$probe[i,] <- c(probes$probe[i,!d],
						rep(NA, l))
					probes$efficiency[i,] <- c(probes$efficiency[i,!d],
						rep(NA, l))
					probes$FAm[i,] <- c(probes$FAm[i,!d],
						rep(NA, l))
					probes$coverage[i,] <- c(probes$coverage[i,!d],
						rep(NA, l))
				}
			}
		}
		
		if (length(w) > 0)
			probes <- probes[-w,]
		d <- dim(probes)[1]
		if (d==0) {
			warning("No probes met the specified constraints: ",id)
			next
		}
		
		# Determine which probes also meet original model constraints
		w <- which(probes$probe != "")
		index <- rep(1:d, maxPermutations)
		index <- index[w]
		
		templates <- rep(template, length(w))
		.substr(templates,
			probes$start[index],
			probes$start[index] + nchar(probes$probe[w]) - 1) <-
			unlist(strsplit(toString(reverseComplement(DNAStringSet(probes$probe[w]))),
				", ",
				fixed=TRUE))
		eff <- .CalculateEfficiencyFISH(as.character(probes$probe[w]),
			templates,
			probes$start[index],
			probes$start[index] + nchar(probes$probe[w]) - 1,
			hybTemp,
			P,
			Na,
			FA=0,
			batchSize,
			type=molecule)
		
		w <- which(!is.na(eff[,"FAm"]))
		if (length(w) > 0)
			w <- w[which(eff[w, "FAm"] < 70 &
				eff[w, "FAm"] > (FA - 15) &
				eff[w, "dG2"] > -3)]
		if (length(w)==0) {
			w <- 1:(dim(probes)[1])
		} else {
			w <- unique(index[-w])
		}
		
		if (length(w) > 0)
			probes <- probes[-w,]
		d <- dim(probes)[1]
		if (d==0) {
			warning("No probes met the specified constraints: ",id)
			next
		}
		
		if (verbose) {
			cat(" (", d, "):\n", sep="")
			flush.console()
			pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
		}
		
		for (k in 1:length(ids)) {
			if (ids[k] == id)
				next
			
			w <- which(tiles$id==ids[k])
			target <- tiles[w,]
			
			if (any(uw != unique(target$width)))
				next
			
			combos <- numeric(d)
			count <- 1
			nontargets <- targets <- character(d*sum(probes$permutations)*max(table(target$start)))
			FAms <- numeric(d*sum(probes$permutations)*max(table(target$start)))
			
			begin <- 1L
			for (i in 1:d) {
				# find all non-target_sites
				w <- .Call("multiMatchUpper", target$start_aligned, probes$start_aligned[i], begin, PACKAGE="DECIPHER")
				
				if (length(w)==0) # target does not overlap nontarget
					break
				
				begin <- w[1]
				target_site <- target$target_site[w]
				
				j <- length(target_site)
				if (j==0) # no nontarget
					next
				
				w <- .Call("multiMatchCharNotNA", probes$probe[i,], PACKAGE="DECIPHER")
				l <- length(w)
				
				if (l > 0) {
					combos[i] <- j*l
					targets[count:(count + combos[i] - 1)] <- rep(probes$probe[i,w],
						each=j)
					nontargets[count:(count + combos[i] - 1)] <- rep(target_site,
						l)
#					FAms[count:(count + combos[i] - 1)] <- rep(probes$FAm[i,w],
#						each=j)
					FAms[count:(count + combos[i] - 1)] <- min(probes$FAm[i,w])
					count <- count + combos[i]
				}
			}
			
			if (count==1) # sites do not overlap
				next
			
#			p <- pairwiseAlignment(nontargets[1:(count - 1)],
#				reverseComplement(DNAStringSet(targets[1:(count - 1)])),
#				type="local-global",
#				gapOpen=-10,
#				gapExtension=-10)
			
			eff <- CalculateEfficiencyFISH(targets[1:(count - 1)],
				nontargets[1:(count - 1)],
				hybTemp,
				P,
				Na,
				FA,
				batchSize)
			
			eff_max <- numeric(d)
			ddG1 <- numeric(d)
			FAm <- numeric(d)
			dFAms <- numeric(d)
			dFAms[] <- -Inf
			index <- integer(d)
			count <- 1
			for (i in 1:d) {
				if (combos[i] != 0) {
					dFAm <- eff[count:(count + combos[i] - 1), "FAm"] - FAms[count:(count + combos[i] - 1)]
					w <- which.max(dFAm)
					dFAms[i] <- dFAm[w]
					eff_max[i] <- eff[count:(count + combos[i] - 1), "HybEff"][w]
					ddG1[i] <- eff[count:(count + combos[i] - 1), "ddG1"][w]
					FAm[i] <- eff[count:(count + combos[i] - 1), "FAm"][w]
					index[i] <- count + w - 1
					count <- count + combos[i]
				}
			}
			
			probes$score <- probes$score - ifelse(dFAms < -20, 0, 0.2 + 1.2^ifelse(dFAms > 0, 0, dFAms)) # score = 0.2*n + 1.2^dFAm
			
			w <- which(dFAms >= -20)
			
			if (length(w)==0)
				next
			
			s1 <- targets[index[w]]
			m <- match(s1, probes$probe)
			s2 <- reverseComplement(DNAStringSet(nontargets[index[w]]))
			p <- pairwiseAlignment(s1,
				paste("----", s2, "----", sep=""),
				type="global-local",
				gapOpen=-10,
				gapExtension=-10)
			s1 <- toString(pattern(p))
			s1 <- unlist(strsplit(s1, ", ", fixed=TRUE))
			s2 <- toString(subject(p))
			s2 <- unlist(strsplit(s2, ", ", fixed=TRUE))
			s2 <- toString(reverseComplement(DNAStringSet(s2)))
			s2 <- unlist(strsplit(s2, ", ", fixed=TRUE))
			
			probes$mismatches[w] <- paste(probes$mismatches[w],
				ids[k],
				" (",
				formatC(100*eff_max[w],
					digits=1,
					width=1,
					format="f"),
				"%,",
				formatC(ddG1[w],
					digits=2,
					width=1,
					format="f"),
				"kcal/mol,",
				formatC(dFAms[w],
					digits=1,
					width=1,
					format="f"),
				"%;",
				s1,
				"/",
				s2,
				") ",
				sep="")
			
			if (worstScore > -Inf) {
				w <- probes$score < worstScore
				if (any(w)) {
					probes$probe[w,] <- NA
					probes$efficiency[w,] <- NA
					probes$FAm[w,] <- NA
					probes$mismatches[w] <- ""
					probes$score[w] <- -Inf
					probes$permutations[w] <- 0
					probes$coverage[w,] <- NA
				}
				w <- which(w)
				if (length(w) > 0)
					probes <- probes[-w,]
				d <- dim(probes)[1]
				if (d==0) {
					warning("All probes were below the worstScore: ",id)
					next
				}
			}
			
			if (verbose)
				setTxtProgressBar(pBar,
					floor(100*k/l_ids))
		}
		
		if (numProbeSets > 0) {
			if (d==1) {
				warning("Not enough probes to design a probe set:", id)
				break
			}
			# order probes by score
			searchSpace <- ifelse(numProbeSets > 100, numProbeSets, 100)
			f <- which(is.finite(probes$score))
			n <- unlist(lapply(strsplit(probes$mismatches, ") ", fixed=TRUE), length))
			o <- order(-1*n[f],
				probes$score[f],
				-1*probes$permutations[f],
				rowSums(as.matrix(probes$coverage[f,]), na.rm=TRUE),
				decreasing=TRUE)
			s <- ifelse(length(f) > searchSpace, searchSpace, length(f))
			if (s < 2) { # need at least 2 probes to form a set
				warning("Not enough probe sets meet the specified constaints:", id)
				
				if (verbose) {
					setTxtProgressBar(pBar, 100)
					close(pBar)
					time.2 <- Sys.time()
					cat("\n")
					print(signif(difftime(time.2,
						time.1,
						units='secs'),
						digits=2))
				}
				
				next
			}
			
			# calculate the set's efficiency
			MMs <- strsplit(probes$mismatches[f[o][1:s]], "mol,", fixed=TRUE)
			ls <- unlist(lapply(MMs, length))
			ls <- ifelse(ls > 0, ls - 1, 0)
			index <- rep(1:length(ls), ls)
			if (length(index) > 0) {
				MMs <- strsplit(unlist(MMs), "%;", fixed=TRUE)
				MMs <- unlist(strsplit(unlist(MMs), ") ", fixed=TRUE))
				effs <- as.numeric(MMs[seq(2, length(MMs), 3)])
				MMs <- MMs[seq(1, length(MMs), 3)]
				MMs <- unlist(strsplit(MMs, " (", fixed=TRUE))
				MMs <- MMs[seq(1, length(MMs), 2)]
			}
			
			m <- matrix(Inf,
				nrow=s,
				ncol=s,
				dimnames=list(f[o][1:s],
					f[o][1:s]))
			n <- matrix(Inf,
				nrow=s,
				ncol=s,
				dimnames=list(f[o][1:s],
					f[o][1:s]))
			for (i in 1:(s - 1)) {
				w1 <- which(index==i)
				if (length(w1)==0) {
					n[i, (i + 1):s] <- 0
					m[i, (i + 1):s] <- 0
					next
				}
				for (j in (i + 1):s) {
					w2 <- which(index==j)
					if (length(w2)==0) {
						n[i, j] <- 0
						m[i, j] <- 0
						next
					}
					overlap_MMs <- match(MMs[w1], MMs[w2])
					w <- which(!is.na(overlap_MMs))
					n[i, j] <- length(w)
					if (length(w) > 0) {
						m[i, j] <- sum(tapply(c(effs[w1[w]], effs[w2[overlap_MMs[w]]]),
							rep(1:length(w), 2),
							function(x) return(1.2^ifelse(min(x) > 0, 0, min(x)))))
					} else {
						m[i, j] <- 0
					}
				}
			}
			
			# choose the best probe sets
			d <- dimnames(m)
			w <- which(pSets$identifier==id)
			p <- outer(probes$permutations[f[o][1:s]],
				probes$permutations[f[o][1:s]],
				FUN="+")
			c <- -1*outer(rowSums(as.matrix(probes$coverage[f[o][1:s],]), na.rm=TRUE),
				ifelse(rep(s, s)==1,
					sum(probes$coverage[f[o][1:s],], na.rm=TRUE),
					rowSums(as.matrix(probes$coverage[f[o][1:s],]), na.rm=TRUE)),
				FUN="*")
			ss <- -1*outer(probes$score[f[o][1:s]],
				probes$score[f[o][1:s]],
				FUN="+")
			o <- order(m + n/5, p, c, ss) # score = 0.2*n + 1.2^dFAm
			j <- 0
			dims <- dim(m)
			for (k in 1:numProbeSets) {
				if ((j + 1) > length(o)) {
					warning("Not enough probe sets meet the specified constaints:", id)
					break
				}
				for (j in (j + 1):length(o)) {
					w_o <- c((o[j] - 1)%%dims[1] + 1,
						(o[j] - 1)%/%dims[1] + 1)
					f <- as.numeric(d[[1]][w_o[1]])
					r <- as.numeric(d[[2]][w_o[2]])
					
					start_F <- probes$start[f]
					start_R <- probes$start[r]
					if (abs(start_F - start_R) > (maxLength + 50)) # 50 nt separation
						break
				}
				
				if (abs(start_F - start_R) > maxLength) {
					pSets$start_one[w[k]] <- probes$start[f]
					pSets$start_two[w[k]] <- probes$start[r]
					pSets$start_aligned_one[w[k]] <- probes$start_aligned[f]
					pSets$start_aligned_two[w[k]] <- probes$start_aligned[r]
					pSets$permutations_one[w[k]] <- probes$permutations[f]
					pSets$permutations_two[w[k]] <- probes$permutations[r]
					pSets$score_one[w[k]] <- probes$score[f]
					pSets$score_two[w[k]] <- probes$score[r]
					pSets$score_set[w[k]] <- -m[o[j]]/100
					pSets$one_probe[w[k],] <- probes$probe[f,]
					pSets$two_probe[w[k],] <- probes$probe[r,]
					pSets$one_efficiency[w[k],] <- probes$efficiency[f,]
					pSets$two_efficiency[w[k],] <- probes$efficiency[r,]
					pSets$one_FAm[w[k],] <- probes$FAm[f,]
					pSets$two_FAm[w[k],] <- probes$FAm[r,]
					pSets$one_coverage[w[k],] <- probes$coverage[f,]
					pSets$two_coverage[w[k],] <- probes$coverage[r,]
					pSets$mismatches_one[w[k]] <- probes$mismatches[f]
					pSets$mismatches_two[w[k]] <- probes$mismatches[r]
					
					w_F <- which(index==w_o[1])
					w_R <- which(index==w_o[2])
					overlap_MMs <- match(MMs[w_F], MMs[w_R])
					w_S <- which(!is.na(overlap_MMs))
					if (length(w_S)==0)
						next
					EFFs <- tapply(c(effs[w_F[w_S]], effs[w_R[overlap_MMs[w_S]]]),
							rep(1:length(w_S), 2),
							min)
					# record set mismatches
					w1 <- which(EFFs >= -20)
					if (length(w1) > 0)
						pSets$mismatches_set[w[k]] <- paste(MMs[w_F][w_S][w1],
							" (",
							formatC(EFFs[w1],
								digits=1,
								width=1,
								format="f"),
							"%)",
							sep="",
							collapse=", ")
				} else {
					warning("Not enough probes to meet the specified constraints: ", id)
					break
				}
			}
		}
		
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			time.2 <- Sys.time()
			cat("\n")
			print(signif(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
		}
		
		if (numProbeSets == 0)
			probes_all <- rbind(probes_all,probes)
	}
	
	if (numProbeSets > 0) {
		w <- which(pSets$start==0)
		if (length(w) > 0)
			pSets <- pSets[-w,]
		return(pSets)
	} else {
		w <- which(probes_all$start==0)
		if (length(w) > 0)
			probes_all <- probes_all[-w,]
		return(probes_all)
	}
}
