Array2Matrix <- function(probes,
	verbose=TRUE) {
	
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (verbose) {
		pBar <- txtProgressBar(style=3)
		time.1 <- Sys.time()
	}
	
	MMs <- strsplit(probes$mismatch, ", ", fixed=TRUE)
	
	hyb_effs <- list()
	l <- length(MMs)
	for (i in 1:l) {
		if (length(MMs[[i]])==0)
			next
		
		MMs[[i]] <- unlist(strsplit(MMs[[i]], " ", fixed=TRUE))
		
		num <- length(MMs[[i]])
		
		hyb_effs[[i]] <- MMs[[i]][2*(1:(num/2))]
		hyb_effs[[i]] <- .Call("replaceChar",
			hyb_effs[[i]],
			"%",
			"",
			PACKAGE="DECIPHER")
		hyb_effs[[i]] <- .Call("replaceChar",
			hyb_effs[[i]],
			"(",
			"",
			PACKAGE="DECIPHER")
		hyb_effs[[i]] <- as.numeric(.Call("replaceChar",
			hyb_effs[[i]],
			")",
			"",
			PACKAGE="DECIPHER"))
		
		MMs[[i]] <- MMs[[i]][2*(1:(num/2)) - 1]
		
		if (verbose)
			setTxtProgressBar(pBar, i/l)
	}
	
	index <- lapply(1:l, function(l) {
		rep(l, length(MMs[[l]]))
	})
	
	p_names <- unique(probes$name)
	MMs <- match(unlist(MMs), p_names)
	
	A <- list(i=c(1:l,unlist(index)),
		j=c(match(probes$name, p_names),MMs),
		x=c(as.numeric(probes$hyb_eff),unlist(hyb_effs))/100,
		dimnames=list(1:l, p_names))
	
	if (verbose) {
		close(pBar)
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(A)
}
