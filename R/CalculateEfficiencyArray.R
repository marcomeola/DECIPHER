CalculateEfficiencyArray <- function(probe,
	target,
	FA=0,
	dGini=1.96,
	Po=10^-2.0021,
	m=0.1731,
	temp=42,
	deltaGrules=NULL) {
	
	# error checking
	if (is(probe, "DNAStringSet"))
		probe <- strsplit(toString(probe), ", ", fixed=TRUE)[[1]]
	if (is(target, "DNAStringSet"))
		target <- strsplit(toString(target), ", ", fixed=TRUE)[[1]]
	if (!is.character(probe))
		stop("probe must be a character vector.")
	if (!is.character(target))
		stop("target must be a character vector.")
	if (any(nchar(target) != nchar(probe)))
		stop("probe and target must be aligned (equal length).")
	if (!is.numeric(Po))
		stop("Po must be a numeric.")
	if (!(Po > 0))
		stop("Po must be greater than zero.")
	if (!is.numeric(m))
		stop("m must be a numeric.")
	if (!(m > 0))
		stop("m must be greater than zero.")
	if (!is.numeric(dGini))
		stop("dGini must be a numeric.")
	if (!is.numeric(FA))
		stop("FA must be a numeric.")
	if (any(FA < 0))
		stop("FA must be greater than or equal to zero.")
	if (!is.numeric(temp))
		stop("temp must be a numeric.")
	if (temp < -273)
		stop("temp must be greater than or equal to absolute zero.")
	
	if (is.null(deltaGrules)) {
		data(deltaGrules, envir=environment())
	} else {
		if (!is.numeric(deltaGrules))
			stop("deltaGrules must be numeric.")
		if (length(deltaGrules)!=390625)
			stop("deltaGrules must be of dimensions 5 x 5 x 5 x 5 x 5 x 5 x 5 x 5.")
	}
	
	l <- length(probe)
	if (l==0)
		stop("No probe specified.")
	if (l!=length(target))
		stop("probe is not the same length as target.")
	
	dG <- dGini + .Call("calculateDeltaG", probe, target, deltaGrules, PACKAGE="DECIPHER")
	
	RT <- .0019871*(273.15 + temp) # [kcal/mol]
	deltaG <- matrix(0, nrow=l, ncol=length(FA), dimnames=list(1:l, paste("dG", FA, sep="_")))
	eff <- matrix(0, nrow=l, ncol=length(FA), dimnames=list(1:l, paste("HybEff", FA, sep="_")))
	for (i in 1:length(FA)) {
		deltaG[, i] <- dG + m*FA[i]
		eff[, i] <- Po*exp(-deltaG[, i]/RT)/(1 + Po*exp(-deltaG[, i]/RT))
	}
	
	return(cbind(eff, deltaG))
}
