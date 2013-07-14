NNLS <- function(A,
	b,
	precision=sqrt(.Machine$double.eps),
	verbose=TRUE) {
	if (length(A) != 4)
		stop("A must have four components: A$i, A$j, A$x, and A$dimnames.")
	if (!is.integer(A$i))
		stop("Rows (i) must be an integer vector.")
	if (!is.integer(A$j))
		stop("Columns (j) must be an integer vector.")
	if (!is.numeric(A$x))
		stop("Values (x) must be a numeric vector.")
	if (length(A$i) != length(A$j))
		stop("The length of columns (j) and rows (i) must be equal.")
	if (length(A$i) != length(A$x))
		stop("The length of rows (i) and values (x) must be equal.")
	if (max(A$j) > length(A$dimnames[[2]]))
		stop("More columns than column names.")
	if (max(A$i) > length(A$dimnames[[1]]))
		stop("More rows than row names.")
	if (!is.numeric(b))
		stop("b must be a numeric vector or matrix.")
	if (!((length(b) %% length(A$dimnames[[1]]))==0))
		stop("The length of b must be a multiple of the number of rows in A.")
	if (class(b)=="matrix")
		if (nrow(b)!=length(A$dimnames[[1]]))
			stop("The number of rows in b must equal the number of rows in A.")
	if (!is.numeric(precision))
		stop("precision must be a numeric.")
	if (precision <= 0)
		stop("precision must be a positive number.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(max=100, style=3)
	} else {
		pBar <- NULL
	}
	
	o <- order(A$i)
	x <- .Call("NNLS",
		A$i[o],
		A$j[o],
		A$x[o],
		length(A$dimnames[[1]]),
		length(A$dimnames[[2]]),
		b,
		precision,
		verbose,
		pBar,
		PACKAGE="DECIPHER")
	b <- matrix(b, ncol=ncol(x))
	res <- b - .Call("sparseMult",
		A$i,
		A$j,
		A$x,
		length(A$dimnames[[1]]),
		length(A$dimnames[[2]]),
		x,
		PACKAGE="DECIPHER")
	
	rownames(x) <- A$dimnames[[2]]
	rownames(res) <- A$dimnames[[1]]
	
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
	
	return(list(x=x, residuals=res))
}
