DB2FASTA <- function(file,
	dbFile,
	tblName="DNA",
	identifier="",
	limit=-1,
	replaceChar=NULL,
	orderBy="row_names",
	append=FALSE,
	comments=TRUE,
	removeGaps="none",
	verbose=TRUE,
	...) {
	
	time.1 <- Sys.time()
	
	# error checking
	if (!is.character(file))
		if (!inherits(file, "connection"))
			stop("file must be a character string or connection.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.numeric(limit))
		stop("limit must be a numeric.")
	if (floor(limit)!=limit)
		stop("limit must be a whole number.")
	if (!is.character(orderBy))
		stop("orderBy must be a character string.")
	if (!is.logical(append))
		stop("append must be a logical.")
	if (!is.logical(comments))
		stop("comments must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps, GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
	
	# initialize database
	driver = dbDriver("SQLite")
	if (is.character(dbFile)) {
		dbConn = dbConnect(driver, dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn = dbFile
		if (!inherits(dbConn,"SQLiteConnection")) 
			stop("'dbFile' must be a character string or SQLiteConnection")
		if (!isIdCurrent(dbConn))
			stop("The connection has expired.")
	}
	
	# build the search expression
	searchExpression <- tblName
	args <- list(...)
	if (identifier!="" ||
		length(args) > 0)
		searchExpression <- paste(searchExpression,
			' where',
			sep="")
	if (identifier!="")
		searchExpression <- paste(searchExpression,
			' id like "',
			identifier,
			'"',
			sep="")
	firstTime <- TRUE
	for (a in args) {
		if (identifier!="" ||
			!firstTime)
			searchExpression <- paste(searchExpression,
				'and')
		searchExpression <- paste(searchExpression,
				a)
		firstTime <- FALSE
	}
	if (orderBy!="row_names") # default ordering is row_names
		searchExpression <- paste(searchExpression,
			'order by',
			orderBy)
	if (limit > 0)
		searchExpression <- paste(searchExpression,
			'limit',
			limit)
	
	searchExpression1 <- paste('select * from ',
		searchExpression,
		sep="")
	
	if (verbose)
		cat("Search Expression:",
			"\n",
			searchExpression1,
			sep="")
	
	rs <- dbSendQuery(dbConn, searchExpression1)
	searchResult <- fetch(rs, n=-1)
	dbClearResult(rs)
	
	searchExpression2 <- paste('select row_names, sequence from _',
		tblName,
		" where row_names in (select row_names from ",
		searchExpression,
		")",
		sep="")
	rs <- dbSendQuery(dbConn, searchExpression2)
	searchResult2 <- fetch(rs, n=-1)
	dbClearResult(rs)
	
	# decompress the resulting sequences
	searchResult2$sequence <- lapply(searchResult2$sequence,
		memDecompress,
		type="gzip",
		asChar=TRUE)
	searchResult2$sequence <- paste(searchResult2$sequence)
	
	if (!is.null(replaceChar))
		# replace all characters not in the DNA_ALPHABET
		searchResult2$sequence <- .Call("replaceChars",
			searchResult2$sequence,
			replaceChar,
			PACKAGE="DECIPHER")
	
	# remove gaps if applicable
	if (removeGaps==2) {
		searchResult2$sequence <- .Call("replaceChar",
			searchResult2$sequence,
			"-",
			"",
			PACKAGE="DECIPHER")
	} else if (removeGaps==3) {
		searchResult2$sequence <- .Call("commonGaps",
			searchResult2$sequence,
			PACKAGE="DECIPHER")
	}
	
	myDNAStringSet <- DNAStringSet(searchResult2$sequence)
	if (comments) {
		names(myDNAStringSet) <- do.call("paste", c(searchResult, sep = "; "))
	} else {
		names(myDNAStringSet) <- searchResult$description
	}
	
	if (verbose)
		cat("\nWriting FASTA file...")
	writeXStringSet(myDNAStringSet,
		file=file,
		append=append)
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\nWrote ",
			length(myDNAStringSet),
			" sequences.",
			"\n",
			sep="")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	return(TRUE)
}