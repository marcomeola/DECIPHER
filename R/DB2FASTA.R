DB2FASTA <- function(file,
	dbFile,
	tblName="DNA",
	identifier="",
	limit=-1,
	replaceChar=NULL,
	orderBy="row_names",
	append=FALSE,
	comments=FALSE,
	removeGaps="none",
	chunkSize=1e5,
	verbose=TRUE,
	...) {
	
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
	if (!is.numeric(chunkSize))
		stop("chunkSize must be a numeric.")
	if (floor(chunkSize)!=chunkSize)
		stop("chunkSize must be a whole number.")
	if (chunkSize <= 0)
		stop("chunkSize must be greater than zero.")
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
			stop("'dbFile' must be a character string or SQLiteConnection.")
		if (!isIdCurrent(dbConn))
			stop("The connection has expired.")
	}
	
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(style=3)
	}
	
	con <- file(file, ifelse(append, "a", "w"))
	
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
	if (limit > 0) {
		searchExpression1 <- paste(searchExpression,
			'limit',
			limit)
	} else {
		searchExpression1 <- searchExpression
	}
	searchExpression1 <- paste('select count(*) from ',
		searchExpression1,
		sep="")
	rs <- dbSendQuery(dbConn, searchExpression1)
	count <- as.numeric(fetch(rs, n=-1))
	dbClearResult(rs)
	
	if (orderBy!="row_names") # default ordering is row_names
		searchExpression <- paste(searchExpression,
			'order by',
			orderBy)
	
	s <- seq(1, count, chunkSize)
	for (i in 1:length(s)) {
		# build the search expression
		searchExpression1 <- paste(searchExpression,
			'limit',
			ifelse(i==length(s), count - s[i] + 1, chunkSize),
			'offset',
			ifelse(i==1, 0, s[i] - 1))
		searchExpression1 <- paste(ifelse(comments,
				'select * from ',
				'select row_names, description from '),
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
		
		m <- match(searchResult$row_names,
			searchResult2$row_names)
		searchResult2 <- searchResult2[m,]
		
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
		
		fileText <- character(2*dim(searchResult2)[1])
		fileText[seq(2, length(fileText), 2)] <- searchResult2$sequence
		if (comments) {
			fileText[seq(1, length(fileText), 2)] <- do.call("paste", c(searchResult, sep = "; "))
		} else {
			fileText[seq(1, length(fileText), 2)] <- searchResult$description
		}
		
		writeLines(fileText, con)
		
		if (verbose)
			setTxtProgressBar(pBar, ifelse(i==length(s), 1, (s[i + 1] - 1)/count))
	}
	close(con)
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\nWrote ",
			count,
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