IdConsensus <- function(dbFile,
	tblName="DNA",
	identifier="",
	colName="cluster",
	add2tbl=FALSE,
	verbose=TRUE,
	...) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(colName))
		stop("colName must be a character string.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	
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
	
	f <- dbListFields(dbConn, tblName)
	w <- which(f==colName)
	if (length(w)==0)
		stop(paste("The colName '", colName, "' does not exist.", sep=""))
	
	searchExpression <- paste('select distinct ',
		colName,
		' from ',
		tblName,
		sep="")
	if (identifier != "")
		searchExpression <- paste(searchExpression,
			" where id like '",
			identifier,
			"'",
			sep="")
	rs <- dbSendQuery(dbConn, searchExpression)
	groups <- fetch(rs, n=-1)[,eval(colName)]
	dbClearResult(rs)
	
	# remove any null groups in database
	w <- which(is.na(groups))
	if (length(w) > 0)
		groups <- groups[-w]
	
	if (length(groups) < 1)
		stop("No groups in which to form consensus.")
	
	# initialize a progress bar
	if (verbose)
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
	
	consensus <- DNAStringSet()
	seqCount <- numeric(length(groups))
	j <- 0L
	for (i in groups) {
		j <- j + 1L
		dna_subset <- SearchDB(dbFile,
			tblName=tblName,
			verbose=FALSE,
			identifier=identifier,
			...=paste(colName,
				"= '",
				gsub("'", "''", i, fixed=TRUE),
				"'",
				sep=""))
		
		if (length(consensus)==0)
			consensus <- ConsensusSequence(dna_subset,
					verbose=FALSE,
					...)
		else
			consensus <-c(consensus,
				 ConsensusSequence(dna_subset,
					verbose=FALSE,
					...))
		
		seqCount[j] <- length(dna_subset)
		if (verbose)
			setTxtProgressBar(pBar, 100*j/length(groups))
	}
	if (identifier=="") {
		names(consensus) <- paste(groups,
			"_",
			seqCount,
			"seqs",
			sep="")
	} else {
		names(consensus) <- paste(groups,
			"_",
			seqCount,
			"seqs",
			"_",
			identifier,
			sep="")
	}
	
	if (is.character(add2tbl) || add2tbl) {
		Seqs2DB(consensus,
			type="D",
			dbFile=dbFile,
			tblName=ifelse(is.character(add2tbl),add2tbl,tblName),
			identifier=identifier,
			verbose=FALSE)
	}
	
	if (verbose) {
		close(pBar)
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded ",
				length(groups),
				" consensus sequences to ",
				ifelse(is.character(add2tbl),add2tbl,tblName),
				".",
				sep="")
		
		cat("\nFound consensus for ",
			length(groups),
			" groups.",
			sep="")
		
		cat("\n")
		time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
		cat("\n")
	}
	return(consensus)
}