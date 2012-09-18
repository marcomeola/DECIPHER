FormGroups <- function(dbFile,
	tblName="DNA",
	goalSize=1000,
	minGroupSize=500,
	maxGroupSize=10000,
	add2tbl=FALSE,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(goalSize) || goalSize <= 0)
		stop("goalSize must be a numeric greater than zero.")
	if (!is.numeric(minGroupSize) || minGroupSize < 0)
		stop("minGroupSize must be a numeric greater than or equal to zero.")
	if (!is.numeric(maxGroupSize) || minGroupSize > maxGroupSize)
		stop("maxGroupSize must be a numeric greater than or equal to minGroupSize.")
	if (verbose)
		time.1 <- Sys.time()
	
	# initialize databases
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
	
	if (is.na(match("rank",
		dbListFields(dbConn,
			tblName))))
		stop("No rank column in table.")
	
	searchExpression <- paste("select distinct rank, count(rank) as count from",
		tblName,
		"group by rank")
	rs <- dbSendQuery(dbConn, searchExpression)
	searchResult <- fetch(rs, n=-1)
	dbClearResult(rs)
	
	searchResult <- searchResult[order(searchResult$count,
		decreasing=FALSE),]
	
	searchResult$id <- ""
	searchResult$oneup <- ""
	searchResult$origin <- ""
	
	if (verbose)
		pBar <- txtProgressBar(style=3)
	
	for (i in 1:length(searchResult$rank)) {
					# if (length(w1) > 0)
						# warning(lineage[j],
							# " replaced ",
							# paste(searchResult$id[w[w1]],
								# collapse=", "),
							# ".",
							#sep="")
							if ((counts + sum(-1*searchResult$count[w2[w3]])) < maxGroupSize) {
								# total group size will be less than than the maximum
							}
					if (substr(origin,
						nchar(origin),
						nchar(origin))=='"')
						origin <- substr(origin,
							1,
							nchar(origin) - 1)
					# w1 <- which(searchResult$id[w]!="" && searchResult$count[w] > 0)
					# if (length(w1) > 0)
						# warning(lineage[j],
							# " replaced ",
							# paste(searchResult$id[w[w1]],
								# collapse=", "),
							# ".",
							# sep="")
							if ((counts + sum(-1*searchResult$count[w2[w3]])) < maxGroupSize) {
								# total group size will be less than than the maximum
								searchResult$count[w2[w3]] <- -1*searchResult$count[w2[w3]]
								w <- c(w, w2[w3])
							}
					if (substr(origin,
						nchar(origin),
						nchar(origin))=='"')
						origin <- substr(origin,
							1,
							nchar(origin) - 1)
		if (verbose)
			setTxtProgressBar(pBar, i/length(searchResult$rank))
	
	if (is.character(add2tbl) || add2tbl) {
		if (verbose)
			cat("\nUpdating column: \"id\"...")
			tblName,
			" set id = (select id from taxa where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank = taxa.rank) where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank in (select rank from taxa)",
			sep="")
		if (is.na(match("origin",
			dbListFields(dbConn,
				ifelse(is.character(add2tbl),add2tbl,tblName))))) {
				ifelse(is.character(add2tbl),add2tbl,tblName),
				" add column origin",
				sep="")
		
		if (verbose)
			cat("\nUpdating column: \"origin\"...")
			ifelse(is.character(add2tbl),add2tbl,tblName),
			" set origin = (select origin from taxa where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank = taxa.rank) where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank in (select rank from taxa)",
			sep="")
	
	
	if (verbose) {
		cat("\nFormed",
			length(unique(searchResult$id)),
			"distinct groups.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to table ",
				ifelse(is.character(add2tbl),add2tbl,tblName),
				": \"id\", \"origin\".",
				sep="")
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	searchResult <- searchResult[,-match("oneup",
		names(searchResult))]
	
	return(searchResult)
}