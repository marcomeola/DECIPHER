IdentifyByRank <- function(dbFile,
	tblName="DNA",
	level=3,
	add2tbl=FALSE,
	verbose=TRUE) {
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(level) || level == 0 || floor(level) != level)
		stop("level must be an integer other than zero.")
	
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
			stop("'dbFile' must be a character string or SQLiteConnection.")
		if (!isIdCurrent(dbConn))
			stop("The connection has expired.")
	}
	
	if (is.na(match("rank",
		dbListFields(dbConn,
			tblName))))
		stop("No rank column in table.")
	
	searchExpression <- paste("select rank from", tblName)
	rs <- dbSendQuery(dbConn, searchExpression)
	searchResult <- fetch(rs, n=-1)
	dbClearResult(rs)
	f <- factor(searchResult$rank)
	x <- data.frame(table(f))
	
	z <- x
	for (j in 1:length(x$f)) {
		a <- strsplit(as.character(x$f[j]), ";")[[1]]
		l <- length(a)
		if (level < 0) {
			temp_level <- l + level
		} else {
			temp_level <- level
		}
		
		if (temp_level > l) {
			id <- as.character(a[l])
		} else if (temp_level < 1) {
			id <- as.character(a[1])
		} else {
			id <- as.character(a[temp_level])
		}
		z$origin[j] <- unlist(strsplit(as.character(x$f[j]),
			id,
			fixed=TRUE))[1]
		id <- .Call("replaceChar", id, ";", "", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, ".", "", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, '"', "", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, "'", "", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, " ", "", PACKAGE="DECIPHER")
		z$id[j] <- as.character(id)
	}
	
	if (is.character(add2tbl) || add2tbl) {
		dbWriteTable(dbConn, "taxa", z)
		
		if (verbose)
			cat("\nUpdating column: \"id\"...")
		
		searchExpression <- paste("update or replace ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			" set id = (select id from taxa where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank = taxa.f) where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank in (select f from taxa)",
			sep="")
		dbGetQuery(dbConn, searchExpression)
		
		searchExpression <- "drop table taxa"
		dbGetQuery(dbConn, searchExpression)
	}
	
	if (verbose) {
		cat("\nFormed",
			length(unique(z$id)),
			"distinct groups.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to table ",
				ifelse(is.character(add2tbl),add2tbl,tblName),
				": \"id\".",
				sep="")
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	names(z)[1] <- "rank"
	return(z)
}
