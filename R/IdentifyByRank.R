IdentifyByRank <- function(dbFile,
	tblName="Seqs",
	level=0,
	add2tbl=FALSE,
	verbose=TRUE) {
	
	# error checking:
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(level) || floor(level) != level)
		stop("level must be an integer.")
	
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
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	
	if (is.na(match("rank",
		dbListFields(dbConn,
			tblName))))
		stop("No rank column in table.")
	
	searchExpression <- paste("select distinct rank from", tblName)
	rs <- dbSendQuery(dbConn, searchExpression)
	x <- dbFetch(rs, n=-1, row.names=FALSE)
	dbClearResult(rs)
	
	.change <- function(id) {
		id <- .Call("replaceChar", id, '"', "", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, "'", "", PACKAGE="DECIPHER")
		id <- gsub("^\\s+|\\s+$", "", id)
		id <- gsub("\\.+$", "", id)
		return(id)
	}
	
	z <- x
	if (level==0) {
		x <- strsplit(x$rank, "\n", fixed=TRUE)
		z$origin <- unlist(lapply(x,
			function (x) {
				x <- paste(x[-1], collapse=" ")
			}))
		z$identifier <- .change(unlist(lapply(x, `[`, 1L)))
	} else {
		x$rank <- unlist(lapply(strsplit(x$rank,
				"\n",
				fixed=TRUE),
			function (x) {
				x <- paste(x[-1], collapse=" ")
			}))
		for (j in seq_along(x$rank)) {
			a <- strsplit(x$rank[j], ";")[[1]]
			l <- length(a)
			if (level < 0) {
				temp_level <- l + level + 1L
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
			z$origin[j] <- unlist(strsplit(as.character(x$rank[j]),
				id,
				fixed=TRUE))[1]
			
			z$identifier[j] <- .change(id)
		}
	}
	
	if (is.character(add2tbl) || add2tbl) {
		dbWriteTable(dbConn, "taxa", z)
		
		if (verbose)
			cat("\nUpdating column: \"identifier\"...")
		
		searchExpression <- paste("update ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			" set identifier = (select identifier from taxa where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank = taxa.rank)",
			sep="")
		dbGetQuery(dbConn, searchExpression)
		
		searchExpression <- "drop table taxa"
		dbGetQuery(dbConn, searchExpression)
	}
	
	if (verbose) {
		cat("\nFormed",
			length(unique(z$identifier)),
			"distinct groups.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to table ",
				ifelse(is.character(add2tbl),add2tbl,tblName),
				": \"identifier\".",
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
