FormGroups <- function(dbFile,
	tblName="Seqs",
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
			stop("'dbFile' must be a character string or SQLiteConnection.")
		if (!dbIsValid(dbConn))
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
	
	searchResult$origin <- ""
	searchResult$identifier <- ""
	
	if (verbose)
		pBar <- txtProgressBar(style=3)
	
	rank <- unlist(lapply(strsplit(searchResult$rank,
			"\n",
			fixed=TRUE),
		function (x) {
			x <- paste(x[-1], collapse=" ")
		}))
	for (i in seq_along(rank)) {		if (searchResult$identifier[i]=="") {			lineage <- unlist(strsplit(as.character(rank[i]),				";",				fixed=TRUE))
			for (j in length(lineage):1) {				w <- which(grepl(paste(lineage[1:j], collapse=";"),
					rank,
					fixed=TRUE))				counts <- sum(abs(searchResult$count[w]))				
				if (counts > goalSize) {					if (counts > maxGroupSize && j < length(lineage)) {						j <- j + 1 # go down one rank						w <- which(grepl(paste(lineage[1:j], collapse=";"),
							rank,
							fixed=TRUE))
						counts <- sum(abs(searchResult$count[w]))					}
					
					if (j > 1) {
						origin <- paste(lineage[1:j], collapse=";")
					} else {
						origin <- ""
					}
										if (counts < minGroupSize) { # mark for later inclusion						searchResult$count[w] <- -abs(searchResult$count[w])					} else if (j > 1) { # try to include a little more						w2 <- which(searchResult$count < 0)						if (length(w2) > 0) {							w3 <- which(origin==searchResult$origin[w2])							if (length(w3) > 0) {
								if ((counts - sum(searchResult$count[w2[w3]])) <= maxGroupSize) {
									searchResult$count[w2[w3]] <- -abs(searchResult$count[w2[w3]])									w <- c(w, w2[w3])
								}							}
						}					}
										searchResult$origin[w] <- origin										# remove some punctuation					lineage[j] <- gsub('"', '', lineage[j], fixed=TRUE)					lineage[j] <- gsub('.', '', lineage[j], fixed=TRUE)					searchResult$identifier[w] <- gsub("^\\s+|\\s+$", "", lineage[j])
					break				} else if (j==1 && counts < goalSize) {
					searchResult$origin[i] <- ""
					
					# remove some punctuation
					lineage[j] <- gsub('"', '', lineage[j], fixed=TRUE)
					lineage[j] <- gsub('.', '', lineage[j], fixed=TRUE)
					searchResult$identifier[i] <- gsub("^\\s+|\\s+$", "", lineage[j])
				}
			}		}
		if (verbose)
			setTxtProgressBar(pBar, i/length(rank))	}
	
	if (is.character(add2tbl) || add2tbl) {		dbWriteTable(dbConn, "taxa", searchResult)		
		if (verbose)
			cat("\nUpdating column: \"identifier\"...")		searchExpression <- paste("update or replace ",
			tblName,
			" set identifier = (select identifier from taxa where ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			".rank = taxa.rank) where ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			".rank in (select rank from taxa)",
			sep="")		dbGetQuery(dbConn, searchExpression)		
		if (is.na(match("origin",
			dbListFields(dbConn,
				ifelse(is.character(add2tbl), add2tbl, tblName))))) {			searchExpression <- paste("alter table ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				" add column origin",
				sep="")			dbGetQuery(dbConn, searchExpression)		}
		
		if (verbose)
			cat("\nUpdating column: \"origin\"...")		searchExpression <- paste("update or replace ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			" set origin = (select origin from taxa where ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			".rank = taxa.rank) where ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			".rank in (select rank from taxa)",
			sep="")		dbGetQuery(dbConn, searchExpression)				searchExpression <- "drop table taxa"		dbGetQuery(dbConn, searchExpression)	}
	
	if (verbose) {
		cat("\nFormed",
			length(unique(searchResult$identifier)),
			"distinct groups.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to table ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				": \"identifier\", \"origin\".",
				sep="")
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	searchResult$count <- abs(searchResult$count)
	
	return(searchResult)
}
