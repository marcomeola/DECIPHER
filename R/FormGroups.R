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
	
	for (i in 1:length(searchResult$rank)) {		if (searchResult$id[i]=="") {			lineage <- unlist(strsplit(as.character(searchResult$rank[i]),				"; ",				fixed=TRUE))			for (j in length(lineage):1) {				if (substr(lineage[j], 1, 1)==" ") {					# remove leading spaces					lineage[j] <- substr(lineage[j],						2,						nchar(lineage[j]))				}								w <- which(grepl(lineage[j], searchResult$rank, fixed=TRUE))				counts <- sum(searchResult$count[w])				if (counts >= maxGroupSize) {					if (j < length(lineage))						j <- j + 1 # go down one rank					w <- which(grepl(lineage[j], searchResult$rank, fixed=TRUE))										# report conflicts					# w1 <- which(searchResult$id[w]!="" && searchResult$count[w] > 0)
					# if (length(w1) > 0)
						# warning(lineage[j],
							# " replaced ",
							# paste(searchResult$id[w[w1]],
								# collapse=", "),
							# ".",
							#sep="")										# mark for later inclusion					counts <- sum(searchResult$count[w])					if (counts < minGroupSize) {						searchResult$count[w] <- -counts						if (j > 3)							searchResult$oneup[w] <- lineage[j - 1]					} else if (j > 3) {						# try to include a little more						w2 <- which(searchResult$count < 0)						w3 <- c()						if (length(w2) > 0)							w3 <- which(grepl(lineage[j - 1], searchResult$oneup[w2], fixed=TRUE))						if (length(w3) > 0) {							searchResult$count[w2[w3]] <- -1*searchResult$count[w2[w3]]							w <- c(w, w2[w3])						}					}										# record origin before lineage					origin <- unlist(strsplit(as.character(searchResult$rank[i]),						lineage[j],						fixed=TRUE))[1]
					if (substr(origin,
						nchar(origin),
						nchar(origin))=='"')
						origin <- substr(origin,
							1,
							nchar(origin) - 1)					searchResult$origin[w] <- origin										# remove some punctuation					lineage[j] <- gsub('"', '', lineage[j], fixed=TRUE)					lineage[j] <- gsub('.', '', lineage[j], fixed=TRUE)					searchResult$id[w] <- lineage[j]					break				} else if (counts > goalSize && counts < maxGroupSize) {					# report conflicts
					# w1 <- which(searchResult$id[w]!="" && searchResult$count[w] > 0)
					# if (length(w1) > 0)
						# warning(lineage[j],
							# " replaced ",
							# paste(searchResult$id[w[w1]],
								# collapse=", "),
							# ".",
							# sep="")										# try to include a little more					if (j > 3) {						w2 <- which(searchResult$count < 0)						w3 <- c()						if (length(w2) > 0)							w3 <- which(grepl(lineage[j - 1], searchResult$oneup[w2], fixed=TRUE))						if (length(w3) > 0) {							searchResult$count[w2[w3]] <- -1*searchResult$count[w2[w3]]							w <- c(w, w2[w3])						}					}										# record origin before lineage					origin <- unlist(strsplit(as.character(searchResult$rank[i]),						lineage[j],						fixed=TRUE))[1]
					if (substr(origin,
						nchar(origin),
						nchar(origin))=='"')
						origin <- substr(origin,
							1,
							nchar(origin) - 1)					searchResult$origin[w] <- origin										# remove some punctuation					lineage[j] <- gsub('"', '', lineage[j], fixed=TRUE)					lineage[j] <- gsub('.', '', lineage[j], fixed=TRUE)					searchResult$id[w] <- lineage[j]					break				}			}		}
		if (verbose)
			setTxtProgressBar(pBar, i/length(searchResult$rank))	}
	
	if (is.character(add2tbl) || add2tbl) {		dbWriteTable(dbConn, "taxa", searchResult)		
		if (verbose)
			cat("\nUpdating column: \"id\"...")		searchExpression <- paste("update or replace ",
			tblName,
			" set id = (select id from taxa where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank = taxa.rank) where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank in (select rank from taxa)",
			sep="")		dbGetQuery(dbConn, searchExpression)		
		if (is.na(match("origin",
			dbListFields(dbConn,
				ifelse(is.character(add2tbl),add2tbl,tblName))))) {			searchExpression <- paste("alter table ",
				ifelse(is.character(add2tbl),add2tbl,tblName),
				" add column origin",
				sep="")			dbGetQuery(dbConn, searchExpression)		}
		
		if (verbose)
			cat("\nUpdating column: \"origin\"...")		searchExpression <- paste("update or replace ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			" set origin = (select origin from taxa where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank = taxa.rank) where ",
			ifelse(is.character(add2tbl),add2tbl,tblName),
			".rank in (select rank from taxa)",
			sep="")		dbGetQuery(dbConn, searchExpression)				searchExpression <- "drop table taxa"		dbGetQuery(dbConn, searchExpression)	}
	
	
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
			digits=1))
		cat("\n")
	}
	
	searchResult <- searchResult[,-match("oneup",
		names(searchResult))]
	
	return(searchResult)
}
