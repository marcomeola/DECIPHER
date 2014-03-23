Add2DB <- function(myData,
	dbFile,
	tblName="DNA",
	clause="",
	verbose=TRUE) {
	
	time.1 <- Sys.time()
	
	# error checking
	if (!is.data.frame(myData))
		stop("myData must be a data frame object.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (substr(tblName, 1, 1) == "_")
		stop("Invalid tblName.")
	if (tblName == "taxa")
		stop("taxa is a reserved tblName.")
	if (!is.character(clause))
		stop("clause must be a character string.")
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
	
	result  <- dbListTables(dbConn)
	w <- which(result==tblName)
	if (length(w)==0) { # need to make the table
		dbWriteTable(dbConn, tblName, myData, row.names=TRUE,
			overwrite=TRUE, append=FALSE)
		if (verbose)
			cat("Created table", tblName, "\n")
	} else {
		result <- dbListFields(dbConn, tblName)
		colIDs <- names(myData)
		
		if (is.na(match("row_names", colIDs)))
			myData$row_names=row.names(myData)
		
		# loop through each of the columns to add
		for (i in 1:length(colIDs)) {
			colName <- colIDs[i]
			
			if (is.na(match(colName, result))) {
				# first add the column if it does not already exist
				expression1 <- paste("alter table ",
					tblName,
					" add column ",
					colName,
					" ",
					sqliteDataType(myData[colIDs[i],1]),
					sep="")
				if (verbose)
					cat("Expression:  ", expression1, "\n", sep="")
				dbGetQuery(dbConn, expression1)
			}
			
			# next update the column with new data
			expression2 <- paste("update or replace ",
				tblName,
				" set ",
				colName,
				" = :",
				colName,
				" where row_names = :row_names",
				sep="")
			if (clause!="")
				expression2 <- paste(expression2,
					clause)
			if (verbose)
				cat("Expression: ", expression2, "\n\n")
			dbBeginTransaction(dbConn)
			dbGetPreparedQuery(dbConn,
				expression2,
				bind.data=myData)
			if (!dbCommit(dbConn)) {
				warning("Unsucessful transaction!")
				invisible(FALSE)
			}
		}
	}
	
	if (verbose) { # print the elapsed time to update table
		cat("Added to table ",
			tblName,
			":  \"",
			paste(names(myData)[-match("row_names",
				names(myData))],
			collapse="\" and \""),
			"\".",
			sep="")
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	invisible(TRUE)
}