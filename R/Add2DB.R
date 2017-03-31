# below function from RSQLite package version 0.11.4
.sqliteDataType <- function(obj, ...) {
	rs.class <- data.class(obj)
	rs.mode <- storage.mode(obj)
	switch(rs.class,
		numeric = if (rs.mode=="integer") "INTEGER" else "REAL",
		character = "TEXT",
		logical = "INTEGER",
		factor = "TEXT",
		ordered = "TEXT",
		## list maps to BLOB. Although not checked, the list must
		## either be empty or contain only raw vectors or NULLs.
		list = "BLOB",
		## attempt to store obj according to its storage mode if it has
		## an unrecognized class.
		switch(rs.mode,
			integer = "INTEGER",
			double = "REAL",
			## you'll get this if class is AsIs for a list column
			## within a data.frame
			list = if (rs.class == "AsIs") "BLOB" else "TEXT",
			"TEXT"))
}

Add2DB <- function(myData,
	dbFile,
	tblName="Seqs",
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
	if (ncol(myData)==0)
		stop("myData contains no columns.")
	if (any(grepl(".", names(myData), fixed=TRUE)))
		stop("Column names cannot contain periods.")
	
	# initialize database
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
	
	result  <- dbListTables(dbConn)
	w <- which(result==tblName)
	if (length(w)==0) { # need to make the table
		dbWriteTable(dbConn, tblName, myData, row.names=TRUE,
			overwrite=TRUE, append=FALSE)
		if (verbose)
			cat("Created table '", tblName, "'.\n", sep="")
	} else {
		result <- dbListFields(dbConn, tblName)
		colIDs <- names(myData)
		
		if (is.na(match("row_names", colIDs)))
			myData$row_names=row.names(myData)
		
		x <- sqliteQuickColumn(dbConn, tblName, "row_names")
		m <- match(myData$row_names, x)
		if (any(is.na(m)))
			stop("row.names of myData are missing from '", tblName, "'.")
		
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
					.sqliteDataType(myData[colIDs[i], 1]),
					sep="")
				if (verbose)
					cat("Expression:\n",
						paste(strwrap(expression1,
								width=getOption("width") - 1L),
							collapse="\n"),
						"\n\n",
						sep="")
				dbGetQuery(dbConn, expression1)
			}
			
			# next update the column with new data
			expression2 <- paste("update ",
				tblName,
				" set ",
				colName,
				" = :",
				colName,
				" where row_names = :row_names",
				sep="")
			if (clause!="")
				expression2 <- paste(expression2,
					" where ",
					clause,
					sep="")
			if (verbose)
				cat("Expression:\n",
					paste(strwrap(expression2,
							width=getOption("width") - 1L),
						collapse="\n"),
					"\n\n",
					sep="")
			rs <- dbSendQuery(dbConn, expression2)
			dbBind(rs, myData[, c("row_names", colIDs[i])])
			dbFetch(rs, n=-1, row.names=FALSE)
			dbClearResult(rs)
		}
		if (verbose) {
			cat("Added to table ",
				tblName,
				":  \"",
				paste(names(myData)[-match("row_names",
					names(myData))],
				collapse="\" and \""),
				"\".\n",
				sep="")
		}
	}
	
	if (verbose) { # print the elapsed time to update table
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	invisible(TRUE)
}
