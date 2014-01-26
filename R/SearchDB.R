SearchDB <- function(dbFile,
	tblName="DNA",
	identifier="",
	type="DNAStringSet",
	limit=-1,
	replaceChar="-",
	nameBy="row_names",
	orderBy="row_names",
	countOnly=FALSE,
	removeGaps="none",
	verbose=TRUE,
	...) {
	
	time.1 <- Sys.time()
	
	# error checking
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps[1], GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet",
		"QualityScaledDNAStringSet", "QualityScaledRNAStringSet", "QualityScaledAAStringSet", "QualityScaledBStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(orderBy))
		stop("orderBy must be a character string.")
	if (!is.logical(countOnly))
		stop("countOnly must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (type==1 || type==5) {
		if (is.na(pmatch(replaceChar, DNA_ALPHABET)) && (replaceChar!=""))
			stop("replaceChar must be a character in the DNA_ALPHABET or empty character.")
	} else if (type==2 || type==6) {
		if (is.na(pmatch(replaceChar, RNA_ALPHABET)) && (replaceChar!=""))
			stop("replaceChar must be a character in the RNA_ALPHABET or empty character.")
	} else if (type==3 || type==7) {
		if (is.na(pmatch(replaceChar, AA_ALPHABET)) && (replaceChar!=""))
			stop("replaceChar must be a character in the AA_ALPHABET or empty character.")
	}
	if (type > 4 && removeGaps > 1)
		stop(paste('removeGaps must be "none" when type is ', TYPES[type], '.', sep=''))
	if (is.numeric(limit)) {
		if (floor(limit)!=limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
	} else {
		if (!grepl("[0-9],[0-9]", limit, perl=T)) {
			limit <- as.numeric(limit)
			if (floor(limit)!=limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
		}
	}
	if (limit > 0 && countOnly)
		stop("limit cannot be specified when countOnly is TRUE.")
	
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
	
	# build the search expression
	if (countOnly) {
		searchExpression <- paste('select count(*) from ',
			tblName,
			sep="")
	} else if (nameBy=="row_names" && orderBy=="row_names") {
		searchExpression <- paste('select row_names, sequence',
			ifelse(type > 4, ', quality', ''),
			' from _',
			tblName,
			" where row_names in (select row_names from ",
			tblName,
			sep="")
	} else {
		searchExpression <- paste('select ',
			ifelse(nameBy=="row_names",
				paste(tblName, ".row_names", sep=""),
				nameBy),
			', _',
			tblName,
			'.sequence',
			ifelse(type > 4,
				paste(', _',
					tblName,
					'.quality',
					sep=""),
				''),
			' from ',
			tblName,
			' join _',
			tblName,
			' on ',
			tblName,
			'.row_names',
			' = _',
			tblName,
			'.row_names where _',
			tblName,
			'.row_names in (select row_names from ',
			tblName,
			sep="")
	}
	
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
	if (!countOnly)
		searchExpression <- paste(searchExpression, ")", sep="")
	
	if (orderBy!="row_names") # default ordering is row_names
		searchExpression <- paste(searchExpression,
			'order by',
			orderBy)
	if (limit > 0)
		searchExpression <- paste(searchExpression,
			'limit',
			limit)
	
	if (verbose)
		cat("Search Expression:","\n",searchExpression,sep="")
	
	rs <- dbSendQuery(dbConn, searchExpression)
	searchResult <- fetch(rs, n=-1)
	dbClearResult(rs)
	
	if (countOnly) {
		count <- as.integer(searchResult)
	} else {
		# decompress the resulting sequences
		searchResult$sequence <- unlist(lapply(searchResult$sequence,
			memDecompress,
			type="gzip",
			asChar=TRUE))
		
		# replace characters that are not in the DNA_ALPHABET
		searchResult$sequence <- .Call("replaceChars",
			searchResult$sequence,
			replaceChar,
			type,
			PACKAGE="DECIPHER")
		
		# remove gaps if applicable
		if (removeGaps==2)
			searchResult$sequence <- .Call("replaceChar",
				searchResult$sequence,
				"-",
				"",
				PACKAGE="DECIPHER")
		else if (removeGaps==3)
			searchResult$sequence <- .Call("commonGaps",
				searchResult$sequence,
				PACKAGE="DECIPHER")
		
		# build an XStringSet based on the database
		if (type==1) {
			myXStringSet <- DNAStringSet(searchResult$sequence)
		} else if (type==2) {
			myXStringSet <- RNAStringSet(searchResult$sequence)
		} else if (type==3) {
			myXStringSet <- AAStringSet(searchResult$sequence)
		} else if (type==4) {
			myXStringSet <- BStringSet(searchResult$sequence)
		} else {
			searchResult$quality <- unlist(lapply(searchResult$quality,
				memDecompress,
				type="gzip",
				asChar=TRUE))
			if (type==5) {
				myXStringSet <- QualityScaledDNAStringSet(DNAStringSet(searchResult$sequence),
					PhredQuality(BStringSet(searchResult$quality)))
			} else if (type==6) {
				myXStringSet <- QualityScaledRNAStringSet(RNAStringSet(searchResult$sequence),
					PhredQuality(BStringSet(searchResult$quality)))
			} else if (type==7) {
				myXStringSet <- QualityScaledAAStringSet(AAStringSet(searchResult$sequence),
					PhredQuality(BStringSet(searchResult$quality)))
			} else { # type==8
				myXStringSet <- QualityScaledBStringSet(BStringSet(searchResult$sequence),
					PhredQuality(BStringSet(searchResult$quality)))
			}
		}
		names(myXStringSet) <- searchResult[, 1]
	}
	
	if (verbose) {
		time.2 <- Sys.time()
		if (countOnly) {
			cat("\n\nCount = ",
				count,
				" matches.\n",
				sep="")	
		} else {
			cat("\n\n",
				TYPES[type],
				" of length: ",
				length(myXStringSet),
				"\n",
				sep="")
		}
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	if (countOnly) {
		return(count)
	} else {
		return(myXStringSet)
	}
}