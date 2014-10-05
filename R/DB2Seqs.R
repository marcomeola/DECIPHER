DB2Seqs <- function(file,
	dbFile,
	tblName="DNA",
	identifier="",
	type="BStringSet",
	limit=-1,
	replaceChar="-",
	nameBy="description",
	orderBy="row_names",
	removeGaps="none",
	append=FALSE,
	width=80,
	chunkSize=1e5,
	clause="",
	verbose=TRUE) {
	
	# error checking
	if (!is.character(file))
		stop("file must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet",
		"QualityScaledDNAStringSet", "QualityScaledRNAStringSet", "QualityScaledAAStringSet", "QualityScaledBStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(limit))
		stop("limit must be a numeric.")
	if (floor(limit)!=limit)
		stop("limit must be a whole number.")
	if (!is.character(nameBy))
		stop("nameBy must be a character string.")
	if (!is.character(orderBy))
		stop("orderBy must be a character string.")
	if (!is.logical(append))
		stop("append must be a logical.")
	if (!is.character(clause))
		stop("clause must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(chunkSize))
		stop("chunkSize must be a numeric.")
	if (floor(chunkSize)!=chunkSize)
		stop("chunkSize must be a whole number.")
	if (chunkSize <= 0)
		stop("chunkSize must be greater than zero.")
	if (!is.numeric(width))
		stop("width must be a numeric.")
	if (floor(width)!=width)
		stop("width must be a whole number.")
	if (width < 1)
		stop("width must be at least 1.")
	if (width > 20001)
		stop("width must be at most 20001.")
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps, GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
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
	
	if (verbose)
		time.1 <- Sys.time()
	
	searchExpression <- tblName
	
	if (identifier!="")
		searchExpression <- paste(searchExpression,
			' where id like "',
			identifier,
			'"',
			sep="")
	if (clause!="")
		searchExpression <- paste(searchExpression,
			clause)
	
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
	
	if (count < 1)
		stop("No sequences matched the specified parameters.")
	if (limit!=-1 && count > limit)
		count <- limit
	if (count > chunkSize && removeGaps==3)
		stop("chunkSize must be at least ",
			count,
			" with this query if removeGaps is 'common'.")
	
	if (orderBy!="row_names") # default ordering is row_names
		searchExpression <- paste(searchExpression,
			'order by',
			orderBy)
	
	if (verbose)
		pBar <- txtProgressBar(style=3)
	s <- seq(1, count, chunkSize)
	for (i in 1:length(s)) {
		# build the search expression
		searchExpression1 <- paste(searchExpression,
			'limit',
			ifelse(i==length(s), count - s[i] + 1, chunkSize),
			'offset',
			ifelse(i==1, 0, s[i] - 1))
		if (all(nameBy=="row_names")) {
			searchExpression1 <- paste("select row_names from",
				searchExpression1)
		} else {
			temp <- 'select row_names'
			for (j in 1:length(nameBy)) {
				if (nameBy[j]=="row_names")
					next
				temp <- paste(temp, ", ", nameBy[j], sep="")
			}
			searchExpression1 <- paste(temp,
				"from",
				searchExpression1)
		}
		
		rs <- dbSendQuery(dbConn, searchExpression1)
		searchResult <- fetch(rs, n=-1)
		dbClearResult(rs)
		
		searchExpression2 <- paste(ifelse(type > 4,
				'select row_names, sequence, quality from _',
				'select row_names, sequence from _'),
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
		
		if (type!=4 && type!=8) {
			# replace characters that are not in the DNA_ALPHABET
			searchResult2$sequence <- .Call("replaceChars",
				searchResult2$sequence,
				replaceChar,
				type,
				PACKAGE="DECIPHER")
		}
		
		# remove gaps if applicable
		if (removeGaps==2) {
			searchResult2$sequence <- .Call("replaceChar",
				searchResult2$sequence,
				"-",
				"",
				PACKAGE="DECIPHER")
			searchResult2$sequence <- .Call("replaceChar",
				searchResult2$sequence,
				".",
				"",
				PACKAGE="DECIPHER")
		} else if (removeGaps==3) {
			searchResult2$sequence <- .Call("commonGaps",
				searchResult2$sequence,
				PACKAGE="DECIPHER")
		}
		
		if (type > 4) {
			# decompress the resulting qualities
			searchResult2$quality <- lapply(searchResult2$quality,
				memDecompress,
				type="gzip",
				asChar=TRUE)
			searchResult2$quality <- paste(searchResult2$quality)
		}
		
		if (type==1 || type==5) {
			myXStringSet <- DNAStringSet(searchResult2$sequence)
		} else if (type==2 || type==6) {
			myXStringSet <- RNAStringSet(searchResult2$sequence)
		} else if (type==3 || type==7) {
			myXStringSet <- AAStringSet(searchResult2$sequence)
		} else if (type==4 || type==8) {
			myXStringSet <- BStringSet(searchResult2$sequence)
		}
		
		if (length(nameBy) > 1) {
			names(myXStringSet) <- do.call(paste,
				c(searchResult[, nameBy],
					sep="::"))
		} else {
			names(myXStringSet) <- searchResult[, nameBy]
		}
		
		if (type > 4) {
			writeXStringSet(myXStringSet,
				file,
				append=ifelse(i==1, append, TRUE),
				format="FASTQ",
				qualities=PhredQuality(searchResult2$quality))
		} else {
			writeXStringSet(myXStringSet,
				file,
				append=ifelse(i==1, append, TRUE),
				format="FASTA",
				width=width)
		}
		
		if (verbose)
			setTxtProgressBar(pBar, ifelse(i==length(s), 1, (s[i + 1] - 1)/count))
	}
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\n\nWrote ",
			count,
			ifelse(count > 1, " sequences.", " sequence."),
			"\n",
			sep="")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	invisible(count)
}