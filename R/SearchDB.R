SearchDB <- function(dbFile,
	tblName="Seqs",
	identifier="",
	type="XStringSet",
	limit=-1,
	replaceChar=NA,
	nameBy="row_names",
	orderBy="row_names",
	countOnly=FALSE,
	removeGaps="none",
	quality="Phred",
	clause="",
	processors=1,
	verbose=TRUE) {
	
	# error checking
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps[1], GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet",
		"QualityScaledDNAStringSet", "QualityScaledRNAStringSet", "QualityScaledAAStringSet", "QualityScaledBStringSet", "XStringSet", "QualityScaledXStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (length(identifier) != 1)
		stop("identifier must be a single character string.")
	if (!is.character(orderBy))
		stop("orderBy must be a character string.")
	if (!is.logical(countOnly))
		stop("countOnly must be a logical.")
	QUALS <- c("Phred", "Solexa", "Illumina")
	quality <- pmatch(quality[1], QUALS)
	if (is.na(quality))
		stop("Invalid quality.")
	if (quality == -1)
		stop("Ambiguous quality.")
	if (quality==1) {
		quality <- PhredQuality
	} else if (quality==2) {
		quality <- SolexaQuality
	} else if (quality==3) {
		quality <- IlluminaQuality
	}
	if (!is.character(clause))
		stop("clause must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors)!=processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	if (type > 4 && type != 9 && removeGaps > 1)
		stop(paste('removeGaps must be "none" when type is ', TYPES[type], '.', sep=''))
	if (is.numeric(limit)) {
		if (floor(limit)!=limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
	} else {
		if (!grepl("[0-9],[0-9]", limit, perl=TRUE)) {
			limit <- as.numeric(limit)
			if (floor(limit)!=limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
		}
	}
	if (limit > 0 && countOnly)
		stop("limit cannot be specified when countOnly is TRUE.")
	
	if (verbose)
		time.1 <- Sys.time()
	
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
	
	# build the search expression
	if (countOnly) {
		searchExpression <- paste('select count(*) from ',
			tblName,
			sep="")
	} else if (nameBy=="row_names" && orderBy=="row_names") {
		searchExpression <- paste('select row_names, sequence',
			ifelse(type > 4 && type != 9, ', quality', ''),
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
			ifelse(type > 4 && type != 9,
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
	
	if (identifier!="")
		searchExpression <- paste(searchExpression,
			' where identifier is "',
			identifier,
			'"',
			sep="")
	if (clause!="")
		searchExpression <- paste(searchExpression,
			ifelse(identifier=="", " where ", " and "),
			clause,
			sep="")
	
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
		cat("Search Expression:",
			strwrap(searchExpression,
					width=getOption("width") - 1L),
				sep="\n")
	
	rs <- dbSendQuery(dbConn, searchExpression)
	searchResult <- fetch(rs, n=-1)
	dbClearResult(rs)
	
	if (countOnly) {
		count <- as.integer(searchResult)
	} else {
		# decompress the resulting sequences
		searchResult$sequence <- Codec(searchResult$sequence,
			processors=processors)
		
		if (type==9 || type==10) {
			# guess the input type of XStringSet
			freqs <- .Call("composition",
				searchResult$sequence,
				PACKAGE="DECIPHER")
			
			if (all(freqs[1:2] < 0.6)) { # not DNA/RNA
				if (freqs[3] > 0.9) { # AA
					if (type==9) {
						type <- 3
					} else {
						type <- 7
					}
				} else {
					if (type==9) {
						type <- 4
					} else {
						type <- 8
					}
				}
			} else if (freqs[1] > freqs[2]) { # DNA
				if (type==9) {
					type <- 1
				} else {
					type <- 5
				}
			} else { # RNA
				if (type==9) {
					type <- 2
				} else {
					type <- 6
				}
			}
		}
		
		if (is.na(replaceChar)) {
			replaceChar <- NA_character_
		} else if (type==1 || type==5) {
			if (is.na(pmatch(replaceChar, DNA_ALPHABET)) && (replaceChar!=""))
				stop("replaceChar must be a character in the DNA_ALPHABET or empty character.")
		} else if (type==2 || type==6) {
			if (is.na(pmatch(replaceChar, RNA_ALPHABET)) && (replaceChar!=""))
				stop("replaceChar must be a character in the RNA_ALPHABET or empty character.")
		} else if (type==3 || type==7) {
			if (is.na(pmatch(replaceChar, AA_ALPHABET)) && (replaceChar!=""))
				stop("replaceChar must be a character in the AA_ALPHABET or empty character.")
		}
		
		if (type!=4 && type!=8) {
			# replace characters that are not in the DNA_ALPHABET
			searchResult$sequence <- .Call("replaceChars",
				searchResult$sequence,
				replaceChar,
				type,
				PACKAGE="DECIPHER")
		}
		
		# remove gaps if applicable
		if (removeGaps==2) {
			searchResult$sequence <- .Call("replaceChar",
				searchResult$sequence,
				"-",
				"",
				PACKAGE="DECIPHER")
			searchResult$sequence <- .Call("replaceChar",
				searchResult$sequence,
				".",
				"",
				PACKAGE="DECIPHER")
		} else if (removeGaps==3) {
			searchResult$sequence <- .Call("commonGaps",
				searchResult$sequence,
				PACKAGE="DECIPHER")
		}
		
		# build an XStringSet based on the database sequences
		if (type==1) {
			myXStringSet <- DNAStringSet(searchResult$sequence)
		} else if (type==2) {
			myXStringSet <- RNAStringSet(searchResult$sequence)
		} else if (type==3) {
			myXStringSet <- AAStringSet(searchResult$sequence)
		} else if (type==4) {
			myXStringSet <- BStringSet(searchResult$sequence)
		} else {
			searchResult$quality <- Codec(searchResult$quality)
			if (type==5) {
				myXStringSet <- QualityScaledDNAStringSet(DNAStringSet(searchResult$sequence),
					quality(searchResult$quality))
			} else if (type==6) {
				myXStringSet <- QualityScaledRNAStringSet(RNAStringSet(searchResult$sequence),
					quality(searchResult$quality))
			} else if (type==7) {
				myXStringSet <- QualityScaledAAStringSet(AAStringSet(searchResult$sequence),
					quality(searchResult$quality))
			} else { # type==8
				myXStringSet <- QualityScaledBStringSet(BStringSet(searchResult$sequence),
					quality(searchResult$quality))
			}
		}
		names(myXStringSet) <- searchResult[, 1]
	}
	
	if (verbose) {
		time.2 <- Sys.time()
		if (countOnly) {
			cat("\nCount = ",
				count,
				" matches.\n",
				sep="")	
		} else {
			cat("\n",
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
