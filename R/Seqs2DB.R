Seqs2DB <- function(seqs,
	type,
	dbFile,
	identifier,
	tblName="DNA",
	chunkSize=1e5,
	replaceTbl=FALSE,
	verbose=TRUE) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	SEQTYPES <- c("FASTA","GenBank","DNAStringSet")
	type <- pmatch(type, SEQTYPES)
	if (is.na(type))
		stop("Invalid seqs type.")
	if (type == -1)
		stop("Ambiguous seqs type.")
	if (length(type) > 1)
		stop("Too many seqs types provided - choose one.")
	if (!is.character(identifier))
		stop("Identifier must be a character!")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (substr(tblName, 1, 1) == "_")
		stop("Invalid tblName.")
	if (tblName == "taxa")
		stop("taxa is a reserved tblName.")
	if (!is.numeric(chunkSize))
		stop("chunkSize must be a numeric.")
	if (floor(chunkSize)!=chunkSize)
		stop("chunkSize must be a whole number.")
	if (chunkSize <= 0)
		stop("chunkSize must be greater than zero.")
	if (!is.logical(replaceTbl))
		stop("replaceTbl must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.character(seqs) &&
		(type==1 || type==2))
		stop("seqs must be a character string.")
	if (!is(seqs, "DNAStringSet") &&
		type==3)
		stop("seqs must be a DNAStringSet.")
	
	# initialize database
	driver = dbDriver("SQLite")
	numSeq <- integer(1)
	if (is.character(dbFile)) {
		dbConn = dbConnect(driver, dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn = dbFile
		if (!inherits(dbConn,"SQLiteConnection")) 
			stop("'dbFile' must be a character string or connection.")
		if (!isIdCurrent(dbConn))
			stop("The connection has expired.")
	}
	result <- dbListTables(dbConn)
	w <- which(result==tblName)
	if (length(w)==1 && length(which(result==paste("_", tblName, sep=""))) != 1)
		stop("Table is corrupted")
	f <- character(0)
	if (length(w)==1 && !replaceTbl) {
		searchExpression <- paste("select max(row_names) from ", tblName, sep="")
		numSeq <- as.integer(dbGetQuery(dbConn, searchExpression))
		
		searchExpression <- paste("select max(row_names) from _", tblName, sep="")
		if (as.integer(dbGetQuery(dbConn, searchExpression)) != numSeq)
			stop("Table is corrupted.")
		f <- dbListFields(dbConn, paste("_", tblName, sep=""))
		
		# make sure all the necessary fields exist
		colIDs <- c("row_names", "sequence")
		types <- c("INTEGER PRIMARY KEY ASC", "BLOB")
		for (i in 1:length(colIDs)) {
			if (is.na(match(colIDs[i], f))) {
				# first add the column if it does not already exist
				expression1 <- paste("alter table _",
					tblName,
					" add column ",
					colIDs[i],
					" ",
					fts[i],
					sep="")
				if (verbose)
					cat("Expression:  ", expression1, "\n", sep="")
				dbGetQuery(dbConn, expression1)
			}
		}
		
		f <- dbListFields(dbConn, tblName)
		
		# make sure all the necessary fields exist
		if (type==1 || type==3) {
			colIDs <- c("row_names", "id", "description")
			types <- c("INTEGER PRIMARY KEY ASC", "TEXT", "TEXT")
		} else if (type==2) {
			colIDs <- c("row_names", "id", "rank", "accession", "description")
			fts <- c("INTEGER PRIMARY KEY ASC", "TEXT", "TEXT", "TEXT", "TEXT")
		}
		for (i in 1:length(colIDs)) {
			if (is.na(match(colIDs[i], f))) {
				# first add the column if it does not already exist
				expression1 <- paste("alter table ",
					tblName,
					" add column ",
					colIDs[i],
					" ",
					fts[i],
					sep="")
				if (verbose)
					cat("Expression:  ", expression1, "\n", sep="")
				dbGetQuery(dbConn, expression1)
			}
		}
	} else {
		numSeq <- 0
		replaceTbl <- TRUE # necessary for field types
	}
	
	if (type==1) { # FASTA
		# scan in the FASTA file
		skipSize <- 1L
		s1 <- character(chunkSize)
		l <- 0L
		newSeqs <- 0
		buffer <- character()
		con <- file(seqs, "r")
		
		while (length(s1)==(chunkSize + l)) {
			# scan piece of file into memory
			if (verbose) {
				digits <- nchar(as.character(floor(skipSize/chunkSize)))
				cat("Reading FASTA file from line ",
					formatC(skipSize, digits=digits),
					" to ",
					formatC(skipSize + chunkSize - 1, digits=digits),
					"\n",
					sep="")
			}
			s1 <- c(buffer, readLines(con=con, n=chunkSize))
			skipSize <- skipSize + chunkSize
			
			# descriptions contains the line index of each sequence
			descriptions <- which(substr(s1, 1L, 1L) == ">")
			
			# create new vectors by adding or subtracting 1
			# from each line index in description
			dp <- descriptions + 1L
			dm <- descriptions - 1L
			end <- c(dm[-1], length(s1))
			
			# remove the last incomplete sequence unless it is the end of file
			l <- length(buffer)
			if (length(s1)==(chunkSize + l)) {
				buffer <- s1[descriptions[length(descriptions)]:length(s1)]
				length(descriptions) <- length(descriptions) - 1L
			} else {
				buffer <- character()
			}
			
			# numF contains the number of sequences for this iteration
			numF <- length(descriptions)
			if (numF == 0)
				stop("No FASTA sequences found, try increasing chunkSize.")
			
			# build the data frame sequence by sequence
			myData <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				id=identifier)
			
			myData$description <- substr(s1[descriptions],
				2L,
				nchar(s1[descriptions]))
			
			sequence <- lapply(seq_len(numF),
				function(i) {
					seq <- paste(s1[dp[i]:end[i]],
						collapse = "")
					seq <- gsub(" ",
						"",
						seq,
						fixed=TRUE)
				}
			)
			myData_ <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				sequence=I(lapply(sequence,
					memCompress,
					type="gzip")))
			
			# add database columns to the data frame
			if (length(f) > 0) {
				for (i in 1:length(f)) {
					if (is.na(match(f[i], names(myData)))) {
						d <- data.frame(rep(NA, numF))
						names(d) <- f[i]
						myData <- data.frame(myData, d)
					}
				}
			}
			
			# numSeq contains the total number of sequences so far
			numSeq <- numSeq + length(descriptions)
			newSeqs <- newSeqs + length(descriptions)
			
			if (replaceTbl) {
				ft <- list(row_names="INTEGER PRIMARY KEY ASC",
					id="TEXT",
					description="TEXT")
				ft_ <- list(row_names="INTEGER PRIMARY KEY ASC",
					sequence="BLOB")
			} else {
				ft <- ft_ <- NULL
			}
			
			dbWriteTable(dbConn,
				tblName,
				myData,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft)
			dbWriteTable(dbConn,
				paste("_", tblName, sep=""),
				myData_,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft_)
			
			replaceTbl <- FALSE
		}
		close(con)
	} else if (type==2) { # GenBank
		# scan in the GenBank file
		skipSize <- 1L
		s1 <- character(chunkSize)
		l <- 0L
		newSeqs <- 0
		buffer <- character()
		con <- file(seqs, "r")
		
		while (length(s1)==(chunkSize + l)) {
			# scan piece of file into memory
			if (verbose) {
				digits <- nchar(as.character(floor(skipSize/chunkSize)))
				cat("Reading GenBank file from line ",
					formatC(skipSize, digits=digits),
					" to ",
					formatC(skipSize + chunkSize - 1, digits=digits),
					"\n",
					sep="")
			}
			s1 <- c(buffer, readLines(con, n=chunkSize))
			skipSize <- skipSize + chunkSize
			
			# descriptions contains the line index of each sequence
			descriptions <- which(substr(s1, 1L, 10L) == "DEFINITION")
			rank <- which(substr(s1, 3L, 10L) == "ORGANISM") + 1
			accession <- which(substr(s1, 1L, 9L) == "ACCESSION")
			seq_start <- which(substr(s1, 1L, 6L) == "ORIGIN") + 1
			seq_end <- which(substr(s1, 1L, 2L) == "//") - 1
			
			# remove the last incomplete sequence unless it is the end of file
			l <- length(buffer)
			if (length(s1)==(chunkSize + l)) {
				buffer <- s1[descriptions[length(descriptions)]:length(s1)]
				length(descriptions) <- length(descriptions) - 1L
			} else {
				buffer <- character()
			}
			
			# numF contains the number of sequences for this iteration
			numF <- length(descriptions)
			if (numF == 0)
				stop("No GenBank records found, try increasing chunkSize.")
			
			# build the data frame sequence by sequence
			myData <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
					id=identifier)
			
			myData$description <- substr(s1[descriptions],
				13L,
				nchar(s1[descriptions]))
			
			myData$accession <- substr(s1[accession[1:length(descriptions)]],
				13L,
				nchar(s1[accession[1:length(descriptions)]]))
			
			sequence <- lapply(seq_len(numF),
				function(i) {
					seq <- paste(
						s1[seq_start[i]:seq_end[i]],
						collapse = "")
					seq <- gsub("[0-9]|[[:space:]]",
						"",
						seq,
						perl=TRUE)
				}
			)
			myData_ <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				sequence=I(lapply(sequence,
					memCompress,
					type="gzip")))
			
			for (i in 1:numF) {
				# parse the ORGANISM line into a rank field
				j <- 0L
				lineage <- FALSE
				if (is.na(rank[i])) {
					myData$rank[i] <- "unknown"
					next
				}
				while (substr(s1[rank[i] + j], 1L, 1L)==" ") {
					if (lineage || grepl(";",
							substr(s1[rank[i] + j],
							13L,
							nchar(s1[rank[i] + j])),
							fixed=TRUE)) {
						if (!lineage)
							myData$rank[i] <- ""
						lineage <- TRUE
						myData$rank[i] <- paste(myData$rank[i],
							substr(s1[rank[i] + j],
							13L,
							nchar(s1[rank[i] + j])),
							sep="")
					}
					j <- j + 1L
				}
			}
			
			# add database columns to the data frame
			if (length(f) > 0) {
				for (i in 1:length(f)) {
					if (is.na(match(f[i], names(myData)))) {
						d <- data.frame(rep(NA, numF))
						names(d) <- f[i]
						myData <- data.frame(myData, d)
					}
				}
			}
			
			# numSeq contains the total number of sequences so far
			numSeq <- numSeq + length(descriptions)
			newSeqs <- newSeqs + length(descriptions)
			
			if (replaceTbl) {
				ft <- list(row_names="INTEGER PRIMARY KEY ASC",
					id="TEXT",
					description="TEXT",
					accession="TEXT",
					rank="TEXT")
				ft_ <- list(row_names="INTEGER PRIMARY KEY ASC",
					sequence="BLOB")
			} else {
				ft <- ft_ <- NULL
			}
			
			dbWriteTable(dbConn,
				tblName,
				myData,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft)
			dbWriteTable(dbConn,
				paste("_", tblName, sep=""),
				myData_,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft_)
			
			replaceTbl <- FALSE
		}
		close(con)
	} else if (type==3) { # DNAStringSet
		# add the sequences to the database
		newSeqs <- length(seqs)
		
		if (verbose)
			cat("Adding", newSeqs, "sequences to the database...\n")
		
		# give the seqs names if they do not have any
		if (is.null(names(seqs)))
			names(seqs) <- (numSeq + 1):(numSeq + newSeqs)
		
		# build the data frame sequence by sequence
		myData <- data.frame(row_names=seq(from=(numSeq + 1),
			to=(numSeq + newSeqs)),
			id=identifier,
			description=names(seqs))
		
		myData_ <- data.frame(row_names=seq(from=(numSeq + 1),
			to=(numSeq + newSeqs)),
			sequence=I(lapply(lapply(seqs,
					toString),
				memCompress,
				type="gzip")))
		
		# add database columns to the data frame
		if (length(f) > 0) {
			for (i in 1:length(f)) {
				if (is.na(match(f[i], names(myData)))) {
					d <- data.frame(rep(NA, newSeqs))
					names(d) <- f[i]
					myData <- data.frame(myData, d)
				}
			}
		}
		
		if (replaceTbl) {
			ft <- list(row_names="INTEGER PRIMARY KEY ASC",
				id="TEXT",
				description="TEXT")
			ft_ <- list(row_names="INTEGER PRIMARY KEY ASC",
				sequence="BLOB")
		} else {
			ft <- ft_ <- NULL
		}
		
		dbWriteTable(dbConn,
			tblName,
			myData,
			row.names=FALSE,
			overwrite=replaceTbl,
			append=!replaceTbl,
			field.types=ft)
		dbWriteTable(dbConn,
			paste("_", tblName, sep=""),
			myData_,
			row.names=FALSE,
			overwrite=replaceTbl,
			append=!replaceTbl,
			field.types=ft_)
	}
	
	searchExpression <- paste("select count(*) from ",
		tblName,
		sep="")
	numSeq <- as.integer(dbGetQuery(dbConn, searchExpression))
	
	if (verbose) { # print the elapsed time to import
		time.2 <- Sys.time()
		if (newSeqs!=numSeq)
			cat("\n",
				"Added ",
				newSeqs,
				" new sequences to the table ",
				tblName,
				".",
				sep="")
		cat("\n",
			numSeq,
			" total sequences in the table ",
			tblName,
			".",
			"\n",
			sep="")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	return(numSeq)
}