AlignDB <- function(dbFile,
	tblName="DNA",
	identifier="",
	type="DNAStringSet",
	add2tbl="DNA",
	batchSize=10000,
	perfectMatch=NULL,
	misMatch=NULL,
	gapOpening=NULL,
	gapExtension=NULL,
	terminalGap=-1,
	substitutionMatrix=NULL,
	processors=NULL,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	print(length(identifier)==2)
	if (!xor(length(identifier)==2, length(tblName)==2))
		stop("tblName or identifier must be length two.")
	if (length(identifier)==2 && identifier[1]==identifier[2])
		stop("The first and second identifier are the same.")
	if (length(tblName)==2 && tblName[1]==tblName[2])
		stop("The first and second tblName are the same.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.character(add2tbl))
		stop("add2tbl must be a character specifying a table name.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize)!=batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
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
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	if (type==1) { # DNAStringSet
		if (is.null(perfectMatch))
			perfectMatch <- 6
		if (is.null(misMatch))
			misMatch <- -3
		if (is.null(gapOpening))
			gapOpening <- -11
		if (is.null(gapExtension))
			gapExtension <- -3
	} else if (type==2) { # RNAStringSet
		if (is.null(perfectMatch))
			perfectMatch <- 8
		if (is.null(misMatch))
			misMatch <- 3
		if (is.null(gapOpening))
			gapOpening <- -9
		if (is.null(gapExtension))
			gapExtension <- -2
	} else { # AAStringSet
		if (is.null(perfectMatch))
			perfectMatch <- 4
		if (is.null(misMatch))
			misMatch <- 0
		if (is.null(gapOpening))
			gapOpening <- -5
		if (is.null(gapExtension))
			gapExtension <- -3
	}
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!all(is.numeric(terminalGap)))
		stop("terminalGap must be a numeric.")
	if (length(terminalGap) > 2 || length(terminalGap) < 1)
		stop("Length of terminalGap must be 1 or 2.")
	if (any(is.infinite(terminalGap)))
		stop("terminalGap must be finite.")
	if (length(terminalGap)==1)
		terminalGap[2] <- terminalGap[1]
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
	
	if (type==3) {
		if (is.null(substitutionMatrix)) {
			substitutionMatrix <- "BLOSUM62"
		} else if (is.character(substitutionMatrix)) {
			if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
		"PAM30", "PAM40", "PAM70", "PAM120", "PAM250")))
				stop("Invalid substitutionMatrix.")
		}
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.matrix(substitutionMatrix)) {
			if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix)))
		}
		subMatrix <- subMatrix[AAs, AAs]
	} else {
		if (!is.null(substitutionMatrix)) {
			if (is.matrix(substitutionMatrix)) {
				bases <- c("A", "C", "G", "T")
				if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
					any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
					stop("substitutionMatrix is incomplete.")
				substitutionMatrix <- substitutionMatrix[bases, bases]
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		}
	}
	
	count1 <- SearchDB(dbFile=dbConn,
		identifier=identifier[1],
		tblName=tblName[1],
		countOnly=TRUE,
		verbose=FALSE)
	if (count1==0)
		stop("No sequences found in dbFile matching the specified criteria.")
	count2 <- SearchDB(dbFile=dbConn,
		identifier=identifier[ifelse(length(identifier)==1, 1, 2)],
		tblName=tblName[ifelse(length(tblName)==1, 1, 2)],
		countOnly=TRUE,
		verbose=FALSE)
	if (count2==0)
		stop("No sequences found in dbFile matching the specified criteria.")
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(max=2*(count1 + count2), style=3)
		totSeqs <- 0L
	}
	
	for (i in 1:ceiling(count1/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[1],
			tblName=tblName[1],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			verbose=FALSE)
		
		temp <- unique(width(myXStringSet))
		if (length(temp) > 1)
			stop("Multiple width sequences found.")
		if (i==1) {
			uw1 <- temp
		} else {
			if (uw1!=temp)
				stop("Multiple width sequences found.")
		}
		
		temp <- .Call(ifelse(is(myXStringSet, "AAStringSet"),
					"consensusProfileAA",
					"consensusProfile"),
				myXStringSet,
				rep(1, length(myXStringSet)),
				PACKAGE="DECIPHER")
		if (i==1) {
			p.profile <- temp*length(myXStringSet)/count1
		} else {
			p.profile <- p.profile + temp*length(myXStringSet)/count1
		}
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	for (i in 1:ceiling(count2/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[ifelse(length(identifier)==1, 1, 2)],
			tblName=tblName[ifelse(length(tblName)==1, 1, 2)],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			verbose=FALSE)
		
		temp <- unique(width(myXStringSet))
		if (length(temp) > 1)
			stop("Multiple width sequences found.")
		if (i==1) {
			uw2 <- temp
			if (uw1*uw2 > 2147483647) # maximum when indexing by signed integer
			stop(paste("Alignment larger (",
				uw1*uw2,
				") than the maximum allowable size (2,147,483,647).",
				sep=""))
		} else {
			if (uw2!=temp)
				stop("Multiple width sequences found.")
		}
		
		temp <- .Call(ifelse(is(myXStringSet, "AAStringSet"),
					"consensusProfileAA",
					"consensusProfile"),
				myXStringSet,
				rep(1, length(myXStringSet)),
				PACKAGE="DECIPHER")
		if (i==1) {
			s.profile <- temp*length(myXStringSet)/count2
		} else {
			s.profile <- s.profile + temp*length(myXStringSet)/count2
		}
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	restrict <- min(perfectMatch,
		misMatch,
		gapOpening,
		gapExtension)*max(uw1, uw2) # unrestricted
	if (type==3) { # AAStringSet
		t <- .Call("alignProfilesAA",
			p.profile,
			s.profile,
			subMatrix,
			perfectMatch,
			misMatch,
			gapOpening,
			gapExtension,
			terminalGap[1],
			terminalGap[2],
			restrict,
			processors,
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		t <- .Call("alignProfiles",
			p.profile,
			s.profile,
			substitutionMatrix,
			perfectMatch,
			misMatch,
			gapOpening,
			gapExtension,
			terminalGap[1],
			terminalGap[2],
			restrict,
			processors,
			PACKAGE="DECIPHER")
	}
	
	inserts <- .align(t, uw1, uw2)
	if (length(identifier)==2) {
		id <- paste(identifier[1],
			identifier[2],
			sep="_")
	} else {
		id <- paste(tblName[1],
			tblName[2],
			sep="_")
	}
	
	for (i in 1:ceiling(count1/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[1],
			tblName=tblName[1],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			verbose=FALSE)
		
		myXStringSet <- .Call("insertGaps",
			myXStringSet,
			as.integer(inserts[[1]]),
			as.integer(inserts[[2]]),
			type,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- paste(names(myXStringSet),
			ifelse(length(identifier)==2, identifier[1], tblName[1]),
			sep="_")
		Seqs2DB(myXStringSet,
			TYPES[type],
			dbConn,
			id,
			add2tbl,
			verbose=FALSE)
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	for (i in 1:ceiling(count2/batchSize)) {
		myXStringSet <- SearchDB(dbFile=dbConn,
			identifier=identifier[ifelse(length(identifier)==1, 1, 2)],
			tblName=tblName[ifelse(length(tblName)==1, 1, 2)],
			type=TYPES[type],
			limit=paste((i - 1)*batchSize,
				",",
				batchSize,
				sep=""),
			verbose=FALSE)
		
		myXStringSet <- .Call("insertGaps",
			myXStringSet,
			as.integer(inserts[[3]]),
			as.integer(inserts[[4]]),
			type,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- paste(names(myXStringSet),
			ifelse(length(identifier)==2, identifier[2], tblName[2]),
			sep="_")
		Seqs2DB(myXStringSet,
			TYPES[type],
			dbConn,
			id,
			add2tbl,
			verbose=FALSE)
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	if (verbose) {
		close(pBar)
		cat("\nAdded ",
			count1 + count2,
			" aligned sequences to table ",
			add2tbl,
			" with identifier '",
			id,
			"'.\n\n",
			sep="")
	}
	
	invisible(count1 + count2)
}