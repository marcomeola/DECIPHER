AlignDB <- function(dbFile,
	tblName="Seqs",
	identifier="",
	type="DNAStringSet",
	add2tbl="Seqs",
	batchSize=10000,
	perfectMatch=5,
	misMatch=0,
	gapOpening=-13,
	gapExtension=-1,
	gapPower=-1,
	terminalGap=-5,
	normPower=1,
	substitutionMatrix=NULL,
	processors=1,
	verbose=TRUE,
	...) {
	
	# error checking
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
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
	if (!all(is.numeric(normPower)))
		stop("normPower must be a numeric.")
	if (normPower < 0)
		stop("normPower must be at least zero.")
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
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	gapOpening <- gapOpening/2 # split into gap opening and closing
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(gapPower))
		stop("gapPower must be a numeric.")
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
	
	if (type==3L) { # AAStringSet
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.null(substitutionMatrix)) {
			# use MIQS
			substitutionMatrix <- matrix(c(3.2,-1.3,-0.4,-0.4,1.5,-0.2,-0.4,0.4,-1.2,-1.3,-1.4,-0.7,-1,-2.3,-0.1,0.8,0.8,-3.6,-2.4,0,-6.1,-1.3,6.2,-0.1,-1.5,-2.7,1.8,-0.7,-1.9,0.9,-2.4,-2.5,3.3,-1.1,-3.3,-1.1,-0.3,-0.9,-3.8,-1.9,-2.3,-6.1,-0.4,-0.1,5.1,2.6,-1.6,0.9,0.8,0.2,1,-3.6,-3.5,0.7,-2.3,-3.5,-1.4,0.9,0,-4.5,-1.5,-2.6,-6.1,-0.4,-1.5,2.6,5.7,-3.7,0.9,2.7,-0.5,0.3,-4.5,-4.6,0.4,-3.3,-5.8,-0.3,0.3,-0.2,-5.3,-3.9,-3.5,-6.1,1.5,-2.7,-1.6,-3.7,11.7,-2.8,-3.2,-1.7,-1.2,0.2,-2.3,-3.2,0.1,-2.8,-2.8,1,0,-6.1,-0.7,1.8,-6.1,-0.2,1.8,0.9,0.9,-2.8,3.6,2.1,-1.6,1.2,-2.2,-1.9,1.7,-0.4,-2.4,-0.4,0.4,0.1,-5.4,-2.8,-1.8,-6.1,-0.4,-0.7,0.8,2.7,-3.2,2.1,4.3,-1.3,-0.2,-3.3,-2.8,1.1,-2.3,-4.1,0,0.4,-0.2,-5.8,-2.4,-2.3,-6.1,0.4,-1.9,0.2,-0.5,-1.7,-1.6,-1.3,7.6,-1.6,-5.4,-4.8,-1.7,-3.6,-4.6,-1.6,0,-1.9,-4.8,-4.5,-3.8,-6.1,-1.2,0.9,1,0.3,-1.2,1.2,-0.2,-1.6,7.5,-2.2,-1.9,0,-2.1,0,-1.5,0,-0.2,-0.3,2.1,-2.3,-6.1,-1.3,-2.4,-3.6,-4.5,0.2,-2.2,-3.3,-5.4,-2.2,4.6,3.1,-2.3,1.7,0.7,-3.7,-2.8,-0.7,-0.7,-0.8,3.3,-6.1,-1.4,-2.5,-3.5,-4.6,-2.3,-1.9,-2.8,-4.8,-1.9,3.1,4.6,-2.4,3.2,2.1,-2.8,-2.9,-1.6,-0.2,0,2,-6.1,-0.7,3.3,0.7,0.4,-3.2,1.7,1.1,-1.7,0,-2.3,-2.4,3.6,-1.1,-3.7,-0.1,0,0,-4,-2.3,-2,-6.1,-1,-1.1,-2.3,-3.3,0.1,-0.4,-2.3,-3.6,-2.1,1.7,3.2,-1.1,5.4,1.4,-2.8,-1.8,-0.8,-2.1,-0.9,1.4,-6.1,-2.3,-3.3,-3.5,-5.8,-2.8,-2.4,-4.1,-4.6,0,0.7,2.1,-3.7,1.4,7.4,-3.7,-2.6,-2.3,4.2,5.2,-0.3,-6.1,-0.1,-1.1,-1.4,-0.3,-2.8,-0.4,0,-1.6,-1.5,-3.7,-2.8,-0.1,-2.8,-3.7,8.4,-0.1,-0.5,-3.6,-4.5,-2.5,-6.1,0.8,-0.3,0.9,0.3,1,0.4,0.4,0,0,-2.8,-2.9,0,-1.8,-2.6,-0.1,3.1,1.6,-3.5,-1.5,-1.4,-6.1,0.8,-0.9,0,-0.2,0,0.1,-0.2,-1.9,-0.2,-0.7,-1.6,0,-0.8,-2.3,-0.5,1.6,3.8,-5.3,-2.1,-0.1,-6.1,-3.6,-3.8,-4.5,-5.3,-6.1,-5.4,-5.8,-4.8,-0.3,-0.7,-0.2,-4,-2.1,4.2,-3.6,-3.5,-5.3,14.8,4.9,-3.3,-6.1,-2.4,-1.9,-1.5,-3.9,-0.7,-2.8,-2.4,-4.5,2.1,-0.8,0,-2.3,-0.9,5.2,-4.5,-1.5,-2.1,4.9,8.3,-1.2,-6.1,0,-2.3,-2.6,-3.5,1.8,-1.8,-2.3,-3.8,-2.3,3.3,2,-2,1.4,-0.3,-2.5,-1.4,-0.1,-3.3,-1.2,3.5,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,1),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (is.character(substitutionMatrix)) {
			if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
		"PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")))
				stop("Invalid substitutionMatrix.")
		}
		if (is.matrix(substitutionMatrix)) {
			if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix, envir=environment())))
		}
		subMatrix <- subMatrix[AAs, AAs]
		subMatrix <- as.numeric(subMatrix)
	} else {
		if (!is.null(substitutionMatrix)) {
			if (is.matrix(substitutionMatrix)) {
				bases <- c("A", "C", "G",
					ifelse(type==2L, "U", "T"))
				if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
					any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
					stop("substitutionMatrix is incomplete.")
				substitutionMatrix <- substitutionMatrix[bases, bases]
				substitutionMatrix <- as.numeric(substitutionMatrix)
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		}
	}
	
	count1 <- SearchDB(dbFile=dbConn,
		identifier=identifier[1],
		tblName=tblName[1],
		countOnly=TRUE,
		processors=processors,
		verbose=FALSE)
	if (count1==0)
		stop("No sequences found in dbFile matching the specified criteria.")
	count2 <- SearchDB(dbFile=dbConn,
		identifier=identifier[ifelse(length(identifier)==1, 1, 2)],
		tblName=tblName[ifelse(length(tblName)==1, 1, 2)],
		countOnly=TRUE,
		processors=processors,
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
			processors=processors,
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
		
		if (type==3L) {
			temp <- .Call("consensusProfileAA",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		} else {
			temp <- .Call("consensusProfile",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		}
		
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
			processors=processors,
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
		
		if (type==3L) {
			temp <- .Call("consensusProfileAA",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		} else {
			temp <- .Call("consensusProfile",
				myXStringSet,
				rep(1, length(myXStringSet)),
				NULL,
				PACKAGE="DECIPHER")
		}
		
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
		inserts <- .Call("alignProfilesAA",
			p.profile,
			s.profile,
			subMatrix,
			numeric(),
			gapOpening,
			gapExtension,
			gapPower,
			normPower,
			terminalGap[1],
			terminalGap[2],
			restrict,
			processors,
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		inserts <- .Call("alignProfiles",
			p.profile,
			s.profile,
			substitutionMatrix,
			numeric(),
			perfectMatch,
			misMatch,
			gapOpening,
			gapExtension,
			gapPower,
			normPower,
			terminalGap[1],
			terminalGap[2],
			restrict,
			processors,
			PACKAGE="DECIPHER")
	}
	
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
			processors=processors,
			verbose=FALSE)
		ns <- names(myXStringSet)
		
		myXStringSet <- .Call("insertGaps",
			myXStringSet,
			as.integer(inserts[[1]]),
			as.integer(inserts[[2]]),
			type,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- paste(ns,
			ifelse(length(identifier)==2, identifier[1], tblName[1]),
			sep="_")
		
		Seqs2DB(myXStringSet,
			TYPES[type],
			dbConn,
			id,
			add2tbl,
			processors=processors,
			verbose=FALSE,
			...)
		
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
			processors=processors,
			verbose=FALSE)
		ns <- names(myXStringSet)
		
		myXStringSet <- .Call("insertGaps",
			myXStringSet,
			as.integer(inserts[[3]]),
			as.integer(inserts[[4]]),
			type,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- paste(ns,
			ifelse(length(identifier)==2, identifier[2], tblName[2]),
			sep="_")
		
		Seqs2DB(myXStringSet,
			TYPES[type],
			dbConn,
			id,
			add2tbl,
			processors=processors,
			verbose=FALSE,
			...)
		
		if (verbose) {
			totSeqs <- totSeqs + length(myXStringSet)
			setTxtProgressBar(pBar,
				totSeqs)
		}
	}
	
	if (verbose) {
		close(pBar)
		cat(strwrap(paste("\nAdded ",
					count1 + count2,
					" aligned sequences to table ",
					add2tbl,
					" with identifier '",
					id,
					"'.\n\n",
					sep="",
					collapse=""),
				width=getOption("width") - 1L),
			sep="\n")
	}
	
	invisible(count1 + count2)
}
