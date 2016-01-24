IdLengths <- function(dbFile,
	tblName="Seqs",
	identifier="",
	type="DNAStringSet",
	add2tbl=FALSE,
	batchSize=10000,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (type > 2)
		stop('type must be one of either "DNAStringSet" or "RNAStringSet".')
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize)!=batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
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
	
	count <- SearchDB(dbFile=dbFile,
		identifier=identifier,
		tblName=tblName,
		countOnly=TRUE,
		processors=processors,
		verbose=FALSE)
	
	if (count > batchSize) {
		
		# initialize a progress bar
		if (verbose)
			pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
		
		for (i in 1:ceiling(count/batchSize)) {
			myXStringSet <- SearchDB(dbFile=dbFile,
				identifier=identifier,
				tblName=tblName,
				type=TYPES[type],
				limit=paste((i - 1)*batchSize,
					",",
					batchSize,
					sep=""),
				processors=processors,
				verbose=FALSE)
			
			numF <- length(myXStringSet)
			
			alphabetTable <- alphabetFrequency(myXStringSet)
			
			if (i==1) {
				lengths <- data.frame(bases=integer(numF),
					nonbases=integer(numF),
					width=integer(numF),
					row.names=names(myXStringSet))
				
				lengths$bases = (alphabetTable[,"A"] +
					alphabetTable[,ifelse(type==1, "T", "U")] +
					alphabetTable[,"G"] +
					alphabetTable[,"C"])
				lengths$nonbases = (alphabetTable[,"M"] +
					alphabetTable[,"R"] +
					alphabetTable[,"W"] +
					alphabetTable[,"S"] +
					alphabetTable[,"Y"] +
					alphabetTable[,"K"] +
					alphabetTable[,"V"] +
					alphabetTable[,"H"] +
					alphabetTable[,"D"] +
					alphabetTable[,"B"] +
					alphabetTable[,"N"])
				lengths$width <- width(myXStringSet)
					
				if (is.character(add2tbl) || add2tbl)
					Add2DB(myData=lengths,
						dbFile=dbFile,
						tblName=ifelse(is.character(add2tbl),add2tbl,tblName),
						verbose=FALSE)
			} else {
				lengths_temp <- data.frame(bases=integer(numF),
					nonbases=integer(numF),
					width=integer(numF),
					row.names=names(myXStringSet))
				
				lengths_temp$bases = (alphabetTable[,"A"] +
					alphabetTable[,ifelse(type==1, "T", "U")] +
					alphabetTable[,"G"] +
					alphabetTable[,"C"])
				lengths_temp$nonbases = (alphabetTable[,"M"] +
					alphabetTable[,"R"] +
					alphabetTable[,"W"] +
					alphabetTable[,"S"] +
					alphabetTable[,"Y"] +
					alphabetTable[,"K"] +
					alphabetTable[,"V"] +
					alphabetTable[,"H"] +
					alphabetTable[,"D"] +
					alphabetTable[,"B"] +
					alphabetTable[,"N"])
				lengths_temp$width <- width(myXStringSet)
				
				lengths <- rbind(lengths, lengths_temp)
				
				if (is.character(add2tbl) || add2tbl)
					Add2DB(myData=lengths_temp,
						dbFile=dbFile,
						tblName=ifelse(is.character(add2tbl),add2tbl,tblName),
						verbose=FALSE)
			}
			
			if (verbose)
				setTxtProgressBar(pBar,
					floor(100*i/ceiling(count/batchSize)))
		}
	} else {
		myXStringSet <- SearchDB(dbFile=dbFile,
			identifier=identifier,
			tblName=tblName,
			type=TYPES[type],
			processors=processors,
			verbose=FALSE)
		
		numF <- length(myXStringSet)
		
		alphabetTable <- alphabetFrequency(myXStringSet)
		lengths <- data.frame(bases=integer(numF),
			nonbases=integer(numF),
			width=integer(numF),
			row.names=names(myXStringSet))
		
		lengths$bases = (alphabetTable[,"A"] +
			alphabetTable[,ifelse(type==1, "T", "U")] +
			alphabetTable[,"G"] +
			alphabetTable[,"C"])
		lengths$nonbases = (alphabetTable[,"M"] +
			alphabetTable[,"R"] +
			alphabetTable[,"W"] +
			alphabetTable[,"S"] +
			alphabetTable[,"Y"] +
			alphabetTable[,"K"] +
			alphabetTable[,"V"] +
			alphabetTable[,"H"] +
			alphabetTable[,"D"] +
			alphabetTable[,"B"] +
			alphabetTable[,"N"])
		lengths$width <- width(myXStringSet)
		
		if (is.character(add2tbl) || add2tbl)
			Add2DB(myData=lengths,
				dbFile=dbFile,
				tblName=ifelse(is.character(add2tbl),add2tbl,tblName),
				verbose=FALSE)
	}
	
	if (verbose) {
		cat("\nLengths counted for ", count, " sequences.", sep="")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				":  \"bases\", \"nonbases\", and \"width\".",
				sep="")
		cat("\n\n")
		
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(lengths)
}
