IdLengths <- function(dbFile,
	tblName="DNA",
	identifier="",
	add2tbl=FALSE,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
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
			stop("'dbFile' must be a character string or SQLiteConnection")
		if (!isIdCurrent(dbConn))
			stop("The connection has expired.")
	}
	
	count <- SearchDB(dbFile=dbFile,
		identifier=identifier,
		tblName=tblName,
		countOnly=TRUE,
		verbose=FALSE)
	
	if (count > 10000) {
		
		# initialize a progress bar
		if (verbose)
			pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
		
		for (i in 1:ceiling(count/10000)) {
			myDNAStringSet <- SearchDB(dbFile=dbFile,
				identifier=identifier,
				tblName=tblName,
				limit=paste((i - 1)*10000,
					",",
					10000,
					sep=""),
				verbose=FALSE)
			
			numF <- length(myDNAStringSet)
			
			alphabetTable <- alphabetFrequency(myDNAStringSet)
			
			if (i==1) {
				lengths <- data.frame(bases=integer(numF),
					nonbases=integer(numF),
					width=integer(numF),
					row.names=names(myDNAStringSet))
				
				lengths$bases = (alphabetTable[,"A"] +
					alphabetTable[,"T"] +
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
				lengths$width <- width(myDNAStringSet)
					
				if (is.character(add2tbl) || add2tbl)
					Add2DB(myData=lengths,
						dbFile=dbFile,
						tblName=ifelse(is.character(add2tbl),add2tbl,tblName),
						verbose=FALSE)
			} else {
				lengths_temp <- data.frame(bases=integer(numF),
					nonbases=integer(numF),
					width=integer(numF),
					row.names=names(myDNAStringSet))
				
				lengths_temp$bases = (alphabetTable[,"A"] +
					alphabetTable[,"T"] +
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
				lengths_temp$width <- width(myDNAStringSet)
				
				lengths <- rbind(lengths, lengths_temp)
				
				if (is.character(add2tbl) || add2tbl)
					Add2DB(myData=lengths_temp,
						dbFile=dbFile,
						tblName=ifelse(is.character(add2tbl),add2tbl,tblName),
						verbose=FALSE)
			}
			
			if (verbose)
				setTxtProgressBar(pBar,
					floor(100*i/ceiling(count/10000)))
		}
	} else {
		myDNAStringSet <- SearchDB(dbFile=dbFile,
			identifier=identifier,
			tblName=tblName,
			verbose=FALSE)
		
		numF <- length(myDNAStringSet)
		
		alphabetTable <- alphabetFrequency(myDNAStringSet)
		lengths <- data.frame(bases=integer(numF),
			nonbases=integer(numF),
			width=integer(numF),
			row.names=names(myDNAStringSet))
		
		lengths$bases = (alphabetTable[,"A"] +
			alphabetTable[,"T"] +
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
		lengths$width <- width(myDNAStringSet)
		
		if (is.character(add2tbl) || add2tbl)
			Add2DB(myData=lengths,
				dbFile=dbFile,
				tblName=ifelse(is.character(add2tbl),add2tbl,tblName),
				verbose=FALSE)
	}
	
	if (verbose) {
		cat("\nLengths counted for ",count," sequences.", sep="")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to ",
				ifelse(is.character(add2tbl),add2tbl,tblName),
				":  \"bases\", \"nonbases\", and \"width\".",
				sep="")
		cat("\n\n")
	}
	
	return(lengths)
}