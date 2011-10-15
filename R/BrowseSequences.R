BrowseSequences <- function(myDNAStringSet,
	htmlFile=paste(tempdir(),"/dna.html",sep=""),
	colorBases=FALSE,
	...) {
	# error checking
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	# check that the file exist
	if (is.character(htmlFile)) {
		htmlfile <- file(htmlFile, "w")
		on.exit(close(htmlfile))
	} else {
		stop("'htmlFile' must be a character string or connection")
	}
	
	# convert the DNAStringSet to character strings
	html <- sapply(myDNAStringSet, toString)
	
	# pad shorter sequences with periods
	maxW <- max(width(myDNAStringSet))
	for (i in 1:length(myDNAStringSet)) {
		html[i] <- paste(html[i],
			paste(rep("~",maxW - nchar(html[i])), collapse=""),
			sep="")
	}
	
	# add the lengths on the right
	alphabetTable <- alphabetFrequency(myDNAStringSet)
	myLengths <- (alphabetTable[,"A"] +
		alphabetTable[,"T"] +
		alphabetTable[,"G"] +
		alphabetTable[,"C"])
	html <- paste(html, myLengths, sep="    ")
	
	# create a legend that gives position
	counter <- 20
	offset <- (counter - 1) - nchar(maxW)
	legend <- paste(paste(rep(" ", offset), collapse=""), format(seq(counter, maxW, by=counter)), collapse="")
	counter <- counter/2
	tickmarks <- paste(paste(rep("'", counter - 1), collapse=""), rep("|", floor(maxW/counter)), collapse="", sep="")
	
	if (colorBases) {
		# build a consensus sequence
		consensusSeq <- paste(toString(ConsensusSequence(myDNAStringSet, ..., verbose=FALSE)), maxW, sep="    ")
		consensusSeq <- gsub("A", "<span class=\"A\">A</span>", consensusSeq)
		consensusSeq <- gsub("T", "<span class=\"T\">T</span>", consensusSeq)
		consensusSeq <- gsub("C", "<span class=\"C\">C</span>", consensusSeq)
		consensusSeq <- gsub("G", "<span class=\"G\">G</span>", consensusSeq)
		
		# add color tags to each nucleotide
		html <- gsub("A", "<span class=\"A\">A</span>", html)
		html <- gsub("T", "<span class=\"T\">T</span>", html)
		html <- gsub("C", "<span class=\"C\">C</span>", html)
		html <- gsub("G", "<span class=\"G\">G</span>", html)
		
		# add the legend and consensus sequence
		html <- paste(format(c("", "", names(myDNAStringSet), "", "Consensus"), justify="right"),		c(legend, tickmarks, html, "", consensusSeq), sep="    ")
		
		styles <- paste("<style type=text/css>",
			"span.A {background: #1E90FF; color: #FFF;}",
			"span.C {background: #32CD32; color: #FFF;}",
			"span.G {background: #9400D3; color: #FFF;}",
			"span.T {background: #000; color: #FFF;}",
			"</style>")
		html <- c("<html>",styles,"<pre>", html, "</pre></html>")
	} else {
		# build a consensus sequence
		consensusSeq <- paste(toString(ConsensusSequence(myDNAStringSet, ..., verbose=FALSE)), maxW, sep="    ")
		
		# add the legend and consensus sequence
		html <- paste(format(c("", "", names(myDNAStringSet), "", "Consensus"), justify="right"),		c(legend, tickmarks, html, "", consensusSeq), sep="    ")
		
		html <- c("<html>","<pre>", html, "</pre></html>")
	}
	writeLines(html, htmlfile)
	browseURL(htmlFile)
	return(TRUE)
}