BrowseSequences <- function(myDNAStringSet,
	htmlFile=paste(tempdir(),"/dna.html",sep=""),
	colorBases=TRUE,
	highlight=0,
	...) {
	# error checking
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	# check that the file exist
	if (is.character(htmlFile)) {
		htmlfile <- file(htmlFile, "w")
		on.exit(close(htmlfile))
	} else {
		stop("'htmlFile' must be a character string or connection.")
	}
	if (!is.logical(colorBases) && !is.numeric(colorBases))
		stop("colorBases must be a logical or numeric.")
	if (is.numeric(colorBases)) {
		if ((length(colorBases) %% 2) == 1 || length(colorBases) == 0)
			stop("colorBases must specify all start and endpoints.")
		if (any((colorBases[seq(2, length(colorBases), 2)] - colorBases[seq(1, length(colorBases), 2)]) <= 0))
			stop("colorBases specifies a negative or zero length range.")
		if (any(colorBases <= 0))
			stop("colorBases must be a positive numeric.")
		if (length(colorBases) > 2)
			if (any((colorBases[seq(3, length(colorBases), 2)] - colorBases[seq(2, length(colorBases) - 2, 2)]) <= 0))
				stop("Ranges specified in colorBases must be non-overlapping.")
	}
	if (is.null(names(myDNAStringSet)))
		names(myDNAStringSet) <- 1:length(myDNAStringSet)
	
	# convert the DNAStringSet to character strings
	html <- sapply(myDNAStringSet, toString)
	if (highlight %in% 1:length(myDNAStringSet)) {
		html <- sapply(html, strsplit, split="", fixed=TRUE)
		for (i in (1:length(html))[-highlight]) {
			l <- min(length(html[[highlight]]), length(html[[i]]))
			w <- which(html[[i]][1:l]==html[[highlight]][1:l])
			if (length(w) > 0)
				html[[i]][w] <- "."
		}
		html <- sapply(html, paste, collapse="")
	}
	
	# pad shorter sequences with periods
	maxW <- max(width(myDNAStringSet))
	for (i in 1:length(myDNAStringSet)) {
		html[i] <- paste(html[i],
			paste(rep("~", maxW - nchar(html[i])), collapse=""),
			sep="")
	}
	
	# add the lengths on the right
	alphabetTable <- alphabetFrequency(myDNAStringSet)
	myLengths <- (alphabetTable[,"A"] +
		alphabetTable[,"T"] +
		alphabetTable[,"G"] +
		alphabetTable[,"C"])
	
	# create a legend that gives position
	counter <- 20
	if (counter > max(width(myDNAStringSet)))
		counter <- max(width(myDNAStringSet))
	offset <- (counter - 1) - nchar(maxW)
	legend <- paste(paste(rep(" ", offset), collapse=""), format(seq(counter, maxW, by=counter)), collapse="")
	counter <- ceiling(counter/2)
	tickmarks <- paste(paste(rep("'", counter - 1), collapse=""), rep("|", floor(maxW/counter)), collapse="", sep="")
	tickmarks <- paste(tickmarks, paste(rep("'", maxW - counter*floor(maxW/counter)), collapse=""), sep="")
	
	if (is.numeric(colorBases) || colorBases) {
		# build a consensus sequence
		consensusSeq <- toString(ConsensusSequence(myDNAStringSet, ..., verbose=FALSE))
		
		if (is.numeric(colorBases)) {
			cs <- substr(consensusSeq, 0, colorBases[1] - 1)
			for (i in seq(1, length(colorBases), 2)) {
				csi <- substr(consensusSeq, colorBases[i], colorBases[i + 1])
				csi <- gsub("A", "<span class=\"A\">A</span>", csi)
				csi <- gsub("T", "<span class=\"T\">T</span>", csi)
				csi <- gsub("C", "<span class=\"C\">C</span>", csi)
				csi <- gsub("G", "<span class=\"G\">G</span>", csi)
				end <- ifelse(i==(length(colorBases) - 1),
					nchar(consensusSeq),
					colorBases[i + 2] - 1)
				cs <- paste(cs, csi, substr(consensusSeq, colorBases[i + 1] + 1, end), sep="")
			}
			consensusSeq <- paste(cs, maxW, sep="    ")
			
			htm <- substr(html, 0, colorBases[1] - 1)
			for (i in seq(1, length(colorBases), 2)) {
				htmi <- substr(html, colorBases[i], colorBases[i + 1])
				htmi <- gsub("A", "<span class=\"A\">A</span>", htmi)
				htmi <- gsub("T", "<span class=\"T\">T</span>", htmi)
				htmi <- gsub("C", "<span class=\"C\">C</span>", htmi)
				htmi <- gsub("G", "<span class=\"G\">G</span>", htmi)
				end <- ifelse(i==(length(colorBases) - 1),
					nchar(html),
					colorBases[i + 2] - 1)
				htm <- paste(htm, htmi, substr(html, colorBases[i + 1] + 1, end), sep="")
			}
			html <- paste(htm, myLengths, sep="    ")
		} else {
			consensusSeq <- gsub("A", "<span class=\"A\">A</span>", consensusSeq)
			consensusSeq <- gsub("T", "<span class=\"T\">T</span>", consensusSeq)
			consensusSeq <- gsub("C", "<span class=\"C\">C</span>", consensusSeq)
			consensusSeq <- gsub("G", "<span class=\"G\">G</span>", consensusSeq)
			consensusSeq <- paste(consensusSeq, maxW, sep="    ")
			
			# add color tags to each nucleotide
			html <- gsub("A", "<span class=\"A\">A</span>", html)
			html <- gsub("T", "<span class=\"T\">T</span>", html)
			html <- gsub("C", "<span class=\"C\">C</span>", html)
			html <- gsub("G", "<span class=\"G\">G</span>", html)
			html <- paste(html, myLengths, sep="    ")
		}
		
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