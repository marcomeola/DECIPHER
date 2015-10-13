BrowseSequences <- function(...) {
	.Defunct("BrowseSeqs",
		msg="'BrowseSequences' is defunct.\nUse 'BrowseSeqs' instead.\nSee help(\"DECIPHER-defunct\").")
}

BrowseSeqs <- function(myXStringSet,
	htmlFile=paste(tempdir(),"/myXStringSet.html",sep=""),
	openURL=interactive(),
	colorPatterns=TRUE,
	highlight=NA,
	patterns=c("-", alphabet(myXStringSet, baseOnly=TRUE)),
	colors=substring(rainbow(length(patterns), v=0.8, start=0.9, end=0.7), 1, 7),
	colWidth=Inf,
	...) {
	
	# error checking
	if (!is(myXStringSet, "XStringSet"))
		stop("myXStringSet must be an XStringSet.")
	if (length(myXStringSet)==0)
		stop("No sequence information to display.")
	type <- switch(class(patterns),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		0L)
	if (type > 0L) {
		patterns <- as.character(patterns)
		if (type==1L) { # DNAStringSet
			patterns <- gsub("M",
				"[ACM]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("R",
				"[AGR]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("W",
				"[ATW]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("S",
				"[CGS]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("Y",
				"[CTY]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("K",
				"[GTK]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("V",
				"[ACGMRSV]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("H",
				"[ACTMWYH]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("D",
				"[AGTRWKD]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("B",
				"[CGTSYKB]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("N",
				"[ACGTMRWSYKVHDBN]",
				patterns,
				fixed=TRUE)
		} else if (type==2L) { # RNAStringSet
			patterns <- gsub("M",
				"[ACM]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("R",
				"[AGR]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("W",
				"[AUW]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("S",
				"[CGS]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("Y",
				"[CUY]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("K",
				"[GUK]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("V",
				"[ACGMRSV]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("H",
				"[ACUMWYH]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("D",
				"[AGURWKD]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("B",
				"[CGUSYKB]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("N",
				"[ACGUMRWSYKVHDBN]",
				patterns,
				fixed=TRUE)
		} else { # AAStringSet
			patterns <- gsub("B",
				"[NDB]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("Z",
				"[QEZ]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("J",
				"[ILJ]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("X",
				"[ARNDCQEGHILKMFPSTWYVUOBJZX]",
				patterns,
				fixed=TRUE)
			patterns <- gsub("*",
				"\\*",
				patterns,
				fixed=TRUE)
		}
	} else {
		if (!is.null(patterns) && !is.character(patterns))
			stop("patterns must be a character vector.")
		if (any(grepl("=|\"|<|>|[1-9]|[a-z]", patterns)))
			stop("patterns cannot contain numbers, lower case characters, or the characters (=, <, >, \").")
	}
	if (any(patterns==""))
		stop("No patterns can be empty.")
	w <- which(patterns %in% c("?", "*", "+", "."))
	if (length(w) > 0)
		patterns[w] <- paste("\\", patterns[w], sep="")
	if (length(colors) != length(patterns) || !is.character(colors))
		stop("colors must be a character vector of the same length as patterns.")
	# check that the file exist
	if (is.character(htmlFile)) {
		htmlfile <- file(htmlFile, "w")
		on.exit(close(htmlfile))
	} else {
		stop("'htmlFile' must be a character string or connection.")
	}
	if (!is.logical(openURL) || is.na(openURL))
		stop("openURL must be TRUE or FALSE.")
	if (!is.logical(colorPatterns) && !is.numeric(colorPatterns))
		stop("colorPatterns must be a logical or numeric.")
	if (is.numeric(colorPatterns)) {
		if ((length(colorPatterns) %% 2) == 1 || length(colorPatterns) == 0)
			stop("colorPatterns must specify all start and endpoints.")
		if (any((colorPatterns[seq(2, length(colorPatterns), 2)] - colorPatterns[seq(1, length(colorPatterns), 2)]) < 0))
			stop("colorPatterns specifies a negative range.")
		if (any(colorPatterns <= 0))
			stop("colorPatterns must be a positive numeric.")
		if (length(colorPatterns) > 2)
			if (any((colorPatterns[seq(3, length(colorPatterns), 2)] - colorPatterns[seq(2, length(colorPatterns) - 2, 2)]) <= 0))
				stop("Ranges specified in colorPatterns must be non-overlapping.")
		if (max(colorPatterns) > max(width(myXStringSet)))
			stop("Ranges specified in colorPatterns are out-of-bounds.")
	}
	if (is.numeric(colorPatterns) & !is.infinite(colWidth))
		stop("colWidth must be Inf if colorPatterns is numeric.")
	if (is.null(names(myXStringSet)))
		names(myXStringSet) <- 1:length(myXStringSet)
	names(myXStringSet) <- gsub("\t", " ", names(myXStringSet), fixed=TRUE)
	if (!is.na(highlight)) {
		if (highlight < 0 || highlight > length(myXStringSet))
			stop("highlight must be 0 or the index of a sequence in myXStringSet.")
		if (highlight != floor(highlight))
			stop("highlight be be a whole number.")
	}
	# add a consensus sequence to myXStringSet
	if (is(myXStringSet, "DNAStringSet") || is(myXStringSet, "RNAStringSet") || is(myXStringSet, "AAStringSet")) {
		myXStringSet <- c(myXStringSet,
			ConsensusSequence(myXStringSet,
				...))
	} else {
		myXStringSet <- c(myXStringSet,
			as(consensusString(myXStringSet,
					...),
				class(myXStringSet)))
	}
	l <- length(myXStringSet)
	names(myXStringSet)[l] <- "Consensus"
	
	# convert the XStringSet to character strings
	html <- sapply(myXStringSet, toString)
	if (!is.na(highlight)) {
		if (highlight==0) {
			highlight <- length(html)
			index <- 1:(length(html) - 1L)
		} else {
			index <- (1:(length(html) - 1L))[-highlight]
		}
		html <- sapply(html, strsplit, split="", fixed=TRUE)
		for (i in index) {
			L <- min(length(html[[highlight]]), length(html[[i]]))
			w <- which(html[[i]][1:L]==html[[highlight]][1:L])
			if (length(w) > 0)
				html[[i]][w] <- "\u00B7"
		}
		html <- sapply(html, paste, collapse="")
	}
	
	# pad shorter sequences with tildes
	maxW <- max(width(myXStringSet))
	if (maxW==0)
		stop("No sequence information to display.")
	if (colWidth > maxW)
		colWidth <- maxW
	for (i in 1:l) {
		html[i] <- paste(html[i],
			paste(rep(" ", maxW - nchar(html[i])), collapse=""),
			sep="")
	}
	
	# create a legend that gives position
	counter <- 20
	if (counter > maxW)
		counter <- maxW
	offset <- (counter - 1) - nchar(maxW)
	if (offset < 0)
		offset <- 0
	legend <- paste(paste(rep(" ", offset), collapse=""), format(seq(counter, maxW, by=counter)), collapse="")
	counter <- ceiling(counter/2)
	tickmarks <- paste(paste(rep("'", counter - 1), collapse=""), rep("|", floor(maxW/counter)), collapse="", sep="")
	tickmarks <- paste(tickmarks, paste(rep("'", maxW - counter*floor(maxW/counter)), collapse=""), sep="")
	
	# split the html into multiple pages
	starts <- seq(1, maxW, colWidth)
	stops <- starts + colWidth - 1L
	stops[length(stops)] <- maxW
	temp <- character((length(html) + 5L)*length(starts))
	count <- 1L
	for (i in 1:length(starts)) {
		temp[count:(count + 1L)] <- substring(c(legend, tickmarks), starts[i], stops[i])
		count <- count + 2L
		temp[c(count:(count + length(html) - 2L), count + length(html))] <- substring(html, starts[i], stops[i])
		count <- count + length(html) + 3L
	}
	html <- temp
	
	# add the cumulative nucleotide lengths on the right
	if (length(starts) > 1) {
		lengths <- numeric(length(starts)*l)
		myLengths <- character(length(starts)*(l + 5L))
		for (i in 1:length(starts)) {
			lengths[((i - 1L)*l + 1L):(i*l)] <- stops[i] - starts[i] + 1L - rowSums(letterFrequency(BStringSet(substring(myXStringSet,
					starts[i],
					stops[i])),
				"-"))
			if (i > 1)
				lengths[((i - 1L)*l + 1L):(i*l)] <- lengths[((i - 1L)*l + 1L):(i*l)] + lengths[((i - 2L)*l + 1L):((i - 1L)*l)]
			myLengths[((i - 1L)*(l + 5L) + 3L):(i*(l + 5L) - 4L)] <- lengths[((i - 1L)*l + 1L):(i*l - 1L)]
			myLengths[(i*(l + 5L) - 2L)] <- lengths[i*l]
		}
	} else {
		lengths <- width(myXStringSet) - rowSums(letterFrequency(myXStringSet, "-"))
		myLengths <- c("", "", lengths[-seq(l, length(lengths), l)], "", lengths[seq(l, length(lengths), l)], "", "")
	}
	
	if (is.numeric(colorPatterns) || colorPatterns) {
		if (is.numeric(colorPatterns)) {
			htm <- substring(html, 0, colorPatterns[1] - 1)
			for (i in seq(1, length(colorPatterns), 2)) {
				htmi <- substring(html, colorPatterns[i], colorPatterns[i + 1])
				for (j in 1:length(colors)) {
					htmi <- gsub(paste("(",
							patterns[j],
							ifelse(nchar(patterns[j]) == 1, "+)", ")"),
							sep=""),
						paste("<span class=\"_",
							j,
							"\">\\1</span>",
							sep=""),
						htmi)
				}
				end <- ifelse(i==(length(colorPatterns) - 1),
					max(nchar(html)),
					colorPatterns[i + 2] - 1)
				
				htm <- paste(htm, htmi, substring(html, colorPatterns[i + 1] + 1, end), sep="")
			}
			html <- paste(htm, myLengths, sep="    ")
		} else {
			# add color tags to each nucleotide
			if (length(colors) > 0) {
				for (j in 1:length(colors)) {
					html <- gsub(paste("(",
							patterns[j],
							ifelse(nchar(patterns[j]) == 1, "+)", ")"),
							sep=""),
						paste("<span class=\"_",
							j,
							"\">\\1</span>",
							sep=""),
						html)
				}
			}
			html <- paste(html, myLengths, sep="    ")
		}
		
		# add the legend and consensus sequence
		html <- paste(format(c("", "", names(myXStringSet)[1:(l - 1)], "", "Consensus", "", ""), justify="right"),
				html,
			sep="    ")
		
		styles <- character()
		if (length(colors) > 0) {
			for (i in 1:length(colors)) {
				styles <- paste(styles,
					"span._", i,
					" {background:", colors[i],
					"; color: #FFF;} ",
					sep="")
			}
		}
		styles <- paste("<style type=text/css> ",
			styles,
			"</style>",
			sep="")
		html <- c('<html><head><meta http-equiv="Content-Type" content="text/html; charset=utf-8"></head>',styles,"<pre>", html, "</pre></html>")
	} else {
		# add the legend and consensus sequence
		html <- paste(format(c("", "", names(myXStringSet)[1:(l - 1)], "", "Consensus", "", ""), justify="right"),	
			html,
			sep="    ")
		
		html <- c("<html>","<pre>", html, "</pre></html>")
	}
	
	# replace unicode 'middle dot' with html entity
	html <- gsub("\u00B7", "&#183;", html, fixed=TRUE)
	
	writeLines(html, htmlfile)
	
	if (openURL)
		browseURL(path.expand(htmlFile))
	
	invisible(htmlFile)
}
