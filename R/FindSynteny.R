FindSynteny <- function(dbFile,
	tblName="Seqs",
	identifier="",
	useFrames=FALSE,
	alphabet=c("MF", "ILV", "A", "C", "WYQHP", "G", "TSN", "RK", "DE"),
	geneticCode=GENETIC_CODE,
	sepCost=-0.01,
	gapCost=-0.2,
	shiftCost=-20,
	codingCost=-3,
	maxSep=5000,
	maxGap=5000,
	minScore=200,
	dropScore=-100,
	maskRepeats=TRUE,
	storage=0.5,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.logical(useFrames))
		stop("useFrames must be a logical.")
	if (!is.numeric(sepCost))
		stop("sepCost must be a single numeric.")
	if (sepCost > 0)
		stop("sepCost must be less than or equal to zero.")
	if (!is.numeric(gapCost))
		stop("gapCost must be a single numeric.")
	if (gapCost > 0)
		stop("gapCost must be less than or equal to zero.")
	if (!is.numeric(shiftCost))
		stop("shiftCost must be a single numeric.")
	if (shiftCost > 0)
		stop("shiftCost must be less than or equal to zero.")
	if (!is.numeric(codingCost))
		stop("codingCost must be a single numeric.")
	if (codingCost > 0)
		stop("codingCost must be less than or equal to zero.")
	if (!is.numeric(maxSep))
		stop("maxSep must be a single numeric.")
	if (maxSep <= 0)
		stop("maxSep must be greater than zero.")
	if (!is.numeric(maxGap))
		stop("maxGap must be a single numeric.")
	if (maxGap < 0)
		stop("maxGap must be at least zero.")
	if (!is.numeric(dropScore))
		stop("dropScore must be a single numeric.")
	if (dropScore > 0)
		stop("dropScore must be less than or equal to zero.")
	if (!is.numeric(minScore))
		stop("minScore must be a single numeric.")
	if (minScore <= 0)
		stop("minScore must be greater than zero.")
	if (!is.logical(maskRepeats))
		stop("maskRepeats must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(storage))
		stop("storage must be a numeric.")
	if (storage < 0)
		stop("storage must be at least zero.")
	storage <- storage*1e9 # convert to bytes
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
	
	ids <- dbGetQuery(dbConn,
		paste("select distinct identifier from",
			tblName))$identifier
	if (identifier[1]=="") {
		identifier <- ids
	} else {
		w <- which(!(identifier %in% ids))
		if (length(w) > 0)
			stop("identifier(s) not found in ",
				tblName,
				": ",
				paste(identifier[w],
					collapse=", "))
	}
	l <- length(identifier)
	if (l < 2)
		stop("At least two identifers are required in ",
			tblName,
			".")
	
	codons <- c('AAA', 'AAC', 'AAG', 'AAT',
		'ACA', 'ACC', 'ACG', 'ACT',
		'AGA', 'AGC', 'AGG', 'AGT',
		'ATA', 'ATC', 'ATG', 'ATT',
		'CAA', 'CAC', 'CAG', 'CAT',
		'CCA', 'CCC', 'CCG', 'CCT',
		'CGA', 'CGC', 'CGG', 'CGT',
		'CTA', 'CTC', 'CTG', 'CTT',
		'GAA', 'GAC', 'GAG', 'GAT',
		'GCA', 'GCC', 'GCG', 'GCT',
		'GGA', 'GGC', 'GGG', 'GGT',
		'GTA', 'GTC', 'GTG', 'GTT',
		'TAA', 'TAC', 'TAG', 'TAT',
		'TCA', 'TCC', 'TCG', 'TCT',
		'TGA', 'TGC', 'TGG', 'TGT',
		'TTA', 'TTC', 'TTG', 'TTT')
	f <- function(gC) {
		m <- match(codons, names(gC))
		if (any(is.na(m)))
			stop("geneticCode must contain all 64 codons.")
		AAStringSet(paste(gC[m],
			collapse=""))
	}
	if (is.list(geneticCode)) {
		if (length(geneticCode)!=l)
			stop("The list geneticCode must have as many items as the number of identifiers.")
		if (!is.null(names(geneticCode))) {
			m <- match(identifier, names(geneticCode))
			if (any(is.na(m)))
				stop("All identifiers must be present in the names of the list geneticCode.")
			geneticCode <- geneticCode[m]
			geneticCode <- lapply(geneticCode, f)
		}
	} else { # all identifiers use the same geneticCode
		geneticCode <- rep(list(f(geneticCode)), l)
	}
	
	if (any(alphabet==""))
		stop("No elements of alphabet can be empty.")
	r <- strsplit(alphabet, "", fixed=TRUE)
	alphabet <- setNames(rep(0L, 20),
		AA_STANDARD)
	for (i in seq_along(r)) {
		w <- which(!(r[[i]] %in% AA_STANDARD))
		if (length(w) > 0)
			stop("Unrecognized letter(s) found in alphabet:  ",
				paste(r[[i]][w], collapse=", "),
				".")
		w <- which(alphabet[r[[i]]] != 0L)
		if (length(w) > 0)
			stop("Repeated amino acids found in alphabet:  ",
				paste(r[[i]][w], collapse=", "),
				".")
		alphabet[r[[i]]] <- i
	}
	w <- which(alphabet==0L)
	if (length(w) > 0)
		stop("Standard amino acids missing from alphabet:  ",
			paste(names(w), collapse=", "),
			".")
	size <- max(alphabet)
	if (size==1)
		stop("More than one grouping of amino acids is required in the alphabet.")
	n <- as.integer(floor(log(4294967295, size)))
	alphabet <- alphabet - 1L
	
	results <- matrix(data.frame(),
		nrow=l,
		ncol=l,
		dimnames=list(identifier, identifier))
	empty_upper <- matrix(integer(),
		nrow=0,
		ncol=8,
		dimnames=list(NULL,
			c("index1",
				"index2",
				"strand",
				"width",
				"start1",
				"start2",
				"frame1",
				"frame2")))
	empty_lower <- matrix(integer(),
		nrow=0,
		ncol=10,
		dimnames=list(NULL,
			c("index1",
				"index2",
				"strand",
				"score",
				"start1",
				"start2",
				"end1",
				"end2",
				"first_hit",
				"last_hit")))
	
	store <- rep(list(list(S=list(),
			E=list(`nt`=list(),
				`nt_rc`=list(),
				`aa`=list(list(), list(), list()),
				`aa_rc`=list(list(), list(), list())),
			O=list(`nt`=list(),
				`nt_rc`=list(),
				`aa`=list(list(), list(), list()),
				`aa_rc`=list(list(), list(), list())))),
		l)
	
	if (verbose) {
		pBar <- txtProgressBar(style=3)
		tot <- l*(l - 1)/2
		its <- 0
	}
	for (g1 in 1:(l - 1)) {
		# remove unnecessary items from store
		store[g1][[1L]][["E"]][["nt_rc"]] <- list()
		store[g1][[1L]][["E"]][["aa_rc"]] <- list()
		store[g1][[1L]][["O"]][["nt_rc"]] <- list()
		store[g1][[1L]][["O"]][["aa_rc"]] <- list()
		
		s1 <- store[g1][[1L]][["S"]]
		if (length(s1)==0) {
			s1 <- SearchDB(dbConn,
				tblName=tblName,
				identifier=identifier[g1],
				type="DNAStringSet",
				removeGaps="all",
				processors=processors,
				verbose=FALSE)
		} else {
			s1 <- s1[[1L]]
			store[g1][[1L]][["S"]] <- list() # clear s1 from store
		}
		w1 <- width(s1)
		INDEX1 <- seq_along(s1)
		results[g1, g1] <- list(w1)
		
		w <- which(w1 <= 32)
		if (length(w) > 0) {
			s1 <- s1[-w]
			INDEX1 <- INDEX1[-w]
		}
		if (length(s1)==0) {
			for (g2 in (g1 + 1):l) {
				results[g1, g2][[1]] <- empty_upper
				results[g2, g1][[1]] <- empty_lower
			}
			
			next
		}
		
		WIDTH1 <- cumsum(width(s1))	
		if (length(s1) > 0) {
			seq1 <- DNAStringSet(unlist(s1))
		} else {
			seq1 <- s1
		}
		
		for (g2 in (g1 + 1):l) {
			s2 <- store[g2][[1L]][["S"]]
			if (length(s2)==0) {
				s2 <- SearchDB(dbConn,
					tblName=tblName,
					identifier=identifier[g2],
					type="DNAStringSet",
					removeGaps="all",
					processors=processors,
					verbose=FALSE)
				if ((object.size(s2) + object.size(store)) < storage)
					store[g2][[1L]][["S"]][[1L]] <- s2
			} else {
				s2 <- s2[[1L]]
			}
			w2 <- width(s2)
			INDEX2 <- seq_along(s2)
			
			w <- which(w2 <= 32)
			if (length(w) > 0) {
				s2 <- s2[-w]
				INDEX2 <- INDEX2[-w]
			}
			if (length(s2)==0) {
				results[g1, g2][[1]] <- empty_upper
				results[g2, g1][[1]] <- empty_lower
				
				next
			}
			
			WIDTH2 <- cumsum(width(s2))
			if (length(s2) > 0) {
				seq2 <- DNAStringSet(unlist(s2))
			} else {
				seq2 <- s2
			}
			
			# calculate the k-mer size with approximately
			# a 1% probability of occurring by chance
			M <- max(width(s1), width(s2))*100
			N <- as.integer(ceiling(log(M, 4)))
			if (N > 16L)
				N <- 16L
			N_AA <- as.integer(ceiling(log(M/3, size)))
			if (N_AA > n)
				N_AA <- n
			
			E1 <- store[g1][[1L]][["E"]][["nt"]][N][[1L]]
			if (is.null(E1)) {
				E1 <- .Call("enumerateSequence",
					seq1,
					N,
					PACKAGE="DECIPHER")[[1]]
				for (i in which(WIDTH1 > (N - 2) & WIDTH1 < length(E1)))
					E1[(WIDTH1[i] - (N - 2)):WIDTH1[i]] <- NA
				if (maskRepeats)
					.Call("maskRepeats",
						E1,
						N,
						7L, # minimum period
						12L, # maximum period
						30L, # minimum length
						PACKAGE="DECIPHER")
				if ((object.size(E1) + object.size(store)) < storage)
					store[g1][[1L]][["E"]][["nt"]][[N]] <- E1
			}
			
			e2 <- store[g2][[1L]][["E"]][["nt"]][N][[1L]]
			if (is.null(e2)) {
				e2 <- .Call("enumerateSequence",
					seq2,
					N,
					PACKAGE="DECIPHER")[[1]]
				for (i in which(WIDTH2 > (N - 2) & WIDTH2 < length(e2)))
					e2[(WIDTH2[i] - (N - 2)):WIDTH2[i]] <- NA
				if (maskRepeats)
					.Call("maskRepeats",
						e2,
						N,
						7L, # minimum period
						12L, # maximum period
						30L, # minimum length
						PACKAGE="DECIPHER")
				if ((object.size(e2) + object.size(store)) < storage)
					store[g2][[1L]][["E"]][["nt"]][[N]] <- e2
			}
			
			O1 <- store[g1][[1L]][["O"]][["nt"]][N][[1L]]
			if (is.null(O1)) {
				O1 <- .Call("radixOrder",
					E1,
					0L,
					PACKAGE="DECIPHER")
				if ((object.size(O1) + object.size(store)) < storage)
					store[g1][[1L]][["O"]][["nt"]][[N]] <- O1
			}
			
			o2 <- store[g2][[1L]][["O"]][["nt"]][N][[1L]]
			if (is.null(o2)) {
				o2 <- .Call("radixOrder",
					e2,
					0L,
					PACKAGE="DECIPHER")
				if ((object.size(o2) + object.size(store)) < storage)
					store[g2][[1L]][["O"]][["nt"]][[N]] <- o2
			}
			
			# match E1 to e2
			m <- .Call("intMatchOnce",
				E1,
				e2,
				O1,
				o2,
				PACKAGE="DECIPHER")
			.Call("fillOverlaps", m, N, PACKAGE="DECIPHER")
			r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
			w <- which(r@values==1)
			widths <- r@lengths[w] + N
			ends <- cumsum(r@lengths)[w] + N
			
			ends1 <- ends
			starts1 <- ends1 - widths + 1L
			ends2 <- m[ends - N] + N
			starts2 <- ends2 - widths + 1L
			
			# re-index the sequences by contig
			index1 <- .Call("indexByContig",
				starts1,
				ends1,
				seq_along(starts1),
				INDEX1,
				WIDTH1)
			o <- .Call("radixOrder",
				starts2,
				1L,
				PACKAGE="DECIPHER")
			index2 <- .Call("indexByContig",
				starts2,
				ends2,
				o,
				INDEX2,
				WIDTH2)
			
			if (useFrames) {
				x.s <- list(starts1)
				x.e <- list(ends1)
				x.i <- list(index1)
				y.s <- list(starts2)
				y.e <- list(ends2)
				y.i <- list(index2)
				x.f <- y.f <- list(rep(0L, length(starts1)))
				weights <- list(widths)
				
				for (rF1 in 1:3) {
					t1 <- .Call("basicTranslate",
						s1,
						geneticCode[[g1]],
						rep(rF1, length(s1)),
						PACKAGE="DECIPHER")
					width1 <- cumsum(width(t1))
					if (length(t1) > 0) {
						t1 <- AAStringSet(unlist(t1))
					}
					
					e1 <- store[g1][[1L]][["E"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(e1)) {
						e1 <- .Call("enumerateSequenceReducedAA",
							t1,
							N_AA,
							alphabet,
							PACKAGE="DECIPHER")[[1]]
						for (i in which(width1 > (N_AA - 2) & width1 < length(e1)))
							e1[(width1[i] - (N_AA - 2)):width1[i]] <- NA
						if (maskRepeats)
							.Call("maskRepeats",
								e1,
								N_AA,
								3L, # minimum period
								11L, # maximum period
								15L, # minimum length
								PACKAGE="DECIPHER")
						if ((object.size(e1) + object.size(store)) < storage)
							store[g1][[1L]][["E"]][["aa"]][[rF1]][[N_AA]] <- e1
					}
					width1 <- width1*3L
					
					o1 <- store[g1][[1L]][["O"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(o1)) {
						o1 <- .Call("radixOrder",
							e1,
							0L,
							PACKAGE="DECIPHER")
						if ((object.size(o1) + object.size(store)) < storage)
							store[g1][[1L]][["O"]][["aa"]][[rF1]][[N_AA]] <- o1
					}
					
					for (rF2 in 1:3) {
						t2 <- .Call("basicTranslate",
							s2,
							geneticCode[[g2]],
							rep(rF2, length(s2)),
							PACKAGE="DECIPHER")
						width2 <- cumsum(width(t2))
						if (length(t2) > 0) {
							t2 <- AAStringSet(unlist(t2))
						}
						
						e2 <- store[g2][[1L]][["E"]][["aa"]][[rF2]][N_AA][[1L]]
						if (is.null(e2)) {
							e2 <- .Call("enumerateSequenceReducedAA",
								t2,
								N_AA,
								alphabet,
								PACKAGE="DECIPHER")[[1]]
							for (i in which(width2 > (N_AA - 2) & width2 < length(e2)))
								e2[(width2[i] - (N_AA - 2)):width2[i]] <- NA
							if (maskRepeats)
								.Call("maskRepeats",
									e2,
									N_AA,
									3L, # minimum period
									11L, # maximum period
									15L, # minimum length
									PACKAGE="DECIPHER")
							if ((object.size(e2) + object.size(store)) < storage)
								store[g2][[1L]][["E"]][["aa"]][[rF2]][[N_AA]] <- e2
						}
						width2 <- width2*3L
						
						o2 <- store[g2][[1L]][["O"]][["aa"]][[rF2]][N_AA][[1L]]
						if (is.null(o2)) {
							o2 <- .Call("radixOrder",
								e2,
								0L,
								PACKAGE="DECIPHER")
							if ((object.size(o2) + object.size(store)) < storage)
								store[g2][[1L]][["O"]][["aa"]][[rF2]][[N_AA]] <- o2
						}
						
						# match e1 to e2
						m <- .Call("intMatchOnce",
							e1,
							e2,
							o1,
							o2,
							PACKAGE="DECIPHER")
						.Call("fillOverlaps", m, N_AA, PACKAGE="DECIPHER")
						r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
						w <- which(r@values==1)
						widths <- r@lengths[w] + N_AA
						ends <- cumsum(r@lengths)[w] + N_AA
						
						ends1 <- ends
						starts1 <- ends1 - widths + 1L
						ends2 <- m[ends - N_AA] + N_AA
						starts2 <- ends2 - widths + 1L
						
						# convert from amino acid to nucleotide positions
						starts1 <- (starts1 - 1L)*3L + rF1
						ends1 <- ends1*3L + rF1 - 1L
						starts2 <- (starts2 - 1L)*3L + rF2
						ends2 <- ends2*3L + rF2 - 1L
						
						# re-index the sequences by contig
						index1 <- .Call("indexByContig",
							starts1,
							ends1,
							seq_along(starts1),
							INDEX1,
							width1)
						o <- .Call("radixOrder",
							starts2,
							1L,
							PACKAGE="DECIPHER")
						index2 <- .Call("indexByContig",
							starts2,
							ends2,
							o,
							INDEX2,
							width2)
						
						x.s <- c(x.s, list(starts1))
						x.e <- c(x.e, list(ends1))
						x.i <- c(x.i, list(index1))
						y.s <- c(y.s, list(starts2))
						y.e <- c(y.e, list(ends2))
						y.i <- c(y.i, list(index2))
						x.f <- c(x.f, list(rep(rF1, length(starts1))))
						y.f <- c(y.f, list(rep(rF2, length(starts1))))
						weights <- c(weights, list(widths*3L))
					}
				}
				
				x.s <- unlist(x.s)
				x.e <- unlist(x.e)
				x.i <- unlist(x.i)
				y.s <- unlist(y.s)
				y.e <- unlist(y.e)
				y.i <- unlist(y.i)
				x.f <- unlist(x.f)
				y.f <- unlist(y.f)
				weights <- unlist(weights)
			} else {
				x.s <- starts1
				x.e <- ends1
				x.i <- index1
				y.s <- starts2
				y.e <- ends2
				y.i <- index2
				x.f <- y.f <- rep(0L, length(starts1))
				weights <- widths
			}
			
			# order by increasing sequence index in g1
			o <- order(x.i, x.s)
			x.s <- x.s[o]
			x.e <- x.e[o]
			x.i <- x.i[o]
			y.s <- y.s[o]
			y.e <- y.e[o]
			y.i <- y.i[o]
			x.f <- x.f[o]
			y.f <- y.f[o]
			weights <- weights[o]
			
			chains <- .Call("chainSegments",
				x.s,
				x.e,
				x.i,
				x.f,
				y.s,
				y.e,
				y.i,
				y.f,
				as.double(weights*ifelse(x.f==0L,
					log(4),
					log(size)/3)),
				sepCost,
				gapCost,
				shiftCost,
				codingCost,
				maxSep,
				maxGap,
				order(x.i, x.e) - 1L,
				minScore,
				PACKAGE="DECIPHER")
			
			scores <- chains[[2]]
			chains <- chains[[1]]
			used <- unlist(chains,
				use.names=FALSE)
			count <- 0L
			for (i in seq_along(chains)) {
				chains[[i]] <- (count + 1L):(count + length(chains[[i]]))
				count <- count + length(chains[[i]])
			}
			results[g2, g1][[1]] <- list(chain=chains,
				score=scores)
			
			result1 <- matrix(c(x.i[used],
					y.i[used],
					rep(0L, length(used)),
					x.e[used] - x.s[used] + 1L,
					x.s[used],
					y.s[used],
					x.f[used],
					y.f[used]),
				nrow=length(used),
				ncol=8,
				dimnames=list(NULL,
					c("index1",
						"index2",
						"strand",
						"width",
						"start1",
						"start2",
						"frame1",
						"frame2")))
			
			if (verbose) {
				its <- its + 0.5
				setTxtProgressBar(pBar, its/tot)
			}
			
			s3 <- reverseComplement(s2)
			if (length(s3) > 0) {
				seq2 <- DNAStringSet(unlist(s3))
			} else {
				seq2 <- s3
			}
			
			e2 <- store[g2][[1L]][["E"]][["nt_rc"]][N][[1L]]
			if (is.null(e2)) {
				e2 <- .Call("enumerateSequence",
					seq2,
					N,
					PACKAGE="DECIPHER")[[1]]
				for (i in which(WIDTH2 > (N - 2) & WIDTH2 < length(e2)))
					e2[(WIDTH2[i] - (N - 2)):WIDTH2[i]] <- NA
				if (maskRepeats)
					.Call("maskRepeats",
						e2,
						N,
						7L, # minimum period
						12L, # maximum period
						30L, # minimum length
						PACKAGE="DECIPHER")
				if ((object.size(e2) + object.size(store)) < storage)
					store[g2][[1L]][["E"]][["nt_rc"]][[N]] <- e2
			}
			
			o2 <- store[g2][[1L]][["O"]][["nt_rc"]][N][[1L]]
			if (is.null(o2)) {
				o2 <- .Call("radixOrder",
					e2,
					0L,
					PACKAGE="DECIPHER")
				if ((object.size(o2) + object.size(store)) < storage)
					store[g2][[1L]][["O"]][["nt_rc"]][[N]] <- o2
			}
			
			# match E1 to e2
			m <- .Call("intMatchOnce",
				E1,
				e2,
				O1,
				o2,
				PACKAGE="DECIPHER")
			.Call("fillOverlaps", m, N, PACKAGE="DECIPHER")
			r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
			w <- which(r@values==1)
			widths <- r@lengths[w] + N
			ends <- cumsum(r@lengths)[w] + N
			
			ends1 <- ends
			starts1 <- ends1 - widths + 1L
			ends2 <- m[ends - N] + N
			starts2 <- ends2 - widths + 1L
			
			# re-index the sequences by contig
			index1 <- .Call("indexByContig",
				starts1,
				ends1,
				seq_along(starts1),
				INDEX1,
				WIDTH1)
			o <- .Call("radixOrder",
				starts2,
				1L,
				PACKAGE="DECIPHER")
			index2 <- .Call("indexByContig",
				starts2,
				ends2,
				o,
				INDEX2,
				WIDTH2)
			
			# correct for reverse complement positioning
			w <- w2[index2]
			temp <- w - starts2 + 1L
			starts2 <- w - ends2 + 1L
			ends2 <- temp
			
			if (useFrames) {
				x.s <- list(starts1)
				x.e <- list(ends1)
				x.i <- list(index1)
				y.s <- list(starts2)
				y.e <- list(ends2)
				y.i <- list(index2)
				x.f <- y.f <- list(rep(0L, length(starts1)))
				weights <- list(widths)
				
				for (rF1 in 1:3) {
					t1 <- .Call("basicTranslate",
						s1,
						geneticCode[[g1]],
						rep(rF1, length(s1)),
						PACKAGE="DECIPHER")
					width1 <- cumsum(width(t1))
					if (length(t1) > 0) {
						t1 <- AAStringSet(unlist(t1))
					}
					
					e1 <- store[g1][[1L]][["E"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(e1)) {
						e1 <- .Call("enumerateSequenceReducedAA",
							t1,
							N_AA,
							alphabet,
							PACKAGE="DECIPHER")[[1]]
						for (i in which(width1 > (N_AA - 2) & width1 < length(e1)))
							e1[(width1[i] - (N_AA - 2)):width1[i]] <- NA
						if (maskRepeats)
							.Call("maskRepeats",
								e1,
								N_AA,
								3L, # minimum period
								11L, # maximum period
								15L, # minimum length
								PACKAGE="DECIPHER")
						if ((object.size(e1) + object.size(store)) < storage)
							store[g1][[1L]][["E"]][["aa"]][[rF1]][[N_AA]] <- e1
					}
					width1 <- width1*3L
					
					o1 <- store[g1][[1L]][["O"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(o1)) {
						o1 <- .Call("radixOrder",
							e1,
							0L,
							PACKAGE="DECIPHER")
						if ((object.size(o1) + object.size(store)) < storage)
							store[g1][[1L]][["O"]][["aa"]][[rF1]][[N_AA]] <- o1
					}
					
					for (rF2 in 1:3) {
						t2 <- .Call("basicTranslate",
							s3,
							geneticCode[[g2]],
							rep(rF2, length(s3)),
							PACKAGE="DECIPHER")
						width2 <- cumsum(width(t2))
						if (length(t2) > 0) {
							t2 <- AAStringSet(unlist(t2))
						}
						
						e2 <- store[g2][[1L]][["E"]][["aa_rc"]][[rF2]][N_AA][[1L]]
						if (is.null(e2)) {
							e2 <- .Call("enumerateSequenceReducedAA",
								t2,
								N_AA,
								alphabet,
								PACKAGE="DECIPHER")[[1]]
							for (i in which(width2 > (N_AA - 2) & width2 < length(e2)))
								e2[(width2[i] - (N_AA - 2)):width2[i]] <- NA
							if (maskRepeats)
								.Call("maskRepeats",
									e2,
									N_AA,
									3L, # minimum period
									11L, # maximum period
									15L, # minimum length
									PACKAGE="DECIPHER")
							if ((object.size(e2) + object.size(store)) < storage)
								store[g2][[1L]][["E"]][["aa_rc"]][[rF2]][[N_AA]] <- e2
						}
						width2 <- width2*3L
						
						o2 <- store[g2][[1L]][["O"]][["aa_rc"]][[rF2]][N_AA][[1L]]
						if (is.null(o2)) {
							o2 <- .Call("radixOrder",
								e2,
								0L,
								PACKAGE="DECIPHER")
							if ((object.size(o2) + object.size(store)) < storage)
								store[g2][[1L]][["O"]][["aa_rc"]][[rF2]][[N_AA]] <- o2
						}
						
						# match e1 to e2
						m <- .Call("intMatchOnce",
							e1,
							e2,
							o1,
							o2,
							PACKAGE="DECIPHER")
						.Call("fillOverlaps", m, N_AA, PACKAGE="DECIPHER")
						r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
						w <- which(r@values==1)
						widths <- r@lengths[w] + N_AA
						ends <- cumsum(r@lengths)[w] + N_AA
						
						ends1 <- ends
						starts1 <- ends1 - widths + 1L
						ends2 <- m[ends - N_AA] + N_AA
						starts2 <- ends2 - widths + 1L
						
						# convert from amino acid to nucleotide positions
						starts1 <- (starts1 - 1L)*3L + rF1
						ends1 <- ends1*3L + rF1 - 1L
						starts2 <- (starts2 - 1L)*3L + rF2
						ends2 <- ends2*3L + rF2 - 1L
						
						# re-index the sequences by contig
						index1 <- .Call("indexByContig",
							starts1,
							ends1,
							seq_along(starts1),
							INDEX1,
							width1)
						o <- .Call("radixOrder",
							starts2,
							1L,
							PACKAGE="DECIPHER")
						index2 <- .Call("indexByContig",
							starts2,
							ends2,
							o,
							INDEX2,
							width2)
						
						# correct for reverse complement positioning
						w <- w2[index2]
						temp <- w - starts2 + 1L
						starts2 <- w - ends2 + 1L
						ends2 <- temp
						
						x.s <- c(x.s, list(starts1))
						x.e <- c(x.e, list(ends1))
						x.i <- c(x.i, list(index1))
						y.s <- c(y.s, list(starts2))
						y.e <- c(y.e, list(ends2))
						y.i <- c(y.i, list(index2))
						x.f <- c(x.f, list(rep(rF1, length(starts1))))
						y.f <- c(y.f, list(rep(rF2, length(starts1))))
						weights <- c(weights, list(widths*3L))
					}
				}
				
				x.s <- unlist(x.s)
				x.e <- unlist(x.e)
				x.i <- unlist(x.i)
				y.s <- unlist(y.s)
				y.e <- unlist(y.e)
				y.i <- unlist(y.i)
				x.f <- unlist(x.f)
				y.f <- unlist(y.f)
				weights <- unlist(weights)
			} else {
				x.s <- starts1
				x.e <- ends1
				x.i <- index1
				y.s <- starts2
				y.e <- ends2
				y.i <- index2
				x.f <- y.f <- rep(0L, length(starts1))
				weights <- widths
			}
			
			# order by increasing sequence index in g1
			o <- order(x.i, x.s)
			x.s <- x.s[o]
			x.e <- x.e[o]
			x.i <- x.i[o]
			y.s <- y.s[o]
			y.e <- y.e[o]
			y.i <- y.i[o]
			x.f <- x.f[o]
			y.f <- y.f[o]
			weights <- weights[o]
			
			if (length(y.s) > 0) {
				d <- max(y.e) - (y.e + y.s)
			} else {
				d <- integer()
			}
			chains <- .Call("chainSegments",
				x.s,
				x.e,
				x.i,
				x.f,
				y.s + d,
				y.e + d,
				y.i,
				y.f,
				as.double(weights*ifelse(x.f==0L,
					log(4),
					log(size)/3)),
				sepCost,
				gapCost,
				shiftCost,
				codingCost,
				maxSep,
				maxGap,
				order(x.i, x.e) - 1L,
				minScore,
				PACKAGE="DECIPHER")
			
			scores <- chains[[2]]
			chains <- chains[[1]]
			used <- unlist(chains,
				use.names=FALSE)
			for (i in seq_along(chains)) {
				chains[[i]] <- (count + 1L):(count + length(chains[[i]]))
				count <- count + length(chains[[i]])
			}
			results[g2, g1][[1]] <- list(chain=c(results[g2, g1][[1]]$chain,
					chains),
				score=c(results[g2, g1][[1]]$score,
					scores))
			
			result2 <- matrix(c(x.i[used],
					y.i[used],
					rep(1L, length(used)),
					x.e[used] - x.s[used] + 1L,
					x.s[used],
					y.e[used],
					x.f[used],
					y.f[used]),
				nrow=length(used),
				ncol=8,
				dimnames=list(NULL,
					c("index1",
						"index2",
						"strand",
						"width",
						"start1",
						"start2",
						"frame1",
						"frame2")))
			
			result <- rbind(result1, result2)
			results[g1, g2][[1]] <- result
			
			# determine ranges
			chains <- results[g2, g1][[1]]$chain
			starts <- unlist(lapply(chains, `[`, 1L))
			ends <- unlist(lapply(chains,
				function(x) x[length(x)]))
			starts1 <- result[starts, "start1"]
			ends1 <- result[ends, "start1"] + result[ends, "width"] - 1L
			strand <- result[starts, "strand"]
			starts2 <- ifelse(strand==0,
				result[starts, "start2"],
				result[ends, "start2"] + 1L - result[ends, "width"])
			ends2 <- ifelse(strand==0,
				result[ends, "start2"] + result[ends, "width"] - 1L,
				result[starts, "start2"])
			index1 <- result[starts, "index1"]
			index2 <- result[starts, "index2"]
			score <- results[g2, g1][[1]]$score
			results[g2, g1][[1]] <- matrix(c(index1,
					index2,
					strand,
					as.integer(score),
					starts1,
					starts2,
					ends1,
					ends2,
					starts,
					ends),
					ncol=10,
					dimnames=list(NULL,
						c("index1",
							"index2",
							"strand",
							"score",
							"start1",
							"start2",
							"end1",
							"end2",
							"first_hit",
							"last_hit")))
			
			# remove overlap
			nblocks <- length(strand)
			if (nblocks > 1) {
				remove <- logical(nblocks)
				hits <- list()
				
				ss1 <- result[, "start1"]
				ee1 <- result[, "start1"] + result[, "width"] - 1L
				ss2 <- ifelse(result[, "strand"]==0,
					result[, "start2"],
					result[, "start2"] + 1L - result[, "width"])
				ee2 <- ifelse(result[, "strand"]==0,
					result[, "start2"] + result[, "width"] - 1L,
					result[, "start2"])
				
				o <- order(index1,
					index2,
					starts1,
					-ends1)
				last <- 1L
				for (i in 2:nblocks) {
					if (starts1[o[i]] <= ends1[o[last]] &&
						index1[o[i]]==index1[o[last]] &&
						index2[o[i]]==index2[o[last]]) {
						if (ends1[o[i]] <= ends1[o[last]]) {
							# completely overlapping
							if (score[o[i]] >= score[o[last]]) {
								remove[o[last]] <- TRUE
								last <- i
							} else {
								remove[o[i]] <- TRUE
							}
						} else {
							# remove partial overlap from last
							chain_l <- chains[[o[last]]]
							over_l <- which(ss1[chain_l] >= (starts1[o[i]] - 2)) # min anchor width of 2
							if (length(over_l)==length(chain_l)) {
								remove[o[last]] <- TRUE
								last <- i
								next
							} else if (length(over_l) > 0) {
								# remove completely overlapping hits
								hits[[length(hits) + 1L]] <- chain_l[over_l]
								chain_l <- chain_l[-over_l]
								chains[[o[last]]] <- chain_l
							}
							
							cl <- chain_l[length(chain_l)]
							overlap <- ee1[cl] - starts1[o[i]] + 1L
							if (overlap > 0) {
								# remove overlap from last hit
								results[g1, g2][[1]][cl, "width"] <- results[g1, g2][[1]][cl, "width"] - overlap
								ee1[cl] <- ee1[cl] - overlap
								if (strand[o[last]]==0) {
									ee2[cl] <- ee2[cl] - overlap
								} else {
									ss2[cl] <- ss2[cl] + overlap
								}
							}
							ends1[o[last]] <- ee1[cl]
							results[g2, g1][[1]][o[last], "end1"] <- ee1[cl]
							if (strand[o[last]]==0) {
								ends2[o[last]] <- ee2[cl]
								results[g2, g1][[1]][o[last], "end2"] <- ee2[cl]
							} else {
								starts2[o[last]] <- ss2[cl]
								results[g2, g1][[1]][o[last], "start2"] <- ss2[cl]
							}
							last <- i
						}
					} else { # non-overlapping
						last <- i
					}
				}
				
				o <- order(index2,
					index1,
					starts2,
					-ends2)
				w <- which(!remove[o])
				last <- w[1]
				for (i in w[-1]) {
					if (starts2[o[i]] <= ends2[o[last]] &&
						index1[o[i]]==index1[o[last]] &&
						index2[o[i]]==index2[o[last]]) {
						if (ends2[o[i]] <= ends2[o[last]]) {
							# completely overlapping
							if (score[o[i]] >= score[o[last]]) {
								remove[o[last]] <- TRUE
								last <- i
							} else {
								remove[o[i]] <- TRUE
							}
						} else {
							# remove partial overlap from last
							chain_l <- chains[[o[last]]]
							over_l <- which(ss2[chain_l] >= (starts2[o[i]] - 2)) # min anchor width of 2
							if (length(over_l)==length(chain_l)) {
								remove[o[last]] <- TRUE
								last <- i
								next
							} else if (length(over_l) > 0) {
								# remove completely overlapping hits
								hits[[length(hits) + 1L]] <- chain_l[over_l]
								chain_l <- chain_l[-over_l]
								chains[[o[last]]] <- chain_l
							}
							if (strand[o[last]]==0) {
								cl <- chain_l[length(chain_l)]
							} else {
								cl <- chain_l[1]
							}
							overlap <- ee2[cl] - starts2[o[i]] + 1L
							if (overlap > 0) {
								# remove overlap from last hit
								results[g1, g2][[1]][cl, "width"] <- results[g1, g2][[1]][cl, "width"] - overlap
								if (strand[o[last]]==0) {
									ee1[cl] <- ee1[cl] - overlap
									ee2[cl] <- ee2[cl] - overlap
								} else {
									ss1[cl] <- ss1[cl] + overlap
									ee2[cl] <- ee2[cl] - overlap
									results[g1, g2][[1]][cl, "start1"] <- results[g1, g2][[1]][cl, "start1"] + overlap
									results[g1, g2][[1]][cl, "start2"] <- results[g1, g2][[1]][cl, "start2"] - overlap
								}
							}
							if (strand[o[last]]==0) {
								results[g2, g1][[1]][o[last], "end1"] <- ee1[cl]
								results[g2, g1][[1]][o[last], "end2"] <- ee2[cl]
							} else {
								results[g2, g1][[1]][o[last], "start1"] <- ss1[cl]
								results[g2, g1][[1]][o[last], "end2"] <- ee2[cl]
							}
							last <- i
						}
					} else { # non-overlapping
						last <- i
					}
				}
				
				widths <- results[g2, g1][[1]][, "end1"] - results[g2, g1][[1]][, "start1"] + 1L
				remove <- which(remove |
					widths*log(4) < minScore) # max possible score is too low
				if (length(remove) > 0 || length(hits) > 0) {
					if (length(remove) > 0)
						results[g2, g1][[1]] <- results[g2, g1][[1]][-remove,, drop=FALSE]
					hits <- c(unlist(hits), unlist(chains[remove]))
					results[g1, g2][[1]] <- results[g1, g2][[1]][-hits,, drop=FALSE]
					if (length(remove) > 0)
						chains <- chains[-remove]
					count <- 0L
					for (i in seq_along(chains)) {
						results[g2, g1][[1]][i, "first_hit"] <- count + 1L
						count <- count + length(chains[[i]])
						results[g2, g1][[1]][i, "last_hit"] <- count
					}
				}
				
				index1 <- results[g2, g1][[1]][, "index1"]
				index2 <- results[g2, g1][[1]][, "index2"]
				starts1 <- results[g2, g1][[1]][, "start1"]
				starts2 <- results[g2, g1][[1]][, "start2"]
			}
			
			# extend blocks
			p <- paste(index1, index2)
			t <- tapply(seq_along(p), p, c, simplify=FALSE)
			index1 <- match(results[g2, g1][[1]][, "index1"],
				INDEX1) - 1L
			index2 <- match(results[g2, g1][[1]][, "index2"],
				INDEX2) - 1L
			for (i in seq_along(t)) {
				s <- t[[i]] # subset with the same indices
				
				o <- order(starts1[s])
				r <- order(o) # rank
				o1p <- r - 1L
				o1p[which.min(o1p)] <- NA # which is zero
				o1p <- o[o1p] # index of previous start
				o1n <- r + 1L
				o1n <- o[o1n] # index of next start
				
				o <- order(starts2[s])
				r <- order(o) # rank
				o2p <- r - 1L
				o2p[which.min(o2p)] <- NA # which is zero
				o2p <- o[o2p] # index of previous start
				o2n <- r + 1L
				o2n <- o[o2n] # index of next start
				
				.Call("extendSegments",
					results[g2, g1][[1]],
					w1,
					w2,
					s1,
					s2,
					o1p - 1L,
					o1n - 1L,
					o2p - 1L,
					o2n - 1L,
					s - 1L,
					dropScore,
					index1,
					index2,
					PACKAGE="DECIPHER")
			}
			
			# extend outer anchors to match extended blocks
			strand <- results[g2, g1][[1]][, "strand"]
			if (length(strand) > 0) {
				first_hit <- results[g2, g1][[1]][, "first_hit"]
				last_hit <- results[g2, g1][[1]][, "last_hit"]
				starts1 <- results[g2, g1][[1]][, "start1"]
				ends1 <- results[g2, g1][[1]][, "end1"]
				starts2 <- results[g2, g1][[1]][, "start2"]
				ends2 <- results[g2, g1][[1]][, "end2"]
				
				# adjust width, start1, and start2 of first_hit
				results[g1, g2][[1]][first_hit, "width"] <- results[g1, g2][[1]][first_hit, "start1"] + results[g1, g2][[1]][first_hit, "width"] - starts1
				results[g1, g2][[1]][first_hit, "start1"] <- starts1
				results[g1, g2][[1]][first_hit, "start2"] <- ifelse(strand==0,
					starts2,
					ends2)
				
				# adjust width and start2 of last_hit
				results[g1, g2][[1]][last_hit, "width"] <- ends1 - results[g1, g2][[1]][last_hit, "start1"] + 1L
				results[g1, g2][[1]][last_hit, "start2"] <- ifelse(strand==0,
					ends2 - results[g1, g2][[1]][last_hit, "width"] + 1L,
					starts2 + results[g1, g2][[1]][last_hit, "width"] - 1L)
			}
			
			if (verbose) {
				its <- its + 0.5
				setTxtProgressBar(pBar, its/tot)
			}
		}
		
		store[[g1]] <- list()
	}
	
	if (!exists("w2",
		envir=environment(),
		inherits=FALSE)) {
		s2 <- SearchDB(dbConn,
			tblName=tblName,
			identifier=identifier[length(identifier)],
			type="DNAStringSet",
			removeGaps="all",
			processors=processors,
			verbose=FALSE)
		w2 <- width(s2)
	}
	results[l, l] <- list(widths=w2)
	
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	class(results) <- "Synteny"
	
	return(results)
}
