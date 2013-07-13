# below function from stats package
.memberDend <- function(x) {
	
    r <- attr(x,"x.member")
    if(is.null(r)) {
	r <- attr(x,"members")
	if(is.null(r)) r <- 1L
    }
    return(r)
}

# below function modified from stats package
to.dendrogram <- function (object,
	hang = -1,
	...) {
	
    z <- list()
    oHgts <- object$lengths
    nMerge <- length(oHgt <- object$height)
    if (nMerge != nrow(object$merge))
	stop("'merge' and 'height' do not fit!")
    hMax <- oHgt[nMerge]
    cMax <- max(object$clusters)
    
    one <- 1L
    two <- 2L
    for (k in 1L:nMerge) {
		x <- object$merge[k, ]# no sort() anymore!
		neg <- x < 0
		if (all(neg)) {			# two leaves
			    zk <- as.list(-x)
			    attr(zk, "members") <- two
			    attr(zk, "midpoint") <- 0.5 # mean( c(0,1) )
			    objlabels <- object$labels[-x]
			    attr(zk[[1L]], "label") <- objlabels[1L]
			    attr(zk[[2L]], "label") <- objlabels[2L]
			    attr(zk[[1L]], "members") <- attr(zk[[2L]], "members") <- one
			    attr(zk[[1L]], "height") <- oHgt[k] - oHgts[k,1]
			    attr(zk[[2L]], "height") <- oHgt[k] - oHgts[k,2]
			    attr(zk[[1L]], "leaf") <- attr(zk[[2L]], "leaf") <- TRUE
			} else if (any(neg)) {		# one leaf, one node
			    X <- as.character(x)
			    ## Originally had "x <- sort(..) above => leaf always left, x[1L];
			    ## don't want to assume this
			    isL <- x[1L] < 0 ## is leaf left?
			    zk <-
				if(isL) list(-x[1L], z[[X[2L]]])
				else	list(z[[X[1L]]], -x[2L])
			    attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + one
			    attr(zk, "midpoint") <-
	        	        (.memberDend(zk[[1L]]) + attr(z[[X[1 + isL]]], "midpoint"))/2
			   	attr(zk[[2 - isL]], "members") <- one
			    attr(zk[[2 - isL]], "height") <- oHgt[k] - oHgts[k,2 - isL]
			    attr(zk[[2 - isL]], "label") <- object$labels[-x[2 - isL]]
			    attr(zk[[2 - isL]], "leaf") <- TRUE
			} else {				# two nodes
		    	x <- as.character(x)
		    	zk <- list(z[[x[1L]]], z[[x[2L]]])
		    	attr(zk, "members") <- attr(z[[x[1L]]], "members") +
				attr(z[[x[2L]]], "members")
				attr(zk, "midpoint") <- (attr(z[[x[1L]]], "members") +
						     attr(z[[x[1L]]], "midpoint") +
						     attr(z[[x[2L]]], "midpoint"))/2
		}
		attr(zk, "height") <- oHgt[k]
		z[[k <- as.character(k)]] <- zk
    }
    z <- z[[k]]
    class(z) <- "dendrogram"
    z
}

.organizeClusters <- function(myClusters,
	dNames,
	o) {
	
	l <- length(dNames)
	clusters <- data.frame(cluster=integer(l),
		row.names=dNames)
	k <- length(myClusters[,1])
	for (i in 1:k) {
		if (myClusters[i,1] < 0)
			clusters$cluster[-1*myClusters[i,1]] <-
				as.integer(myClusters[i,9])
		if (myClusters[i,2] < 0)
			clusters$cluster[-1*myClusters[i,2]] <-
				as.integer(myClusters[i,10])
	}
	
	# order the cluster numbers to match
	# the order of the dendrogram
	temp <- 0
	l <- max(clusters$cluster)
	clustersTemp <- clusters
	v <- vector(mode="numeric",length=l)
	j <- 0
	for (i in 1:length(o)) {
		if (clusters$cluster[o[i]] != temp &
			length(which(v==clusters$cluster[o[i]]))==0) {
			temp <- clusters$cluster[o[i]]
			j <- j + 1
			v[j] <- temp
		}
	}
	for (k in 1:l) {
		w <- which(clusters$cluster == v[k])
		clustersTemp$cluster[w] <- k
	}
	clusters <- clustersTemp
	
	return(clusters)
}

.organizeClustersFast <- function(myClusters,
	dNames) {
	
	l <- length(dNames)
	clusters <- data.frame(cluster=integer(l),
		row.names=dNames)
	k <- length(myClusters[,1])
	for (i in 1:k) {
		if (myClusters[i,1] < 0)
			clusters$cluster[-1*myClusters[i,1]] <-
				as.integer(myClusters[i,9])
		if (myClusters[i,2] < 0)
			clusters$cluster[-1*myClusters[i,2]] <-
				as.integer(myClusters[i,10])
	}
	
	return(clusters)
}

.rates <- function(alpha, nBins) {
	
	# Determine rates based on alpha and the number of bins
	# bins roots normalized to 1 of the General Laguerre Quadrature
	# first nBins elements are rates with mean 1
	# second nBins elements are probabilities with sum 1
	
	findRoots <- function(alpha, nBins) {
		
		# Determine rates based on Gamma's alpha and the number of bins
		# bins roots normalized to 1 of the General Laguerre Polynomial (GLP)
		
		coeff  <- integer(nBins + 1)
		for (i in 0:nBins) {
			coeff[i + 1] <- (-1)^i*choose(nBins + alpha, nBins - i)/factorial(i)
		}
		
		return(sort(Re(polyroot(coeff))))
	}
	
	roots <- findRoots(alpha - 1, nBins)
	
	Laguerre <- function(x, alpha, degree) {
		y <- 0
		for (i in 0:degree) {
			y <- y + (-1)^i*choose(degree + alpha, degree - i)*x^i/factorial(i)
		}
		return(y)
	}
	
	weights <- numeric(nBins)
	f <- prod(1 + (alpha - 1)/(1:nBins))
	
	for (i in 1:nBins) {
		weights[i] <- f*roots[i]/((nBins + 1)^2*Laguerre(roots[i],
			alpha - 1,
			nBins + 1)^2)
	}
	
	roots <- roots/alpha
	
	return(c(roots, weights))
}

.optimizeModel <- function(myClusters,
	model,
	myDNAStringSet,
	N) {
	
	verbose <- FALSE
	pBar <- NULL
	
	if (model=="JC69") {
		LnL <- .Call("clusterML",
			myClusters,
			myDNAStringSet,
			c(0.25, 0.25, 0.25, 0.25, 1, 1, 1, 1),
			verbose,
			pBar,
			PACKAGE="DECIPHER")
		K <- 2*dim(myClusters)[1] - 1
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(NA, NA, NA, NA, NA, NA, NA, LnL, AICc, BIC))
	} else if (model=="JC69+G") {
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(0.25, 0.25, 0.25, 0.25, 1, 1, .rates(params, 6)),
				verbose,
				pBar,
				PACKAGE="DECIPHER")
		}
		o <- optimize(f, c(0.001, 500), tol=1e-4)
		K <- 2*dim(myClusters)[1]
		AICc <- 2*K + 2*o$objective + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*o$objective + K*log(N)
		return(c(NA, NA, NA, NA, NA, NA, o$minimum, o$objective, AICc, BIC))
	} else if (model=="K80" || model=="K80+G") {
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(0.25, 0.25, 0.25, 0.25, params, params, 1, 1),
				verbose,
				pBar,
				PACKAGE="DECIPHER")
		}
		o <- optimize(f, c(0, 10), tol=1e-4)
		
		if (model=="K80+G") {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(0.25, 0.25, 0.25, 0.25, o$minimum, o$minimum, .rates(params, 6)),
					verbose,
					pBar,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 1
		} else {
			a <- NA
			LnL <- o$objective
			K <- 2*dim(myClusters)[1]
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(NA, NA, NA, NA, rep(o$minimum, 2), a, LnL, AICc, BIC))
	} else if (model=="F81" || model=="F81+G") {
		f <- function(params) {
			if (sum(params) > 1)
				return(1e9)
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(params[1], params[2], params[3], 1 - sum(params), 1, 1, 1, 1),
				verbose,
				pBar,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(rep(0.25, 3),
			f,
			upper=rep(1, 3),
			lower=rep(0, 3),
			control=list(rel.tol=1e-4))
		
		if (model=="F81+G") {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(o$par[1], o$par[2], o$par[3], 1 - sum(o$par), 1, 1, .rates(params, 6)),
					verbose,
					pBar,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 3
		} else {
			K <- 2*dim(myClusters)[1] + 2
			a <- NA
			LnL <- o$objective
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(o$par, 1-sum(o$par), NA, NA, a, LnL, AICc, BIC))
	} else if (model=="HKY85" || model=="HKY85+G") {
		f <- function(params) {
			if (sum(params[1:3]) > 1)
				return(1e9)
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(params[1], params[2], params[3], 1 - sum(params[1:3]), params[4], params[4], 1, 1),
				verbose,
				pBar,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(c(rep(0.25, 3), 1),
			f,
			upper=c(rep(1, 3), 10),
			lower=rep(0, 4),
			control=list(rel.tol=1e-4))
		
		if (model=="HKY85+G") {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(o$par[1], o$par[2], o$par[3], 1 - sum(o$par[1:3]), o$par[4], o$par[4], .rates(params, 6)),
					verbose,
					pBar,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 4
		} else {
			a <- NA
			LnL <- o$objective
			K <- 2*dim(myClusters)[1] + 3
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(o$par[1:3], 1-sum(o$par[1:3]), o$par[4], o$par[4], a, LnL, AICc, BIC))
	} else if (model=="T92" || model=="T92+G") {
		f <- function(params) {
			if (params[1] > 0.5)
				return(1e9)
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(params[1], rep((1 - 2*params[1])/2, 2), params[1], params[2], params[2], 1, 1),
				verbose,
				pBar,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(c(0.25, 1),
			f,
			upper=c(1, 10),
			lower=c(0, 0),
			control=list(rel.tol=1e-4))
		
		if (model=="T92+G") {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(o$par[1], rep((1 - 2*o$par[1])/2, 2), o$par[1], o$par[2], o$par[2], .rates(params, 6)),
					verbose,
					pBar,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 2
		} else {
			a <- NA
			LnL <- o$objective
			K <- 2*dim(myClusters)[1] + 1
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(o$par[1], rep((1- 2*o$par[1])/2, 2), o$par[1], o$par[2], o$par[2], a, LnL, AICc, BIC))
	} else if (model=="TN93" || model=="TN93+G") {
		f <- function(params) {
			if (sum(params[1:3]) > 1)
				return(1e9)
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(params[1], params[2], params[3], 1 - sum(params[1:3]), params[4], params[5], 1, 1),
				verbose,
				pBar,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(c(rep(0.25, 3), 1, 1),
			f,
			upper=c(rep(1, 3), 10, 10),
			lower=c(rep(0, 3), 0, 0),
			control=list(rel.tol=1e-4))
		
		if (model=="TN93+G") {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(o$par[1], o$par[2], o$par[3], 1 - sum(o$par[1:3]), o$par[4], o$par[5], .rates(params, 6)),
					verbose,
					pBar,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 5
		} else {
			a <- NA
			LnL <- o$objective
			K <- 2*dim(myClusters)[1] + 4
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(o$par[1:3], 1-sum(o$par[1:3]), o$par[4], o$par[5], a, LnL, AICc, BIC))
	}
}

MODELS <- c("JC69",
	"JC69+G",
	"K80",
	"K80+G",
	"F81",
	"F81+G",
	"HKY85",
	"HKY85+G",
	"T92",
	"T92+G",
	"TN93",
	"TN93+G")

IdClusters <- function(myDistMatrix,
	method="UPGMA",
	cutoff=-Inf,
	showPlot=FALSE,
	asDendrogram=FALSE,
	myDNAStringSet=NULL,
	model=MODELS,
	add2tbl=FALSE,
	dbFile=NULL,
	verbose=TRUE) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	METHODS <- c("NJ","UPGMA", "ML", "complete", "single", "average")
	method <- pmatch(method, METHODS)
	if (is.na(method))
		stop("Invalid method.  Choose either ML, NJ, complete, single, or UPGMA (average).")
	if (method == -1)
		stop("Ambiguous method.  Choose either ML, NJ, complete, single, or UPGMA (average).")
	if (method == 6)
		method <- 2
	if (method==3 && length(model) < 1)
		stop("No model(s) specified.")
	if (method==3 && !is.character(model))
		stop("model must be a character vector.")
	if (method==3 && !all(model %in% MODELS))
		stop(paste("Available models are:",
			paste(MODELS, collapse=", "),
			collapse=" "))
	model <- unique(model)
	if (!is.numeric(cutoff))
		stop("cutoff must be a numeric.")
	if (!is.logical(showPlot))
		stop("showPlot must be a logical.")
	if (!is.logical(asDendrogram))
		stop("asDendrogram must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (is.character(add2tbl) || add2tbl)
		if (!is.character(dbFile) &&
			!inherits(dbFile,"SQLiteConnection"))
				stop("dbFile must be a character string or connection if add2tbl.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	dim <- dim(myDistMatrix)
	if (dim[2]!=dim[1])
		stop("\nYour distance matrix is not square.")
	dim <- dim[1]
	if (dim < 2)
		stop("\nYour distance matrix is too small.")
	if (method == 3) {
		if (is.null(myDNAStringSet))
			stop("A DNAStringSet must be provided when method is ML.")
		if (!is(myDNAStringSet, "DNAStringSet"))
			stop("myDNAStringSet must be a DNAStringSet.")
		if (length(myDNAStringSet)!=dim(myDistMatrix)[1])
			stop("myDistMatrix must have as many rows as the number of sequences.")
	}
	
	w1 <- which(is.infinite(myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)]))
	if (length(w1) > 0) {
		if (verbose)
			warning("\n\nDistance Matrix contains infinite values.\n",
				"Replaced infinite values with max distance >= 1.\n")
		myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)][w1] <- NA
	}
	
	w2 <- which(is.na(myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)]))
	if (length(w2) > 0) {
		if (verbose &&
			length(w2) > length(w1))
			warning("\n\nDistance Matrix contains NA values.\n",
				"Replaced NA values with max distance >= 1.\n")
		max.dist <- max(myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)], na.rm=TRUE)
		if (max.dist <= 1) {
			myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)][w2] <- 1
		} else {
			myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)][w2] <- max.dist
		}
	}
	
	if ((length(cutoff) > 1) && (showPlot || asDendrogram)) {
		cutoff <- cutoff[1]
		warning("Only the first cutoff used since showPlot or asDendrogram.")
	}
	
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=3)
	} else {
		pBar <- NULL
	}
	
	if (method==1)
		myClusters <- .Call("clusterNJ",
			myDistMatrix,
			cutoff[1],
			verbose,
			pBar,
			PACKAGE="DECIPHER")
	
	if (method==2 ||
		method==4 ||
		method==5 ||
		method==6)
		myClusters <- .Call("clusterUPGMA",
			myDistMatrix,
			cutoff[1],
			method,
			verbose,
			pBar,
			PACKAGE="DECIPHER")
	
	if (method==3) {
		# create NJ for use as guide tree
		myClusters <- .Call("clusterNJ",
			myDistMatrix,
			cutoff=-Inf,
			verbose,
			pBar,
			PACKAGE="DECIPHER")
		
		m <- matrix(NA,
			nrow=length(model),
			ncol=10,
			dimnames=list(model,
				c("FreqA", "FreqC", "FreqG", "FreqT",
					"A2G", "C2T", "alpha",
					"-LnL", "AICc", "BIC")))
		
		if (!(length(model)==1 && model[1]=="JC69")) {
			N <- sum(ifelse(apply(consensusMatrix(myDNAStringSet),
						2,
						function(x)
							return(length(which(x[1:14] > 0)))) > 1,
					1,
					0))
			if (verbose)
				cat("\n")
			for (i in 1:length(model)) {
				m[model[i],] <- .optimizeModel(myClusters,
					model[i],
					myDNAStringSet,
					N)
				if (verbose)
					cat("\n", model[i],
						":", paste(rep(" ",
								max(nchar(rownames(m))) - nchar(model[i]) + 1),
							collapse=""),
						"-ln(L)=", round(m[model[i], "-LnL"], 0),
						",\tAICc=", round(m[model[i], "AICc"], 0),
						",\tBIC=", round(m[model[i], "BIC"], 0),
						sep="")
			}
		}
		
		if (length(unique(model)) > 1) { # choose the best model
			w <- which.min(m[,"BIC"])
			l <- logical(nrow(m))
			l[w] <- TRUE
			m <- subset(m, l)
		}
		
		if (verbose) {
			params <- formatC(round(m, 3),
				digits=3,
				format="f")
			cat("\n\nThe selected model was:  ",
				rownames(m),
				ifelse(any(grepl("NA", params[1], fixed=TRUE)),
					"",
					"\n\nWith optimized parameters:"),
				ifelse(grepl("NA", params[1], fixed=TRUE),
					"",
					ifelse(grepl("T92", rownames(m), fixed=TRUE),
						paste0("\nFrequency(A) = Frequency(T) = ", params[1],
							"\nFrequency(C) = Frequency(G) = ", params[2]),
						paste0("\nFrequency(A) = ", params[1],
							"\nFrequency(C) = ", params[2],
							"\nFrequency(G) = ", params[3],
							"\nFrequency(T) = ", params[4]))),
				ifelse(grepl("NA", params[5], fixed=TRUE),
					"",
					ifelse(grepl("TN93", rownames(m), fixed=TRUE),
						paste0("\nRate A <-> G = ", params[5],
							"\nRate C <-> T = ", params[6],
							"\nTransversion rates = 1"),
						paste0("\nTransition rates = ", params[5],
							"\nTransversion rates = 1"))),
				ifelse(grepl("NA", params[7], fixed=TRUE),
					"",
					paste0("\nAlpha = ", params[7])),
				sep="")
		}
		
		model_params <- as.numeric(m[1:7])
		w <- which(is.na(model_params))
		if (length(w) > 0) {
			model_params[w] <- c(0.25, 0.25, 0.25, 0.25, 1, 1, NA)[w]
		}
		if (is.na(model_params[7])) {
			model_params <- c(model_params[1:6], 1, 1)
		} else {
			model_params <- c(model_params[1:6], .rates(model_params[7], 6))
		}
		
		# neglect clusters of nearly identical sequences
		w <- which(myClusters[, 4] > 0.0001 |
			myClusters[, 5] > 0.0001)
		
		# given log(1 + branch length) return adjusted lengths
		adjustTree <- function(x) {
			x <- exp(x) - 1
			myClusters[w, 4] <- head(x, n=length(x)/2)
			myClusters[w, 5] <- tail(x, n=length(x)/2)
			return(myClusters)
		}
		
		# given myClusters return adjusted heights
		adjustTreeHeights <- function(myClusters) {
			cumHeight <- numeric(max(myClusters[, 3]))
			for (i in 1:dim(myClusters)[1]) {
				if (myClusters[i, 1] < 0 && myClusters[i, 2] < 0) {
					cumHeight[myClusters[i, 3]] <- max(myClusters[i, 4], myClusters[i, 5])
					myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
				} else if (myClusters[i, 1] > 0 && myClusters[i, 2] > 0) {
					cumHeight[myClusters[i, 3]] <- max(myClusters[i, 4] + cumHeight[myClusters[i, 1]],
						myClusters[i, 5] + cumHeight[myClusters[i, 2]])
					myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
				} else if (myClusters[i, 1] > 0) {
					cumHeight[myClusters[i, 3]] <- cumHeight[myClusters[i, 1]] + myClusters[i, 4]
					if (myClusters[i, 5] > cumHeight[myClusters[i, 3]])
						cumHeight[myClusters[i, 3]] <- myClusters[i, 5]
					myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
				} else {
					cumHeight[myClusters[i, 3]] <- cumHeight[myClusters[i, 2]] + myClusters[i, 5]
					if (myClusters[i, 4] > cumHeight[myClusters[i, 3]])
						cumHeight[myClusters[i, 3]] <- myClusters[i, 4]
					myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
				}
			}
			
			myClusters <- .Call("adjustHeights",
				myClusters,
				PACKAGE="DECIPHER")
			return(myClusters)
		}
		
		# print progress of likelihood maximization
		.startLnL <- Inf
		.bestLnL <- Inf
		if (verbose) {
			setTxtProgressBar(pBar,100)
			close(pBar)
			cat("\nMaximizing Likelihood of Tree:\n")
			flush.console()
			printLine <- function(value, percentComplete) {
				cat("\r-ln(Likelihood) = ",
					formatC(round(value, 0),
						digits=0,
						format="f"),
					" (",
					formatC(round(-100*(value - .startLnL)/.startLnL,
							2),
						digits=2,
						format="f"),
					"% increase)",
					sep="")
				flush.console()
				invisible(value)
			}
		}
		
		# given log(1 + branch length) return -LnL
		maximizeLikelihood <- function(x) {
			LnL <- .Call("clusterML",
				adjustTree(x),
				myDNAStringSet,
				model_params,
				verbose,
				pBar,
				PACKAGE="DECIPHER")
			
			if (LnL < .bestLnL) {
				.bestLnL <<- LnL
				params <<- x
				if (verbose) {
					if (is.infinite(.startLnL))
						.startLnL <<- LnL
					printLine(LnL)
				}
			}
			
			return(LnL)
		}
		
		# maximize likelihood of branch lengths
		if (length(w) > 0) {
			params <- log(c(myClusters[w, 4], myClusters[w, 5]) + 1)
			#for (i in 1:length(params)) {
			#	f <- function(l) {
			#		params[i] <- l
			#		maximizeLikelihood(params)
			#	}
			#	params[i] <- optimize(f,
			#		c(0, ceiling(params[i] + 0.5)),
			#		tol=.0001)$minimum
			#	myClusters <- adjustTree(params)
			#}
			o <- nlminb(params,
				maximizeLikelihood,
				control=list(rel.tol=1e-2,
					xf.tol=1e-2),
				lower=rep(0.000001, length(w)),
				upper=ceiling(params + 0.5))
			myClusters <- adjustTree(params)
			myClusters <- adjustTreeHeights(myClusters)
			myClusters <- .Call("reclusterNJ",
				myClusters,
				cutoff[1],
				PACKAGE="DECIPHER")
		}
		
		if (verbose) {
			printLine(.bestLnL, 100)
			cat("\n")
			flush.console()
		}
	}
	
	m <- max(myClusters[,9:10])
	
	if (showPlot || asDendrogram) {
		# create a dendrogram object
		myClustersList <- list()
		if (is.null(dimnames(myDistMatrix)[[1]])) {
			myClustersList$labels <- 1:(dim(myClusters)[1] + 1)
		} else {
			myClustersList$labels <- dimnames(myDistMatrix)[[1]]
		}
		myClustersList$merge <- matrix(myClusters[,7:8], ncol=2)
		myClustersList$height <- matrix(myClusters[,6], ncol=1)
		myClustersList$lengths <- matrix(myClusters[,4:5], ncol=2)
		myClustersList$clusters <- matrix(myClusters[,9:10], ncol=2)
		if (dim > 100)
			fontSize <- .6
		else if (dim > 70)	
			fontSize <- .7
		else if (dim > 40)
			fontSize <- .8
		else
			fontSize <- .9
		if (dim > 300) {
			leaves <- "none"
		} else {
			leaves <- "perpendicular"
		}
		d <- to.dendrogram(myClustersList)
		
		# specify the order of clusters that
		# will match the plotted dendrogram
		orderDendrogram <- order.dendrogram(d)
		
		myClusters <- .organizeClusters(myClusters, myClustersList$labels, orderDendrogram)
		
		# create a visibily different vector of colors
		cl <- colors()
		v1 <- c(117,254,73,69,152,51,26,450,503,596,652,610,563,552,97)
		r <- cl[v1]
		
		# color edges by cluster
		colEdge <- function(n, myClusters, colors) {
			if (is.leaf(n)) {
				a <- attributes(n)
				num <- myClusters$cluster[which(myClustersList$labels==as.character(a$label))]
				attr(n, "edgePar") <- list(col=colors[num %% 15 + 1])
			}
			n
		}
		if (is.finite(cutoff))
			d <- dendrapply(d, colEdge, myClusters, r)
	} else { # not plotting
		if (is.null(dimnames(myDistMatrix)[[1]])) {
			dNames <- 1:(dim(myClusters)[1] + 1)
		} else {
			dNames <- dimnames(myDistMatrix)[[1]]
		}
		
		# do not number clusters by order of appearance
		c <- .organizeClustersFast(myClusters, dNames)
		if (length(cutoff) > 1) {
			names(c) <- paste("cluster",
				gsub("\\.","_",100*cutoff[1]),
				METHODS[method],
				sep="")
			for (i in 2:length(cutoff)) {
				if (method==2 ||
					method==4 ||
					method==5 ||
					method==6)
					myClusters <- .Call("reclusterUPGMA",
						myClusters,
						cutoff[i],
						PACKAGE="DECIPHER")
				if (method==1 ||
					method==3)
					myClusters <- .Call("reclusterNJ",
						myClusters,
						cutoff[i],
						PACKAGE="DECIPHER")
				x <- .organizeClustersFast(myClusters, dNames)
				names(x) <- paste("cluster",
					gsub("\\.","_",100*cutoff[i]),
					METHODS[method],
					sep="")
				c <- cbind(c, x)
				m <- c(m, max(myClusters[,9:10]))
			}
		}
		myClusters <- c
	}
	
	if (showPlot)
		plot(d,
			horiz=FALSE,
			leaflab=leaves,	
			nodePar = list(lab.cex=fontSize, pch = NA))
	
	if (is.character(add2tbl) || add2tbl)
		Add2DB(myData=myClusters,
			dbFile=dbFile,
			tblName=ifelse(is.character(add2tbl),add2tbl,"DNA"),
			verbose=FALSE)
	
	if (verbose) {
		if (method != 3) { # already closed pBar
			setTxtProgressBar(pBar,100)
			close(pBar)
		}
		#cat("\nGrouped into",
		#	paste(unique(m), collapse=", "),
		#	"clusters.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to ",
				ifelse(is.character(add2tbl),add2tbl,"DNA"),
				":  \"",
				ifelse(asDendrogram,
					"cluster",
					paste(names(myClusters),
						sep="",
						collapse="\", \"")),
				"\".",
				sep="")
		
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	if (asDendrogram)
		return(d)
	else
		return(myClusters)
}