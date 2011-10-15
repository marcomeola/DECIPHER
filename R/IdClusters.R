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

IdClusters <- function(myDistMatrix,
	method="UPGMA",
	cutoff=-Inf,
	showPlot=FALSE,
	asDendrogram=FALSE,
	myDNAStringSet=NULL,
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
	
	w1 <- which(is.infinite(myDistMatrix))
	if (length(w1) > 0) {
		if (verbose)
			warning("\n\nDistance Matrix contains infinite values.\n",
				"Replaced infinite values with max distance >= 1.\n")
		myDistMatrix[w1] <- NA
	}
	
	w2 <- which(is.na(myDistMatrix))
	if (length(w2) > 0) {
		if (verbose &&
			length(w2) > length(w1))
			warning("\n\nDistance Matrix contains NA values.\n",
				"Replaced NA values with max distance >= 1.\n")
		max.dist <- max(myDistMatrix, na.rm=TRUE)
		if (max.dist <= 1)
			myDistMatrix[w2] <- 1
		else
			myDistMatrix[w2] <- max.dist
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
			cutoff[1],#-Inf
			verbose,
			pBar,
			PACKAGE="DECIPHER")
		
		myClusters <- .Call("clusterML",
			myClusters,
			myDNAStringSet,
			cutoff[1],
			verbose,
			pBar,
			PACKAGE="DECIPHER")
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
				m <- paste(m, ", ", max(myClusters[,9:10]), sep="")
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
		setTxtProgressBar(pBar,100)
		close(pBar)
		cat("\nGrouped into",
			m,
			"clusters.")
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
			digits=1))
		cat("\n")
	}
	if (asDendrogram)
		return(d)
	else
		return(myClusters)
}