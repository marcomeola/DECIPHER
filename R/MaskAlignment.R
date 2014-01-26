.centerPoint <- function(x, size) {
	
	# returns the center-point moving average of x
	# using a window size of size elements
	
	l <- length(x)
	index <- 1:l
	boundsL <- index - size
	boundsL[1:size] <- 1
	boundsL[(l - size + 1):l] <- boundsL[(l - size + 1):l] - 1:size
	boundsR <- index + size
	boundsR[(l - size + 1):l] <- l
	boundsR[1:size] <- boundsR[1:size] + size:1
	
	avgs <- numeric(length(x))
	for (i in index) {
		avgs[i] <- mean(x[boundsL[i]:boundsR[i]])
	}
	return(avgs)
}

MaskAlignment <- function(myXStringSet,
	windowSize=5,
	threshold=1.0,
	maxFractionGaps=0.2,
	showPlot=FALSE) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet") && !is(myXStringSet, "AAStringSet"))
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	if (!is.numeric(windowSize))
		stop("windowSize must be a numeric.")
	if (floor(windowSize)!=windowSize)
		stop("windowSize must be a whole number.")
	if (windowSize < 1)
		stop("windowSize must be greater than zero.")
	u <- unique(width(myXStringSet))
	if (length(u) > 1)
		stop("myXStringSet must be aligned.")
	if (!is.numeric(maxFractionGaps))
		stop("maxFractionGaps must be a numeric.")
	if (maxFractionGaps > 1 || maxFractionGaps < 0)
		stop("maxFractionGaps must be between 0 and 1 inclusive.")
	
	f <- function(x) {
		w <- which(x != 0)
		if (length(w) > 0) {
			x <- x[w]/sum(x[w])
			bits <- 2 + sum(x*log(x, 2))
		} else {
			bits <- 0
		}
		return(bits)
	}
	
	if (is(myXStringSet, "AAStringSet")) {
		pwm <- .Call("consensusProfileAA",
			myXStringSet,
			rep(1, length(myXStringSet)),
			PACKAGE="DECIPHER")
		a <- apply(pwm[1:25,], 2, f)
	} else {
		pwm <- .Call("consensusProfile",
			myXStringSet,
			rep(1, length(myXStringSet)),
			PACKAGE="DECIPHER")
		a <- apply(pwm[1:4,], 2, f)
	}
	
	cm <- consensusMatrix(myXStringSet, as.prob=TRUE)["-",]
	gaps <- which(cm > maxFractionGaps)
	
	if (windowSize*2 + 1 > length(a) - length(gaps))
		stop("windowSize is too large.")
	
	if (length(gaps) > 0) {
		a2 <- a[-gaps]
		c <- .centerPoint(a2, windowSize)
	} else {
		a2 <- a
		c <- .centerPoint(a2, windowSize)
	}
	
	W <- which(c < threshold)
	if (length(W) > 0) {
		if (length(W) > 1) {
			w <- which((W[2:length(W)] - 1) != W[1:(length(W) - 1)])
		} else {
			w <- integer()
		}
		index <- W[c(1, w + 1)]
		
		indicies <- list()
		below_threshold <- as.integer(a2 < threshold)
		rev_below_threshold <- rev(below_threshold)
		for (i in 1:length(index)) {
			rights <- .Call("multiMatch",
				below_threshold,
				1L,
				index[i])
			lefts <- length(below_threshold) - .Call("multiMatch",
				rev_below_threshold,
				1L,
				length(below_threshold) - index[i] + 1L) + 1
			indicies[[i]] <- c(lefts, rights)
		}
		
		W <- c(W, unlist(indicies))
		
		if (length(gaps) > 0) {
			W <- ((1:length(a))[-gaps])[W]
		}
		W <- c(W, gaps)
		W <- unique(sort(W))
		
		if (length(W) > 1) {
			w <- which((W[2:length(W)] - 1) != W[1:(length(W) - 1)])
		} else {
			w <- integer()
		}
		starts <- W[c(1, w + 1)]
		ends <- W[c(w, length(W))]
		
		if (is(myXStringSet, "DNAStringSet")) {
			myXStringSet <- DNAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(starts, ends), "NormalIRanges"))
		} else if (is(myXStringSet, "RNAStringSet")) {
			myXStringSet <- RNAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(starts, ends), "NormalIRanges"))
		} else { # AAStringSet
			myXStringSet <- AAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(starts, ends), "NormalIRanges"))
		}
	} else if (length(gaps) > 0) {
		if (length(gaps) > 1) {
			w <- which((gaps[2:length(gaps)] - 1) != gaps[1:(length(gaps) - 1)])
		} else {
			w <- integer()
		}
		starts <- gaps[c(1, w + 1)]
		ends <- gaps[c(w, length(gaps))]
		
		if (is(myXStringSet, "DNAStringSet")) {
			myXStringSet <- DNAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(starts, ends), "NormalIRanges"))
		} else if (is(myXStringSet, "RNAStringSet")) {
			myXStringSet <- RNAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(starts, ends), "NormalIRanges"))
		} else { # AAStringSet
			myXStringSet <- AAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(starts, ends), "NormalIRanges"))
		}
	} else {
		if (is(myXStringSet, "DNAStringSet")) {
			myXStringSet <- DNAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(), "NormalIRanges"))
		} else if (is(myXStringSet, "RNAStringSet")) {
			myXStringSet <- RNAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(), "NormalIRanges"))
		} else { # AAStringSet
			myXStringSet <- AAMultipleAlignment(myXStringSet,
				rowmask=as(IRanges(), "NormalIRanges"),
				colmask=as(IRanges(), "NormalIRanges"))
		}
	}
	
	if (showPlot) {
		if (length(gaps) > 0) {
			x <- (1:length(a))[-gaps]
		} else {
			x <- 1:length(a)
		}
		if (length(gaps) > 0) {
			org_mar <- par()$mar
			par(mar=org_mar + c(0, 0, 0, 2))
		}
		plot(x,
			c,
			type="l",
			ylim=c(0,2),
			ylab="Information Content (bits)",
			xlab="Column Position in Alignment")
		w <- which(!(W %in% gaps))
		if (length(w) > 0)
			points(W[w],
				a[W[w]],
				pch=20,
				col="red")
		if (length(W) < length(a)) {
			if (length(W) > 0) {
				x <- (1:length(a))[-W]
				y <- a[-W]
			} else {
				x <- 1:length(a)
				y <- a
			}
			points(x,
				y,
				pch=20,
				col="green")
		}
		if (length(gaps) > 0) {
			abline(h=maxFractionGaps,
				lty=2,
				col="orange")
			points(gaps,
				2*cm[gaps],
				pch=20,
				col="orange")
			axis(4,
				seq(0, 2, 0.5),
				labels=c("0", "25", "50", "75", "100"),
				col.axis="orange")
			mtext("Percent Gaps (%)", 4, line=3, col="orange")
			par(mar=org_mar)
		}
		abline(h=threshold,
			lty=2,
			col="red")
		legend("bottomright",
			c(ifelse(length(gaps) > 0, "Thresholds", "Threshold"),
				"Moving Avg.",
				"Kept Position",
				"Masked Pos.",
				"Masked Gaps"),
			bg="white",
			lty=c(2, 1, 0, 0, 0),
			pch=c(NA, NA, 20, 20, 20),
			col=c(ifelse(length(gaps) > 0, "black", "red"), "black", "green", "red", "orange"))
	}
	
	return(myXStringSet)
}
