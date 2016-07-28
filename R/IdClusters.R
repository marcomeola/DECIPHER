# below function modified from stats package
.memberDend <- function(x) {
	
	r <- attr(x, "members")
	if(is.null(r))
		r <- 1L
	
	return(r)
}

# below function modified from stats package
to.dendrogram <- function(object) {
	z <- list()
	oHgts <- object$lengths
	oHgt <- object$height
	nMerge <- length(oHgt)
	if (nMerge != nrow(object$merge))
		stop("'merge' and 'height' do not fit!")
	hMax <- oHgt[nMerge]
	cMax <- max(object$clusters)
	
	one <- 1L
	two <- 2L
	for (k in 1L:nMerge) {
		x <- as.integer(object$merge[k, ])
		neg <- x < 0
		if (all(neg)) { # two leaves
			zk <- as.list(-x)
			attr(zk, "members") <- two
			attr(zk, "midpoint") <- 0.5 # mean(c(0,1))
			objlabels <- object$labels[-x]
			attr(zk[[1L]], "label") <- objlabels[1L]
			attr(zk[[2L]], "label") <- objlabels[2L]
			attr(zk[[1L]], "members") <- attr(zk[[2L]], "members") <- one
			attr(zk[[1L]], "height") <- oHgt[k] - oHgts[k, 1]
			attr(zk[[2L]], "height") <- oHgt[k] - oHgts[k, 2]
			attr(zk[[1L]], "leaf") <- attr(zk[[2L]], "leaf") <- TRUE
		} else if (any(neg)) { # one leaf, one node
			X <- as.character(x)
			isL <- x[1L] < 0 # is leaf left?
			zk <-
			if (isL) {
				list(-x[1L], z[[X[2L]]])
			} else {
				list(z[[X[1L]]], -x[2L])
			}
			attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + one
			attr(zk, "midpoint") <- (.memberDend(zk[[1L]]) + attr(z[[X[1 + isL]]], "midpoint"))/2
			attr(zk[[2 - isL]], "members") <- one
			attr(zk[[2 - isL]], "height") <- oHgt[k] - oHgts[k, 2 - isL]
			attr(zk[[2 - isL]], "label") <- object$labels[-x[2 - isL]]
			attr(zk[[2 - isL]], "leaf") <- TRUE
		} else { # two nodes
			x <- as.character(x)
			zk <- list(z[[x[1L]]], z[[x[2L]]])
			attr(zk, "members") <- attr(z[[x[1L]]], "members") + attr(z[[x[2L]]], "members")
			attr(zk, "midpoint") <- (attr(z[[x[1L]]], "members") + attr(z[[x[1L]]], "midpoint") + attr(z[[x[2L]]], "midpoint"))/2
		}
		attr(zk, "height") <- oHgt[k]
		k <- as.character(k)
		z[[k]] <- zk
	}
	z <- z[[k]]
	class(z) <- "dendrogram"
	z
}

.collapse <- function(dend) {
	if (is.leaf(dend))
		return(dend)
	
	dend[[1]] <- .collapse(dend[[1]])
	dend[[2]] <- .collapse(dend[[2]])
	
	if (!is.leaf(dend[[1]])) {
		h1 <- attr(dend[[1]], "height")
	} else {
		h1 <- -Inf
	}
	if (!is.leaf(dend[[2]])) {
		h2 <- attr(dend[[2]], "height")
	} else {
		h2 <- -Inf
	}
	
	h <- attr(dend, "height")
	
	if (h==h1 || h==h2) { # make multifurcating
		m1 <- attr(dend[[1]], "members")
		m2 <- attr(dend[[2]], "members")
		m <- m1 + m2
		if (h==h1 && h==h2) {
			l1 <- length(dend[[1]])
			l2 <- length(dend[[2]])
			x <- vector("list", l1 + l2)
			for (i in seq_len(l1))
				x[i] <- dend[[1]][i]
			for (i in seq_len(l2))
				x[i + l1] <- dend[[2]][i]
		} else if (h==h1) {
			l <- length(dend[[1]])
			x <- vector("list", l + 1)
			for (i in seq_len(l))
				x[i] <- dend[[1]][i]
			x[l + 1] <- dend[-1]
		} else if (h==h2) {
			l <- length(dend[[2]])
			x <- vector("list", l + 1)
			x[1] <- dend[-2]
			for (i in seq_len(l))
				x[i + 1] <- dend[[2]][i]
		}
		dend <- x
		attr(dend, "height") <- h
		attr(dend, "members") <- m
		attr(dend, "midpoint") <- (m - 1)/2
		class(dend) <- "dendrogram"
	}
	
	return(dend)
}

.organizeClusters <- function(myClusters,
	dNames,
	o) {
	
	l <- length(dNames)
	clusters <- data.frame(cluster=integer(l),
		row.names=dNames)
	w <- which(myClusters[, 1] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 1]] <- as.integer(myClusters[w, 9])
	w <- which(myClusters[, 2] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 2]] <- as.integer(myClusters[w, 10])
	
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
	w <- which(myClusters[, 1] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 1]] <- as.integer(myClusters[w, 9])
	w <- which(myClusters[, 2] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 2]] <- as.integer(myClusters[w, 10])
	
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
	N,
	processors=1) {
	
	rates <- as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", model))
	model <- sub("([^+]*)(\\+G(\\d+))?", "\\1", model)
	
	if (model=="JC69" && is.na(rates)) {
		LnL <- .Call("clusterML",
			myClusters,
			myDNAStringSet,
			c(0.25, 0.25, 0.25, 0.25, 1, 1, 1, 1),
			integer(),
			numeric(),
			processors,
			PACKAGE="DECIPHER")
		K <- 2*dim(myClusters)[1] - 1
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(NA, NA, NA, NA, NA, NA, NA, LnL, AICc, BIC))
	} else if (model=="JC69") { # rates is an integer
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(0.25, 0.25, 0.25, 0.25, 1, 1, .rates(params, rates)),
				integer(),
				numeric(),
				processors,
				PACKAGE="DECIPHER")
		}
		o <- optimize(f, c(0.001, 500), tol=1e-4)
		K <- 2*dim(myClusters)[1]
		AICc <- 2*K + 2*o$objective + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*o$objective + K*log(N)
		return(c(NA, NA, NA, NA, NA, NA, o$minimum, o$objective, AICc, BIC))
	} else if (model=="K80") {
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(0.25, 0.25, 0.25, 0.25, params, params, 1, 1),
				integer(),
				numeric(),
				processors,
				PACKAGE="DECIPHER")
		}
		o <- optimize(f, c(0, 10), tol=1e-4)
		
		if (!is.na(rates)) {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(0.25, 0.25, 0.25, 0.25, o$minimum, o$minimum, .rates(params, rates)),
					integer(),
					numeric(),
					processors,
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
	} else if (model=="F81") {
		baseFreqs <- colSums(alphabetFrequency(myDNAStringSet))[1:4]
		baseFreqs <- baseFreqs/sum(baseFreqs)
#		f <- function(params) {
#			if (sum(params) > 0.99)
#				return(1e9)
#			LnL <- .Call("clusterML",
#				myClusters,
#				myDNAStringSet,
#				c(params[1], params[2], params[3], 1 - sum(params), 1, 1, 1, 1),
#				integer(),
#				numeric(),
#				processors,
#				PACKAGE="DECIPHER")
#			o <- optimize(adjustTreeHeight,
#				c(0.001, 1000),
#				tol=1e-2,
#				params=params)
#			return(o$minimum)
#		}
#		o <- nlminb(rep(0.25, 3),
#			f,
#			upper=rep(0.99, 3),
#			lower=rep(0.01, 3),
#			control=list(rel.tol=1e-4))
		
		if (!is.na(rates)) {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, 1, 1, .rates(params, rates)),
					integer(),
					numeric(),
					processors,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 3
		} else {
			K <- 2*dim(myClusters)[1] + 2
			a <- NA
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(baseFreqs, 1, 1, 1, 1),
				integer(),
				numeric(),
				processors,
				PACKAGE="DECIPHER")
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(baseFreqs, NA, NA, a, LnL, AICc, BIC))
	} else if (model=="HKY85") {
		baseFreqs <- colSums(alphabetFrequency(myDNAStringSet))[1:4]
		baseFreqs <- baseFreqs/sum(baseFreqs)
#		f <- function(params) {
#			if (sum(params) > 0.99)
#				return(1e9)
#			LnL <- .Call("clusterML",
#				myClusters,
#				myDNAStringSet,
#				c(params[1], params[2], params[3], 1 - sum(params[1:3]), 1, 1, 1, 1),
#				integer(),
#				numeric(),
#				processors,
#				PACKAGE="DECIPHER")
#		}
#		baseFreqs <- nlminb(rep(0.25, 3),
#			f,
#			upper=rep(0.99, 3),
#			lower=rep(0.01, 3),
#			control=list(rel.tol=1e-4))$par
		
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(baseFreqs, params, params, 1, 1),
				integer(),
				numeric(),
				processors,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(1,
			f,
			upper=10,
			lower=0.01,
			control=list(rel.tol=1e-4))
		
		if (!is.na(rates)) {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, o$par[1], o$par[1], .rates(params, rates)),
					integer(),
					numeric(),
					processors,
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
		return(c(baseFreqs, o$par[1], o$par[1], a, LnL, AICc, BIC))
	} else if (model=="T92") {
		baseFreqs <- colSums(alphabetFrequency(myDNAStringSet))[1:4]
		baseFreqs <- c((baseFreqs[1] + baseFreqs[4])/(2*sum(baseFreqs)),
			(baseFreqs[2] + baseFreqs[3])/(2*sum(baseFreqs)))
		baseFreqs <- c(baseFreqs[1], baseFreqs[2], baseFreqs[2], baseFreqs[1])
#		f <- function(params) {
#			if (params[1] > 0.49)
#				return(1e9)
#			LnL <- .Call("clusterML",
#				myClusters,
#				myDNAStringSet,
#				c(params, rep((1 - 2*params)/2, 2), params, 1, 1, 1, 1),
#				integer(),
#				numeric(),
#				processors,
#				PACKAGE="DECIPHER")
#		}
#		o <- nlminb(0.25,
#			f,
#			upper=0.49,
#			lower=0.01,
#			control=list(rel.tol=1e-4))
#		baseFreqs <- c(o$par, rep((1 - 2*o$par)/2, 2), o$par)
		
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(baseFreqs, params, params, 1, 1),
				integer(),
				numeric(),
				processors,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(1,
			f,
			upper=10,
			lower=0.01,
			control=list(rel.tol=1e-4))
		
		if (!is.na(rates)) {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, o$par[1], o$par[1], .rates(params, rates)),
					integer(),
					numeric(),
					processors,
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
		return(c(baseFreqs, o$par[1], o$par[1], a, LnL, AICc, BIC))
	} else if (model=="TN93") {
		baseFreqs <- colSums(alphabetFrequency(myDNAStringSet))[1:4]
		baseFreqs <- baseFreqs/sum(baseFreqs)
#		f <- function(params) {
#			if (sum(params[1:3]) > 0.99)
#				return(1e9)
#			LnL <- .Call("clusterML",
#				myClusters,
#				myDNAStringSet,
#				c(params[1], params[2], params[3], 1 - sum(params[1:3]), 1, 1, 1, 1),
#				integer(),
#				numeric(),
#				processors,
#				PACKAGE="DECIPHER")
#		}
#		o <- nlminb(rep(0.25, 3),
#			f,
#			upper=rep(0.99, 3),
#			lower=rep(0.01, 3),
#			control=list(rel.tol=1e-4))
#		baseFreqs <- c(o$par, 1 - sum(o$par))
		
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(baseFreqs, params[1], params[2], 1, 1),
				integer(),
				numeric(),
				processors,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(c(1, 1),
			f,
			upper=c(10, 10),
			lower=c(0.01, 0.01),
			control=list(rel.tol=1e-4))
		
		if (!is.na(rates)) {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, o$par[1], o$par[2], .rates(params, rates)),
					integer(),
					numeric(),
					processors,
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
		return(c(baseFreqs, o$par[1], o$par[2], a, LnL, AICc, BIC))
	}
}

.simultaneousBrent <- function(f, # scalar function
	a, # lower bounds
	b, # best guesses
	c, # upper bounds
	tol=1e-5, # accuracy of params change in x
	overallTol=1e0, # accuracy of overall scalar
	relTol=1e-2) { # accuracy of params change in y
	
	# performs optimization using Brent's method
	# simultaneously optimizes nearly-independent parameters
	
	phi <- (3 - sqrt(5))/2
	
	x <- w <- v <- b
	fx <- numeric(length(x))
	baseline <- f(x)
	fx[] <- baseline
	fw <- fv <- fx
	
	b <- ifelse(a > c, a, c)
	
	W <- (1:length(x))[(c - a) > tol]
	if (length(W)==0)
		return(x)
	
	e <- xm <- tol1 <- tol2 <- numeric(length(W))
	q <- r <- p <- numeric(length(W))
	
	iteration <- 0L
	maxIterations <- 1000L
	delta <- 0
	f_delta <- rep(Inf, length(W))
	while(iteration < maxIterations &&
		(iteration < 3 || delta > overallTol)) {
		iteration <- iteration + 1L
		
		xm[W] <- (a[W] + b[W])/2
		tol1[W] <- tol*abs(x[W]) + 1e-10
		tol2[W] <- 2*(tol1[W])
		
		w1 <- which(abs(x[W] - xm[W]) <= (tol2[W] - (b[W] - a[W])/2) |
			f_delta < relTol)
		if (length(w1) > 0) {
			W <- W[-w1]
			if (length(W)==0)
				break
		}
		
		w1 <- which(abs(e[W]) > tol1[W])
		w2 <- which(!(1:length(W) %in% w1))
		d <- numeric(length(W))
		
		if (length(w1) > 0) {
			r[W[w1]] <- (x[W[w1]] - w[W[w1]])*(fx[W[w1]] - fv[W[w1]])
			q[W[w1]] <- (x[W[w1]] - v[W[w1]])*(fx[W[w1]] - fw[W[w1]])
			p[W[w1]] <- (x[W[w1]] - v[W[w1]])*q[W[w1]] - (x[W[w1]] - w[W[w1]])*r[W[w1]]
			q[W[w1]] <- 2*(q[W[w1]] - r[W[w1]])
			p[W[w1]] <- ifelse(q[W[w1]] > 0,  p[W[w1]], -p[W[w1]])
			q[W[w1]] <- abs(q[W[w1]])
			
			etemp <- e[W[w1]]
			e[W[w1]] <- d[W[w1]]
			
			w3 <- which(abs(p[W[w1]]) >= abs(q[W[w1]]*etemp/2) |
				p[W[w1]] <= q[W[w1]]*(a[W[w1]] - x[W[w1]]) |
				p[W[w1]] >= q[W[w1]]*(b[W[w1]] - x[W[w1]]))
			w4 <- which(!(1:length(w1) %in% w3))
			if (length(w3) > 0) {
				e[W[w1[w3]]] <- ifelse(x[W[w1[w3]]] >= xm[W[w1[w3]]],
					a[W[w1[w3]]] - x[W[w1[w3]]],
					b[W[w1[w3]]] - x[W[w1[w3]]])
				d[W[w1[w3]]] <- phi*e[W[w1[w3]]]
			}
			if (length(w4) > 0) {
				d[W[w1[w4]]] <- p[W[w1[w4]]]/q[W[w1[w4]]]
				u <- x[W[w1[w4]]] + d[W[w1[w4]]]
				d[W[w1[w4]]] <- ifelse(u - a[W[w1[w4]]] < tol2[W[w1[w4]]] |
						b[W[w1[w4]]] - u < tol2[W[w1[w4]]],
					ifelse(xm[W[w1[w4]]] - x[W[w1[w4]]] >= 0,
						abs(tol1[W[w1[w4]]]),
						-abs(tol1[W[w1[w4]]])),
					d[W[w1[w4]]])
			}
		}
		if (length(w2) > 0) {
			e[W[w2]] <- ifelse(x[W[w2]] >= xm[W[w2]],
				a[W[w2]] - x[W[w2]],
				b[W[w2]] - x[W[w2]])
			d[W[w2]] <- phi*e[W[w2]]
		}
		
		# perform one function call per iteration
		u <- ifelse(abs(d[W]) >= tol1[W],
			x[W] + d[W],
			x[W] + ifelse(d[W] > 0,
				abs(tol1[W]),
				-abs(tol1[W])))
		fu <- f(x, W, u) # provide all alternative parameters
		
		newBaseline <- fu[1]
		fu <- fu[2:length(fu)]
		f_delta <- abs(fu - fx[W])
		offset <- fx[W] - newBaseline
		w1 <- which(offset > 0)
		if (length(w1) > 0) {
			fx[W[w1]] <- fx[W[w1]] - offset[w1]
			fw[W[w1]] <- fw[W[w1]] - offset[w1]
			fv[W[w1]] <- fv[W[w1]] - offset[w1]
		}
		delta <- abs(baseline - newBaseline)
		baseline <- newBaseline
		
		w1 <- which(fu <= fx[W])
		w2 <- which(!(1:length(W) %in% w1))
		if (length(w1) > 0){
			w3 <- which(u[w1] >= x[W[w1]])
			if (length(w3) > 0)
				a[W[w1[w3]]] <- x[W[w1[w3]]]
			w4 <- which(u[w1] < x[W[w1]])
			if (length(w4) > 0)
				b[W[w1[w4]]] <- x[W[w1[w4]]]
			v[W[w1]] <- w[W[w1]]
			w[W[w1]] <- x[W[w1]]
			x[W[w1]] <- u[w1]
			fv[W[w1]] <- fw[W[w1]]
			fw[W[w1]] <- fx[W[w1]]
			fx[W[w1]] <- fu[w1]
		}
		if (length(w2) > 0) {
			w3 <- which(u[w2] < x[W[w2]])
			if (length(w3) > 0)
				a[W[w2[w3]]] <- u[w2[w3]]
			w4 <- which(u[w2] >= x[W[w2]])
			if (length(w4) > 0)
				b[W[w2[w4]]] <- u[w2[w4]]
			
			w3 <- which(fu[w2] <= fw[W[w2]] | w[W[w2]]==x[W[w2]])
			if (length(w3) > 0) {
				v[W[w2[w3]]] <- w[W[w2[w3]]]
				w[W[w2[w3]]] <- u[w2[w3]]
				fv[W[w2[w3]]] <- fw[W[w2[w3]]]
				fw[W[w2[w3]]] <- fu[w2[w3]]
			}
			
			w4 <- which(fu[w2] <= fv[W[w2]] | v[W[w2]]==x[W[w2]] | v[W[w2]]==w[W[w2]])
			if (length(w4) > 0) {
				v[W[w2[w4]]] <- u[w2[w4]]
				fv[W[w2[w4]]] <- fu[w2[w4]]
			}
		}
	}
	
	return(x)
}

.reorderClusters <- function(myClusters, all=FALSE) {
	
	# order clusters by branching pattern
	repeat {
		a <- apply(as.matrix(myClusters[, 7:8], ncol=2), 1, max)
		w <- which(a > 1:dim(myClusters)[1])[1]
		if (is.na(w))
			break
		temp <- myClusters[w, c(4, 5, 7, 8)]
		myClusters[w:(a[w] - 1), c(4, 5, 7, 8)] <- myClusters[(w + 1):a[w], c(4, 5, 7, 8)]
		myClusters[a[w], c(4, 5, 7, 8)] <- temp
		w1 <- which(myClusters[w:dim(myClusters)[1], 7:8] %in% (w + 1):a[w])
		w2 <- which(myClusters[(a[w] + 1):dim(myClusters)[1], 7:8] %in% w)
		if (length(w1) > 0)
			myClusters[w:dim(myClusters)[1], 7:8][w1] <- myClusters[w:dim(myClusters)[1], 7:8][w1] - 1
		if (length(w2) > 0)
			myClusters[(a[w] + 1):dim(myClusters)[1], 7:8][w2] <- a[w]
	}
	
	if (all) { # also renumber columns 1 to 3
		count <- 0L
		myClusters[, 1:2] <- myClusters[, 7:8]
		for (i in 1:dim(myClusters)[1]) {
			if (myClusters[i, 7] > 0 && myClusters[i, 8] > 0) {
				myClusters[i, 1] <- myClusters[myClusters[i, 7], 3]
				myClusters[i, 2] <- myClusters[myClusters[i, 8], 3]
				count <- count + 1L
				myClusters[i, 3] <- count
			} else if (myClusters[i, 7] > 0) {
				myClusters[i, 1] <- myClusters[myClusters[i, 7], 3]
				myClusters[i, 3] <- myClusters[i, 1]
			} else if (myClusters[i, 8] > 0) {
				myClusters[i, 2] <- myClusters[myClusters[i, 8], 3]
				myClusters[i, 3] <- myClusters[i, 2]
			} else {
				count <- count + 1L
				myClusters[i, 3] <- count
			}
		}
	}
	
	return(myClusters)
}

.swapBranches <- function(myClusters, r1, c1, r2, c2) {
	
	# swap branch [r1, c1] with [r2, c2]
	
	temp <- myClusters[r1, c1]
	myClusters[r1, c1] <- myClusters[r2, c2]
	myClusters[r2, c2] <- temp
	temp <- myClusters[r1, c1 - 3]
	myClusters[r1, c1 - 3] <- myClusters[r2, c2 - 3]
	myClusters[r2, c2 - 3] <- temp
	
	myClusters <- .reorderClusters(myClusters)
	
	return(myClusters)
}

.NNI <- function(myClusters, bestLnL, NNIs, maximizeLikelihood, tol=1e-1) {
	
	if (dim(myClusters)[1]==1)
		return(myClusters)
	
	# perform rounds of nearest-neighbor interchanges
	W <- dim(myClusters)[1]:2
	while (length(W) > 0) {
		i <- W[1]
		W <- W[-1]
		
		if (i==dim(myClusters)[1]) { # swap nodes
			if (all(myClusters[i, 7:8] > 0)) { # two nodes
				myClustersTemp1 <- .swapBranches(myClusters,
					myClusters[i, 8], 7,
					myClusters[i, 7], 7)
				tempLnL1 <- maximizeLikelihood(myClustersTemp1, NNIs + 1, tol)
				myClustersTemp2 <- .swapBranches(myClusters,
					myClusters[i, 8], 7,
					myClusters[i, 7], 8)
				tempLnL2 <- maximizeLikelihood(myClustersTemp2, NNIs + 1, tol)
				if (tempLnL1 < bestLnL - tol &&
					tempLnL1 < tempLnL2) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp1
					bestLnL <- tempLnL1
				} else if (tempLnL2 < bestLnL - tol) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp2
					bestLnL <- tempLnL2
				}
			}
		} else {
			w <- which(myClusters[i, 7:8] > 0)
			if (length(w)==2) {
				myClustersTemp1 <- .swapBranches(myClusters,
					i, 8,
					myClusters[i, 7], 7)
				tempLnL1 <- maximizeLikelihood(myClustersTemp1, NNIs + 1, tol)
				myClustersTemp2 <- .swapBranches(myClusters,
					i, 8,
					myClusters[i, 7], 8)
				tempLnL2 <- maximizeLikelihood(myClustersTemp2, NNIs + 1, tol)
				myClustersTemp3 <- .swapBranches(myClusters,
					i, 7,
					myClusters[i, 8], 7)
				tempLnL3 <- maximizeLikelihood(myClustersTemp3, NNIs + 1, tol)
				myClustersTemp4 <- .swapBranches(myClusters,
					i, 7,
					myClusters[i, 8], 8)
				tempLnL4 <- maximizeLikelihood(myClustersTemp4, NNIs + 1, tol)
				if (tempLnL1 < bestLnL - tol &&
					tempLnL1 < tempLnL2 &&
					tempLnL1 < tempLnL3 &&
					tempLnL1 < tempLnL4) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp1
					bestLnL <- tempLnL1
				} else if (tempLnL2 < bestLnL - tol &&
					tempLnL2 < tempLnL3 &&
					tempLnL2 < tempLnL4) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp2
					bestLnL <- tempLnL2
				} else if (tempLnL3 < bestLnL - tol &&
					tempLnL3 < tempLnL4) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp3
					bestLnL <- tempLnL3
				} else if (tempLnL4 < bestLnL - tol) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp4
					bestLnL <- tempLnL4
				}
			} else if (length(w)==1) {
				myClustersTemp1 <- .swapBranches(myClusters,
					i, ifelse(w==1, 8, 7),
					myClusters[i, w + 6], 7)
				tempLnL1 <- maximizeLikelihood(myClustersTemp1, NNIs + 1, tol)
				myClustersTemp2 <- .swapBranches(myClusters,
					i, ifelse(w==1, 8, 7),
					myClusters[i, w + 6], 8)
				tempLnL2 <- maximizeLikelihood(myClustersTemp2, NNIs + 1, tol)
				if (tempLnL1 < bestLnL - tol &&
					tempLnL1 < tempLnL2) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp1
					bestLnL <- tempLnL1
				} else if (tempLnL2 < bestLnL - tol) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp2
					bestLnL <- tempLnL2
				}
			}
		}
	}
	
	# return the highest likelihood tree
	return(myClusters)
}

MODELS <- c("JC69",
	"JC69+G4",
	"K80",
	"K80+G4",
	"F81",
	"F81+G4",
	"HKY85",
	"HKY85+G4",
	"T92",
	"T92+G4",
	"TN93",
	"TN93+G4")

.splitClusters <- function(x, y) {
	clusterNum <- 0L
	X <- integer(length(x))
	u.y <- unique(y)
	for (i in u.y) {
		w.y <- which(y==i)
		u.x <- unique(x[w.y])
		for (j in u.x) {
			clusterNum <- clusterNum + 1L
			w.x <- which(x[w.y]==j)
			X[w.y[w.x]] <- clusterNum
		}
	}
	return(X)
}

.midpointRoot <- function(dendrogram) {
	
	# mid-point root the tree based on which
	# leaf is furthest from the leaf at zero
	
	# find the tip that is at zero height
	.zeroFound <- FALSE
	.containsZero <- function(x) {
		if (is.leaf(x)) {
			if (!.zeroFound && isTRUE(all.equal(attr(x, "height"), 0))) {
				.zeroFound <<- TRUE
				attr(x, "containsZero") <- TRUE
			} else {
				attr(x, "containsZero") <- FALSE
			}
		} else {
			x[[1]] <- .containsZero(x[[1]])
			x[[2]] <- .containsZero(x[[2]])
			if (attr(x[[1]], "containsZero") || attr(x[[2]], "containsZero")) {
				attr(x, "containsZero") <- TRUE
			} else {
				attr(x, "containsZero") <- FALSE
			}
		}
		
		return(x)
	}
	
	dendrogram <- .containsZero(dendrogram)
	
	# find the maximum distance to the zeroth tip
	.maxDist <- 0
	.maxLeaf <- NA
	.findMax <- function(x, height=NA) {
		if (is.leaf(x)) {
			dist <- 2*height - attr(x, "height")
			if (!is.na(dist) && dist >= .maxDist) {
				.maxDist <<- dist
				.maxLeaf <<- as.integer(x)
			}
		} else {
			if (is.na(height)) {
				if (attr(x[[1]], "containsZero")) {
					x[[1]] <- .findMax(x[[1]])
				} else {
					x[[1]] <- .findMax(x[[1]], attr(x, "height"))
				}
				if (attr(x[[2]], "containsZero")) {
					x[[2]] <- .findMax(x[[2]])
				} else {
					x[[2]] <- .findMax(x[[2]], attr(x, "height"))
				}
			} else {
				x[[1]] <- .findMax(x[[1]], height)
				x[[2]] <- .findMax(x[[2]], height)
			}
		}
		
		return(x)
	}
	
	dendrogram <- .findMax(dendrogram)
	
	# find paths that contain the most distant tips
	.containsMax <- function(x) {
		if (is.leaf(x)) {
			if (x==.maxLeaf) {
				attr(x, "containsMax") <- TRUE
			} else {
				attr(x, "containsMax") <- FALSE
			}
		} else {
			x[[1]] <- .containsMax(x[[1]])
			x[[2]] <- .containsMax(x[[2]])
			if (attr(x[[1]], "containsMax") || attr(x[[2]], "containsMax")) {
				attr(x, "containsMax") <- TRUE
			} else {
				attr(x, "containsMax") <- FALSE
			}
		}
		
		return(x)
	}
	
	dendrogram <- .containsMax(dendrogram)
	
	
	# midpoint root the tree
	.findMidpoint <- function(x) {
		if (is.leaf(x)) {
			if (attr(x, "containsZero") &&
				attr(x, "height") >= .maxDist/2) {
				attr(x, "containsRoot") <- TRUE
				attr(x, "isRoot") <- TRUE
			} else {
				attr(x, "containsRoot") <- FALSE
			}
		} else if (attr(x, "height") >= .maxDist/2 &&
			attr(x, "containsZero")) {
			x[[1]] <- .findMidpoint(x[[1]])
			x[[2]] <- .findMidpoint(x[[2]])
			if (attr(x[[1]], "height") < .maxDist/2 &&
				attr(x[[1]], "containsZero") &&
				!attr(x[[1]], "containsMax")) {
				attr(x, "containsRoot") <- TRUE
				attr(x, "isRoot") <- TRUE
			} else if (attr(x[[2]], "height") < .maxDist/2 &&
				attr(x[[2]], "containsZero") &&
				!attr(x[[2]], "containsMax")) {
				attr(x, "containsRoot") <- TRUE
				attr(x, "isRoot") <- TRUE
			} else {
				attr(x, "containsRoot") <- attr(x[[1]], "containsRoot") ||
				attr(x[[2]], "containsRoot")
			}
		} else {
			attr(x, "containsRoot") <- FALSE
		}
		
		return(x)
	}
	
	dendrogram <- .findMidpoint(dendrogram)
	
	lower <- list()
	.adjustHeight <- function(x, diff=NA) {
		if (!is.null(attr(x, "containsRoot")) &&
			attr(x, "containsRoot"))
			diff <- 2*attr(x, "height") - .maxDist
		if (!is.na(diff)) {
			attr(x, "height") <- attr(x, "height") - diff
			if (!is.leaf(x)) {
				if (is.null(attr(x, "isRoot"))) {
					x[[1]] <- .adjustHeight(x[[1]], diff)
					if (length(x) > 1) {
						x[[2]] <- .adjustHeight(x[[2]], diff)
					} else { # first child was removed
						x[[1]] <- .adjustHeight(x[[1]], diff)
					}
				} else if (attr(x[[1]], "containsZero")) {
					x[[2]] <- .adjustHeight(x[[2]], diff)
				} else {
					x[[1]] <- .adjustHeight(x[[1]], diff)
				}
			}
		}
		if (!is.null(attr(x, "isRoot"))) {
			if (is.leaf(x)) {
				lower <<- x
				x <- NULL
			} else if (attr(x[[1]], "containsZero")) {
				lower <<- x[[1]]
				x[[1]] <- NULL
			} else {
				lower <<- x[[2]]
				x[[2]] <- NULL
			}
		}
		return(x)
	}
	
	upper <- .adjustHeight(dendrogram)
	
	# initalize a tree from the lower half
	branch <- 0L
	tree <- list(lower, list())
	attr(tree, "height") <- .maxDist/2
	attr(tree, "class") <- "dendrogram"
	total <- attr(dendrogram, "members")
	attr(tree, "members") <- total
	total <- total - attr(lower, "members")
	# reverse edges between original and new root
	.revertList <- function(x) {
		temp <- NULL
		if (!is.null(attr(x[[1]], "containsRoot")) &&
			attr(x[[1]], "containsRoot")) {
			.revertList(x[[1]]) # continue descending
			temp <- x[[2]]
		} else if (length(x) > 1 &&
			!is.null(attr(x[[2]], "containsRoot")) &&
			attr(x[[2]], "containsRoot")) {
			.revertList(x[[2]]) # continue descending
			temp <- x[[1]]
		} else if (!is.null(attr(x, "isRoot")) &&
			attr(x, "isRoot")) {
			temp <- x # last node
			attr(temp, "members") <- attr(temp, "members") - attr(lower, "members")
		}
		if (!is.null(temp)) {
			if (!is.leaf(temp) && length(temp)==1)
				temp <- temp[[1L]] # remove node
			branch <<- branch + 1L # go to next branch
			index <- rep(2L, branch)
			# initialize the subtree
			tree[[index]] <<- list()
			class(tree[[index]]) <<- "dendrogram"
			attr(tree[[index]], "height") <<- attr(x, "height")
			attr(tree[[index]], "members") <<- total
			if (is.leaf(temp)) {
				total <<- total - 1L
			} else {
				total <<- total - attr(temp, "members")
			}
			if (total > 0) {
				tree[[c(index, 1L)]] <<- temp
			} else { # last branch
				tree[[index]] <<- temp
			}
		}
		invisible(NULL)
	}
	
	.revertList(upper)
	if (length(tree[[2]])==0)
		tree[[2]] <- upper[[1]]
	
	# correct midpoints on new tree
	.adjustMidpoints <- function(x) {
		if (is.leaf(x[[1]]) && is.leaf(x[[2]])) {
			attr(x, "midpoint") <- 0.5
		} else if (is.leaf(x[[1]])) {
			x[[2]] <- .adjustMidpoints(x[[2]])
			attr(x, "midpoint") <- (1L + attr(x[[2]], "midpoint"))/2
		} else if (is.leaf(x[[2]])) {
			x[[1]] <- .adjustMidpoints(x[[1]])
			attr(x, "midpoint") <- (attr(x[[1]], "midpoint") + attr(x[[1]], "members"))/2
		} else {
			x[[1]] <- .adjustMidpoints(x[[1]])
			x[[2]] <- .adjustMidpoints(x[[2]])
			attr(x, "midpoint") <- (attr(x[[1]], "members") + attr(x[[1]], "midpoint") + attr(x[[2]], "midpoint"))/2
		}
		
		if (!is.null(attr(x, "isRoot")))
			attr(x, "isRoot") <- NULL
		if (!is.null(attr(x, "containsRoot")))
			attr(x, "containsRoot") <- NULL
		if (!is.null(attr(x, "containsZero")))
			attr(x, "containsZero") <- NULL
		if (!is.null(attr(x, "containsMax")))
			attr(x, "containsMax") <- NULL
		
		return(x)
	}
	
	rooted <- .adjustMidpoints(tree)
	
	return(rooted)
}

IdClusters <- function(myDistMatrix=NULL,
	method="UPGMA",
	cutoff=-Inf,
	showPlot=FALSE,
	type="clusters",
	myXStringSet=NULL,
	model=MODELS,
	processors=1,
	verbose=TRUE) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	METHODS <- c("NJ","UPGMA", "ML", "complete", "single", "WPGMA", "inexact")
	method <- pmatch(method, METHODS)
	if (is.na(method) && !is.null(myDistMatrix)) {
		stop("Invalid method.  Choose either ML, NJ, complete, single, WPGMA, or UPGMA.")
	} else if (is.na(method)) {
		stop("Invalid method.")
	}
	if (method==-1 && !is.null(myDistMatrix)) {
		stop("Ambiguous method.  Choose either ML, NJ, complete, single, WPGMA, or UPGMA.")
	} else if (method==-1) {
		stop("Ambiguous method.")
	}
	if (method==3 && length(model) < 1)
		stop("No model(s) specified.")
	if (method==3 && !is.character(model))
		stop("model must be a character vector.")
	w <- which(!(model %in% MODELS))
	if (method==3 && length(w) > 0) {
		submodels <- sub("([^+]*)(\\+G(\\d+))?", "\\1", model[w])
		if (!all(submodels %in% sub("\\+G(\\d+)$", "", MODELS)))
			stop(paste("Available models are:",
				paste(MODELS, collapse=", "),
				collapse=" "))
		rates <- as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", model[w]))
		if (any(is.na(rates)) || any(floor(rates)!=rates))
			stop("The number rates in the discrete Gamma distribution (i.e., +G4) should be an integer value.")
		if (any(rates > 10))
			stop("Up to 10 rates are allowed in the discrete Gamma distribution (i.e., +G10).")
		if (any(rates < 2))
			stop("A minimum of two rates are required for the discrete Gamma distribution (i.e., +G2).")
	}
	model <- unique(model)
	if (!is.numeric(cutoff))
		stop("cutoff must be a numeric.")
	if (is.integer(cutoff))
		cutoff <- as.numeric(cutoff)
	if (!is.logical(showPlot))
		stop("showPlot must be a logical.")
	TYPES <- c("clusters", "dendrogram", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (method==7 && showPlot)
		stop("showPlot must be FALSE if method is 'inexact'")
	if (method==7 && type > 1)
		stop("type must be 'clusters' when method is 'inexact'")
	if (length(cutoff) > 1 && type > 1)
		warning("More than one cutoff specified when type is ", TYPES[type], ".")
	if (method==7 && any(cutoff < 0))
		stop("cutoff must be at least zero when method is 'inexact'.")
	if (method==7 && any(cutoff >= 1))
		stop("cutoff must be less than one when method is 'inexact'.")
	ASC <- TRUE
	if (is.unsorted(cutoff)) {
		if (is.unsorted(rev(cutoff))) {
			stop("cutoff must be sorted.")
		} else {
			ASC <- FALSE
		}
	}
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
	if (method != 7 && !is(myDistMatrix, "matrix")) {
		stop(paste("myDistMatrix must be a matrix for method '", METHODS[method], "'.", sep=""))
	} else if (method != 7 && method != 8) {
		dim <- dim(myDistMatrix)
		if (dim[2]!=dim[1])
			stop("myDistMatrix is not square.")
		dim <- dim[1]
		if (dim < 2)
			stop("myDistMatrix is too small.")
		if (typeof(myDistMatrix)=="integer")
			myDistMatrix[] <- as.numeric(myDistMatrix)
	}
	if (method == 3) {
		if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
			stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
		if (length(myXStringSet)!=dim(myDistMatrix)[1])
			stop("myDistMatrix must have as many rows as the number of sequences.")
		if (length(unique(width(myXStringSet))) != 1)
			stop("All sequences in myXStringSet must be the same width (aligned).")
		if (!is.null(attr(myDistMatrix, "correction")) &&
			attr(myDistMatrix, "correction") == "none")
			stop('myDistMatrix must have a correction with method="ML".')
	}
	
	if (method == 7) {
		type <- switch(class(myXStringSet),
			`DNAStringSet` = 1L,
			`RNAStringSet` = 2L,
			`AAStringSet` = 3L,
			stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet."))
		a <- vcountPattern("-", myXStringSet)
		if (any(a > 0))
			stop("Gap characters ('-') must be removed before inexact clustering.")
		a <- vcountPattern("+", myXStringSet)
		if (any(a > 0))
			stop("Mask characters ('+') must be removed before inexact clustering.")
		a <- vcountPattern(".", myXStringSet)
		if (any(a > 0))
			stop("Unknown characters ('.') must be removed before inexact clustering.")
		if (all(width(myXStringSet)==0L))
			stop("All sequences in myXStringSet are zero width.")
		
		if (verbose) {
			lastValue <- 0
			pBar <- txtProgressBar(style=3)
		}
		if (type==3L) { # AAStringSet
			wordSize <- ceiling(log(100*mean(width(myXStringSet)), 20))
			if (wordSize > 7)
				wordSize <- 7
			if (wordSize < 1)
				wordSize <- 1
			words <- 20^wordSize
		} else { # DNAStringSet or RNAStringSet
			wordSize <- ceiling(log(100*mean(width(myXStringSet)), 4))
			if (wordSize > 15)
				wordSize <- 15
			if (wordSize < 1)
				wordSize <- 1
			words <- 4^wordSize
		}
		
		l <- length(myXStringSet)
		if (l==0)
			stop("myXStringSet contains no sequences.")
		lc <- length(cutoff)
		if (is.null(names(myXStringSet))) {
			dNames <- 1:l
		} else {
			dNames <- names(myXStringSet)
			w <- which(duplicated(dNames))
			if (length(w) > 0) {
				warning("Duplicated names of myXStringSet appended with index.")
				dNames[w] <- paste(dNames[w],
					w,
					sep="_")
			}
		}
		if (lc > 1) {
			cNames <- paste("cluster",
				gsub("\\.", "_", cutoff),
				METHODS[method],
				sep="")
		} else {
			cNames <- "cluster"
		}
		c <- matrix(0L,
			nrow=l,
			ncol=lc,
			dimnames=list(dNames,
				cNames))
		
		# identify duplicated sequences
		x <- selfmatch(myXStringSet)
		u <- unique(x)
		l <- length(u)
		t <- tabulate(x, length(x))[u]
		
		cutoff <- cutoff/2 # cluster radius is half the diameter
		for (i in seq_len(lc)) {
			if (!ASC && i > 1) {
				o <- u[order(c[u, i - 1],
					width(myXStringSet)[u],
					t, # frequency
					decreasing=TRUE)]
			} else {
				o <- u[order(width(myXStringSet)[u],
					t, # frequency
					decreasing=TRUE)]
			}
			
			if (type==3L) { # AAStringSet
				v <- .Call("enumerateSequenceAA",
					.subset(myXStringSet, o[1:ifelse(l > 999, 999, l)]),
					wordSize,
					PACKAGE="DECIPHER")
			} else { # DNAStringSet or RNAStringSet
				v <- .Call("enumerateSequence",
					.subset(myXStringSet, o[1:ifelse(l > 999, 999, l)]),
					wordSize,
					PACKAGE="DECIPHER")
			}
			v <- lapply(v,
				sort.int,
				method="radix")
			
			cluster_num <- 1L
			offset <- 0L
			c[o[cluster_num], i] <- cluster_num
			
			nGroups <- 1L
			seeds.reps <- v[1L]
			seeds.nums <- list(v[1L])
			seeds.seqs <- list(.subset(myXStringSet, o[1]))
			seeds.clust <- list(1L)
			
			index <- 1L
			for (j in seq_along(o)[-1]) {
				if (verbose) {
					value <- round(((lc - i)*l + j)/lc/l, 2)
					if (value > lastValue) {
						lastValue <- value
						setTxtProgressBar(pBar, value)
					}
				}
				
				if (j %% 1000 == 0) {
					index <- 1L
					if (type==3L) { # AAStringSet
						v <- .Call("enumerateSequenceAA",
							.subset(myXStringSet, o[j:ifelse(j + 999 > l, l, j + 999)]),
							wordSize,
							PACKAGE="DECIPHER")
					} else { # DNAStringSet or RNAStringSet
						v <- .Call("enumerateSequence",
							.subset(myXStringSet, o[j:ifelse(j + 999 > l, l, j + 999)]),
							wordSize,
							PACKAGE="DECIPHER")
					}
					v <- lapply(v,
						sort.int,
						method="radix")
				} else {
					index <- index + 1L
				}
				
				if (!ASC && i > 1) {
					if (c[o[j], i - 1] != c[o[j - 1], i - 1]) {
						# different clusters in last cutoff
						offset <- offset + cluster_num
						cluster_num <- 1L
						c[o[j], i] <- cluster_num + offset
						nGroups <- 1L
						seeds.reps <- v[index]
						seeds.nums <- list(v[index])
						seeds.seqs <- list(.subset(myXStringSet, o[j]))
						seeds.clust <- list(cluster_num)
						next
					}
				}
				
				# Determine to which group does the sequence belong
				m <- .Call("matchListsDual",
					v[index],
					seeds.reps,
					FALSE, # verbose
					NULL, # pBar
					processors,
					PACKAGE="DECIPHER")
				
				# expected number of occurrences of a random word
				probs <- lengths(seeds.reps)/words
				probs <- ifelse(probs > 1, 1, probs)
				# probability of >= `q` matches in `size` trials
				prob <- 1 - pbinom(q=m*length(v[[index]]),
					size=length(v[[index]]),
					prob=probs)
				group <- which(prob < 0.01) # < 1% likelihood by chance
				
				if (length(group)==0) {
					# form a new group
					cluster_num <- cluster_num + 1L
					c[o[j], i] <- cluster_num
					nGroups <- nGroups + 1L
					seeds.reps[nGroups] <- v[index]
					seeds.nums[[nGroups]] <- v[index]
					seeds.seqs[[nGroups]] <- .subset(myXStringSet, o[j])
					seeds.clust[[nGroups]] <- cluster_num
				} else { # part of an existing group
					if (length(group) > 1)
						group <- group[which.max(m[group])]
					
					m <- .Call("matchListsDual",
						v[index],
						seeds.nums[[group]],
						FALSE, # verbose
						NULL, # pBar
						processors,
						PACKAGE="DECIPHER")
					
					w <- which((1 - m) <= cutoff[i])
					if (length(w) > 0) {
						w <- w[which.max(m[w])]
						c[o[j], i] <- seeds.clust[[group]][w] + offset
					} else {
						pattern <- .subset(myXStringSet, o[j])
						
						weights <- m^3 # emphasize closest sequences
						weights <- weights/mean(weights)
						if (any(!is.finite(weights)))
							weights <- 1
						
						temp <- AlignProfiles(pattern,
							seeds.seqs[[group]],
							s.weight=weights,
							processors=processors)
						d <- .Call("distMatrix",
							temp,
							type,
							FALSE, # includeTerminalGaps
							FALSE, # penalizeGapGapMatches
							TRUE, # penalizeGapLetterMatches
							FALSE, # full matrix
							FALSE, # verbose
							NULL, # pBar
							processors,
							PACKAGE="DECIPHER")[-1L]
						w <- which(d <= cutoff[i])
						if (length(w) > 0) {
							w <- w[which.min(d[w])]
							c[o[j], i] <- seeds.clust[[group]][w] + offset
						} else { # form a new cluster
							cluster_num <- cluster_num + 1L
							c[o[j], i] <- cluster_num + offset
							seeds.nums[[group]] <- c(v[index],
								seeds.nums[[group]])
							seeds.seqs[[group]] <- temp
							seeds.clust[[group]] <- c(cluster_num,
								seeds.clust[[group]])
						}
					}
				}
			}
			
			c[, i] <- c[x, i]
		}
		myClusters <- as.data.frame(c)
	} else {
		w1 <- which(is.infinite(myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)]))
		if (length(w1) > 0) {
			if (verbose)
				warning("myDistMatrix contains infinite values.\n",
					"Replaced infinite values with max distance >= 1.\n")
			myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)][w1] <- NA
		}
		
		w2 <- which(is.na(myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)]))
		if (length(w2) > 0) {
			if (verbose &&
				length(w2) > length(w1))
				warning("myDistMatrix contains NA values.\n",
					"Replaced NA values with max distance >= 1.\n")
			max.dist <- max(myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)], na.rm=TRUE)
			if (max.dist <= 1) {
				myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)][w2] <- 1
			} else {
				myDistMatrix[lower.tri(myDistMatrix, diag=FALSE)][w2] <- max.dist
			}
		}
		
		# initialize a progress bar
		if (verbose) {
			if (method==3) {
				cat("Constructing initial neighbor-joining tree:\n")
				flush.console()
			}
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
				processors,
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
				processors,
				PACKAGE="DECIPHER")
		
		if (method==3) {
			# create NJ for use as initial tree
			myClusters <- .Call("clusterNJ",
				myDistMatrix,
				cutoff=-Inf,
				verbose,
				pBar,
				processors,
				PACKAGE="DECIPHER")
			myClusters <- .reorderClusters(myClusters)
			
			m <- matrix(NA,
				nrow=length(model),
				ncol=10,
				dimnames=list(model,
					c("FreqA", "FreqC", "FreqG", "FreqT",
						"A2G", "C2T", "alpha",
						"-LnL", "AICc", "BIC")))
			
			if (!(length(model)==1 && model[1]=="JC69")) {
				N <- sum(ifelse(apply(consensusMatrix(myXStringSet),
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
						myXStringSet,
						N,
						processors=processors)
					if (verbose)
						cat("\n", model[i],
							":", paste(rep(" ",
									max(nchar(rownames(m))) - nchar(model[i]) + 1),
								collapse=""),
							"-ln(L)=", round(m[model[i], "-LnL"], 0),
							", AICc=", round(m[model[i], "AICc"], 0),
							", BIC=", round(m[model[i], "BIC"], 0),
							sep="")
				}
			}
			
			if (length(model) > 1) { # choose the best model
				w <- which.min(m[,"BIC"])
				l <- logical(nrow(m))
				l[w] <- TRUE
				m <- subset(m, l)
			}
			
			if (verbose) {
				if (length(model) > 1)
					cat("\n\nThe selected model was:  ",
						rownames(m),
						sep="")
			}
			
			.giveParams <- function(model_params) {
				w <- which(is.na(model_params))
				if (length(w) > 0) {
					model_params[w] <- c(0.25, 0.25, 0.25, 0.25, 1, 1, NA)[w]
				}
				if (is.na(model_params[7])) {
					model_params <- c(model_params[1:6], 1, 1)
				} else {
					model_params <- c(model_params[1:6], .rates(model_params[7], as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", rownames(m)))))
				}
			}
			model_params <- .giveParams(as.numeric(m[1:7]))
			
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
			.NNIs <- 0L
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
						"% improvement), ",
						.NNIs,
						ifelse(.NNIs==1, " NNI  ", " NNIs "),
						sep="")
					flush.console()
					invisible(value)
				}
			}
			
			# given branch lengths return -LnL
			maximizeLikelihood <- function(x, branches=integer(), lengths=numeric()) {
				myClusters[, 4:5] <- x
				LnL <- .Call("clusterML",
					myClusters,
					myXStringSet,
					model_params,
					branches,
					lengths,
					processors,
					PACKAGE="DECIPHER")
				
				w <- which.min(LnL)
				if (LnL[w] < .bestLnL) {
					if (w > 1) {
						x[branches[w - 1]] <- lengths[w - 1]
						params <<- x
					} else {
						params <<- x
					}
					.bestLnL <<- LnL[w]
					if (verbose) {
						if (is.infinite(.startLnL))
							.startLnL <<- LnL[1]
						printLine(LnL[w])
					}
				}
				
				return(LnL)
			}
			
			maximizeLikelihood2 <- function(myClusters, NNIs, tol) {
				LnL <- .Call("clusterML",
					myClusters,
					myXStringSet,
					model_params,
					integer(),
					numeric(),
					processors,
					PACKAGE="DECIPHER")
				
				if (LnL < .bestLnL - tol) {
					.bestLnL <<- LnL
					if (verbose) {
						.NNIs <<- NNIs
						if (is.infinite(.startLnL))
							.startLnL <<- LnL
						printLine(LnL)
					}
				}
				
				return(LnL)
			}
			
			# maximize likelihood of tree
			index <- rep(TRUE, 2*dim(myClusters)[1])
			repeat {
				currentLnL <- .bestLnL
				currentNNIs <- .NNIs
				params <- as.numeric(myClusters[, 4:5])
				tempClusters <- paste(myClusters[, 7], myClusters[, 8])
				.simultaneousBrent(maximizeLikelihood,
					ifelse(index, 0, params),
					params,
					ifelse(index, 10*params, params))
				myClusters[, 4:5] <- params
				
				myClusters <- .NNI(myClusters,
					.bestLnL,
					.NNIs,
					maximizeLikelihood2)
				
				if (abs(.bestLnL - currentLnL) < 1e0 &&
					currentNNIs==.NNIs) {
					temp_params <- model_params
				} else {
					m[1,] <- .optimizeModel(myClusters,
						rownames(m),
						myXStringSet,
						N,
						processors=processors)
					temp_params <- .giveParams(as.numeric(m[1:7]))
				}
				
				if ((abs(.bestLnL - currentLnL) < 1e0 || # negligible improvement in likelihood
					currentNNIs==.NNIs) && # no new nearest neighbor interchanges (NNIs)
					all(abs(temp_params - model_params) < 0.01)) # negligible change in model_params
					break
				model_params <- temp_params
				index <- rep(!(paste(myClusters[, 7], myClusters[, 8]) %in% tempClusters), 2)
			}
			
			myClusters <- .reorderClusters(myClusters, all=TRUE)
			myClusters <- adjustTreeHeights(myClusters)
			myClusters <- .Call("reclusterNJ",
				myClusters,
				cutoff[1],
				PACKAGE="DECIPHER")
			
			if (verbose) {
				printLine(.bestLnL, 100)
				flush.console()
				params <- formatC(round(m, 3),
					digits=3,
					format="f")
				cat(ifelse(any(grepl("NA", params[1], fixed=TRUE)),
						"",
						"\n\nModel parameters:"),
					ifelse(grepl("NA", params[1], fixed=TRUE),
						"",
						ifelse(grepl("T92", rownames(m), fixed=TRUE),
							paste("\nFrequency(A) = Frequency(T) = ", params[1],
								"\nFrequency(C) = Frequency(G) = ", params[2],
								sep=""),
							paste("\nFrequency(A) = ", params[1],
								"\nFrequency(C) = ", params[2],
								"\nFrequency(G) = ", params[3],
								"\nFrequency(T) = ", params[4],
								sep=""))),
					ifelse(grepl("NA", params[5], fixed=TRUE),
						"",
						ifelse(grepl("TN93", rownames(m), fixed=TRUE),
							paste("\nRate A <-> G = ", params[5],
								"\nRate C <-> T = ", params[6],
								"\nTransversion rates = 1",
								sep=""),
							paste("\nTransition rates = ", params[5],
								"\nTransversion rates = 1",
								sep=""))),
					ifelse(grepl("NA", params[7], fixed=TRUE),
						"",
						paste("\nAlpha = ", params[7], sep="")),
					sep="")
				cat("\n")
			}
		}
		
		if (showPlot || type > 1) {
			# create a dendrogram object
			myClustersList <- list()
			if (is.null(dimnames(myDistMatrix)[[1]])) {
				myClustersList$labels <- 1:(dim(myClusters)[1] + 1)
			} else {
				myClustersList$labels <- dimnames(myDistMatrix)[[1]]
				
				w <- which(duplicated(dimnames(myDistMatrix)[[1]]))
				if (length(w) > 0) {
					warning("Duplicated dimnames in dendrogram appended with index.")
					myClustersList$labels[w] <- paste(myClustersList$labels[w],
						w,
						sep="_")
				}
			}
			
			myClustersList$merge <- matrix(myClusters[,7:8], ncol=2)
			myClustersList$height <- matrix(myClusters[,6], ncol=1)
			myClustersList$lengths <- matrix(myClusters[,4:5], ncol=2)
			myClustersList$clusters <- matrix(myClusters[,9:10], ncol=2)
			if (dim > 100) {
				fontSize <- .6
			} else if (dim > 70) {
				fontSize <- .7
			} else if (dim > 40) {
				fontSize <- .8
			} else {
				fontSize <- .9
			}
			if (dim > 300) {
				leaves <- "none"
			} else {
				leaves <- "perpendicular"
			}
			
			d <- to.dendrogram(myClustersList)
			
			# need to maximize the recursion depth temporarily
			org.options <- options(expressions=5e5)
			on.exit(options(org.options))
			
			# convert bifurcating tree to multifurcating
			d <- .collapse(d)
			
			# midpoint root the dendrogram
			if (method==1 || method==3)
				d <- .midpointRoot(d)
			
			# specify the order of clusters that
			# will match the plotted dendrogram
			orderDendrogram <- order.dendrogram(d)
			
			c <- .organizeClusters(myClusters, myClustersList$labels, orderDendrogram)
			
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
			if (is.finite(cutoff[1])) {
				d <- dendrapply(d, colEdge, c, r)
				
				.reorder <- function(dend) {
					l <- length(dend)
					if (l > 1) {
						for (i in seq_len(l))
							dend[[i]] <- .reorder(dend[[i]])
						
						members <- lapply(dend, unlist)
						# sort tree by ascending cluster number
						o <- sort.list(sapply(members,
								function(x)
									min(c[x, 1])))
						dend[] <- dend[o]
					} else if (!is.leaf(dend)) {
						dend[[1]] <- .reorder(dend[[1]])
					}
					return(dend)
				}
				d <- .reorder(d)
			}
		}
		if (type==1 || type==3) {
			if (is.null(dimnames(myDistMatrix)[[1]])) {
				dNames <- 1:(dim(myClusters)[1] + 1)
			} else {
				dNames <- dimnames(myDistMatrix)[[1]]
				
				w <- which(duplicated(dNames))
				if (length(w) > 0) {
					warning("Duplicated dimnames in myDistMatrix appended with index.")
					dNames[w] <- paste(dNames[w],
						w,
						sep="_")
				}
			}
			
			if (type==1) # do not number clusters by order of appearance
				c <- .organizeClustersFast(myClusters, dNames)
			
			if (length(cutoff) > 1) {
				names(c) <- paste("cluster",
					gsub("\\.", "_", cutoff[1]),
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
					if ((method==1 || method==3) && !ASC) # ensure clusters are subsets
						x[, 1] <- .splitClusters(x[, 1], c[, dim(c)[2]])
					names(x) <- paste("cluster",
						gsub("\\.", "_", cutoff[i]),
						METHODS[method],
						sep="")
					c <- cbind(c, x)
				}
			}
			myClusters <- c
		}
		
		if (showPlot)
			plot(d,
				horiz=FALSE,
				leaflab=leaves,	
				nodePar=list(lab.cex=fontSize, pch = NA))
	}
	
	if (verbose) {
		if (method != 3) { # already closed pBar
			setTxtProgressBar(pBar, 100)
			close(pBar)
		}
		
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	if (type==1) {
		return(myClusters)
	} else if (type==2) {
		return(d)
	} else {
		return(list(myClusters, d))
	}
}
