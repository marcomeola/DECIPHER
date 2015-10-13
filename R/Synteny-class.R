`[.Synteny` <- function(x, i, j, ...) {
	ans <- NextMethod("[", x)
	
	if (missing(j))
		return(ans)
	
	d <- dim(x)
	I <- seq_len(d[1])
	J <- seq_len(d[1])
	d <- dimnames(x)
	names(I) <- d[[1]]
	names(J) <- d[[1]]
	I <- I[i]
	J <- J[j]
	if (length(I) >= 2 &&
		length(I)==length(J) &&
		all(I==J) &&
		!any(duplicated(I)))
		class(ans) <- "Synteny"
	return(ans)
}

print.Synteny <- function(x,
	quote=FALSE,
	right=TRUE,
	...) {
	
	d <- dim(x)
	if (is.null(d))
		stop("x must be a square object of class 'Synteny'.")
	m <- matrix("",
		nrow=d[1],
		ncol=d[2],
		dimnames=dimnames(x))
	l <- s <- integer(d[1])
	for (i in seq_len(d[1])) {
		l[i] <- length(x[i, i][[1]])
		s[i] <- sum(x[i, i][[1]])
	}
	for (i in seq_len(d[1])) {
		for (j in seq_len(d[2])) {
			if (i==j) {
				m[i, j] <- paste(l[i],
					ifelse(l[i]==1,
						"seq",
						"seqs"))
			} else if (i > j) {
				t <- dim(x[i, j][[1]])[1]
				m[i, j] <- paste(t,
					ifelse(t==1,
						"block",
						"blocks"))
			} else { # i < j
				t <- sum(x[i, j][[1]][, "width"])
				m[i, j] <- paste(format(ifelse(s[i] < s[j],
							100*t/s[i],
							100*t/s[j]),
						digits=2),
					"% hits",
					sep="")
			}
		}
	}
	
	print(m,
		quote=quote,
		right=right,
		...)
}

pairs.Synteny <- function(x,
	bounds=TRUE,
	boxBlocks=FALSE,
	labels=abbreviate(rownames(x), 9),
	gap=0.5,
	line.main=3,
	cex.labels=NULL,
	font.labels=1,
	...) {
	
	# error checking:
	if (!is.logical(bounds))
		stop("bounds must be a logical.")
	if (!is.logical(boxBlocks))
		stop("boxBlocks must be a logical.")
	
	d <- dim(x)
	
	dots <- list(...)
	if ("main" %in% names(dots)) {
		main <- dots$main
	} else {
		main <- NULL
	}
	if ("oma" %in% names(dots)) {
		oma <- dots$oma
	} else {
		oma <- c(4, 4, ifelse(is.null(main), 4, 6), 4)
	}
	if (length(labels) != d[1])
		stop("labels does not match the dimensions of x.")
	
	opar <- par(mfrow=d,
		mar=rep.int(gap/2, 4),
		oma=oma)
	on.exit(par(opar))
	
	dev.hold()
	on.exit(dev.flush(),
		add=TRUE)
	
	for (i in seq_len(d[1])) {
		for (j in seq_len(d[2])) {
			c1 <- cumsum(c(0, x[i, i][[1]]))
			c2 <- cumsum(c(0, x[j, j][[1]]))
			
			plot(NA,
				xlim=c(1, c2[length(c2)]),
				ylim=c(1, c1[length(c1)]),
				axes=FALSE,
				xlab="",
				ylab="",
				xaxs="i",
				yaxs="i")
			box()
			
			if (j >= i && i==1 && (j %% 2)==0) {
				axis(3)
			}
			if (j >= i && j==d[2] && (i %% 2)==1) {
				axis(4)
			}
			if (i >= j && j==1 && (i %% 2)==0) {
				axis(2)
			}
			if (i >= j && i==d[2] && (j %% 2)==1) {
				axis(1)
			}
			
			if (i==j) {
				if (is.null(cex.labels)) {
					cex <- max(0.8,
						min(2,
							0.8/max(strwidth(labels, "figure"))))
				} else {
					cex <- cex.labels
				}
				
				text((c2[length(c2)] - 1)/2,
					(c1[length(c1)] - 1)/2,
					labels[i],
					cex=cex,
					font=font.labels)
			} else if (j > i) {
				if (bounds) {
					segments(c2 + 0.5, 1,
						c2 + 0.5, c1[length(c1)],
						lwd=0.2, col="gray")
					segments(1, c1 + 0.5,
						c2[length(c2)], c1 + 0.5,
						lwd=0.2, col="gray")
				}
				
				s <- x[i, j][[1]]
				if (dim(s)[1]==0)
					next
				s1 <- s[, "start1"] + c1[s[, "index1"]]
				s2 <- s[, "start2"] + c2[s[, "index2"]]
				width <- s[, "width"]
				e1 <- s1 + width - 1L
				strand <- s[, "strand"]
				e2 <- s2 + ifelse(strand,
					1L - width,
					width - 1L)
				
				segments(s2, s1, e2, e1,
					col=strand + 1)
				
				if (boxBlocks) {
					s <- x[j, i][[1]]
					rect(s[, "start2"] + c2[s[, "index2"]],
						s[, "start1"] + c1[s[, "index1"]],
						s[, "end2"] + c2[s[, "index2"]],
						s[, "end1"] + c1[s[, "index1"]],
						lwd=0.5)
				}
			} else { # i > j
				if (bounds) {
					segments(c2 + 0.5, 1,
						c2 + 0.5, c1[length(c1)],
						lwd=0.2, col="gray")
					segments(1, c1 + 0.5,
						c2[length(c2)], c1 + 0.5,
						lwd=0.2, col="gray")
				}
				
				s <- x[j, i][[1]]
				if (dim(s)[1]==0)
					next
				s1 <- s[, "start1"] + c2[s[, "index1"]]
				s2 <- s[, "start2"] + c1[s[, "index2"]]
				width <- s[, "width"]
				e1 <- s1 + width - 1L
				strand <- s[, "strand"]
				e2 <- s2 + ifelse(strand,
					1L - width,
					width - 1L)
				
				t <- x[i, j][[1]]
				s1[t[, "first_hit"]] <- t[, "start1"] + c2[t[, "index1"]]
				e1[t[, "last_hit"]] <- t[, "end1"] + c2[t[, "index1"]]
				s2[t[, "first_hit"]] <- ifelse(t[, "strand"]==0L,
					t[, "start2"],
					t[, "end2"]) + c1[t[, "index2"]]
				e2[t[, "last_hit"]] <- ifelse(t[, "strand"]==1L,
					t[, "start2"],
					t[, "end2"]) + c1[t[, "index2"]]
				
				o <- rank(x[i, j][[1]][, "score"],
					ties.method="first")
				l <- mapply(function(x, y) y - x + 1L,
					t[, "first_hit"],
					t[, "last_hit"])
				chains <- lapply(seq_along(l),
					function (x) {
						rep(o[x],
							l[x])
					})
				chains <- unlist(chains)
				if (length(chains) > 0) {
					m <- max(chains)
					chains <- m - chains + 1L
					cols <- rainbow(m,
						start=0.4,
						end=0.9,
						v=0.7)
					cols <- cols[chains]
				} else {
					cols <- character()
				}
				
				segments(s1, s2, e1, e2,
					col=cols)
				s1 <- c(s1, NA)
				s2 <- c(s2, NA)
				e1 <- c(NA, e1)
				e2 <- c(NA, e2)
				cols <- c(cols, NA)
				cols[cumsum(l) + 1L] <- NA
				segments(e1, e2, s1, s2,
					col=cols)
			}
		}
	}
	
	if (!is.null(main)) {
		font <- ifelse("font.main" %in% names(dots),
			dots$font.main,
			par("font.main"))
		cex <- ifelse("cex.main" %in% names(dots),
			dots$cex.main,
			par("cex.main"))
		mtext(main,
			3,
			line.main,
			outer=TRUE,
			at=0.5,
			cex=cex,
			font=font)
	}
	
	invisible(NULL)
}

plot.Synteny <- function(x,
	labels=abbreviate(rownames(x), 9),
	...) {
	
	d <- dim(x)
	if (length(labels) != d[1])
		stop("labels does not match the dimensions of x.")
	c <- lapply(diag(x), cumsum)
	
	dev.hold()
	on.exit(dev.flush(),
		add=TRUE)
	
	plot(NA,
		xlim=c(0.5, d[1] + 0.5),
		ylim=c(0, max(unlist(c)) + 1), # allow space
		xaxt='n',
		yaxs='i',
		xaxs='i',
		xlab="",
		ylab="Cumulative Position (nucleotides)",
		...)
	usr <- par()$usr
	text(x=seq_len(d[1]),
		y=usr[3] - 0.03*(usr[4] - usr[3]),
		labels=labels,
		srt=45,
		adj=c(1, 1),
		xpd=TRUE)
	
	for (i in seq_len(d[1])) {
		c1 <- c(0L, c[[i]])
		if (i < d[1]) {
			j <- i + 1L
			c2 <- c(0L, c[[j]])
			
			s <- x[j, i][[1]]
			if (dim(s)[1] > 0) {
				s1 <- s[, "start1"] + c1[s[, "index1"]]
				s2 <- s[, "start2"] + c2[s[, "index2"]]
				e1 <- s[, "end1"] + c1[s[, "index1"]]
				e2 <- s[, "end2"] + c2[s[, "index2"]]
				w <- which(s[, "strand"]==1)
				if (length(w) > 0) {
					temp <- e2[w]
					e2[w] <- s2[w]
					s2[w] <- temp
				}
				
				cols <- rainbow(length(s1),
					v=0.9,
					alpha=0.5)
				o <- order(s1)
				cols <- cols[o]
			} else {
				s1 <- integer()
			}
		} else if (d[1] > 2) {
			j <- 1L
			c2 <- c(0L, c[[j]])
			
			s <- x[i, j][[1]]
			if (dim(s)[1] > 0) {
				s1 <- s[, "start2"] + c1[s[, "index2"]]
				s2 <- s[, "start1"] + c2[s[, "index1"]]
				e1 <- s[, "end2"] + c1[s[, "index2"]]
				e2 <- s[, "end1"] + c2[s[, "index1"]]
				w <- which(s[, "strand"]==1)
				if (length(w) > 0) {
					temp <- e1[w]
					e1[w] <- s1[w]
					s1[w] <- temp
				}
				
				cols <- rainbow(length(s1),
					v=0.8,
					alpha=0.5)
				o <- order(s1)
				cols <- cols[o]
				
				segments(j - 0.1,
					s2,
					j - 0.9,
					s1,
					col=cols)
				segments(j - 0.1,
					e2,
					j - 0.9,
					e1,
					col=cols)
			} else {
				s1 <- integer()
			}
		}
		
		if (length(s1) > 0 &&
			(i < d[1] || d[1] > 2)) {
			segments(i + 0.1,
				s1,
				i + 0.9,
				s2,
				col=cols)
			segments(i + 0.1,
				e1,
				i + 0.9,
				e2,
				col=cols)
			
			rect(i,
				s1,
				i + 0.1,
				e1,
				col=cols,
				border=NA)
			rect(j,
				s2,
				j - 0.1,
				e2,
				col=cols,
				border=NA)
		}
		
		rect(i - 0.1,
			0,
			i + 0.1,
			c1[length(c1)])
		rect(i - 0.1,
			c1,
			i + 0.1,
			c1 + 1L,
			col="black")
	}
	
	invisible(NULL)
}
