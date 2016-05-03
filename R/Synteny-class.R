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
	colorBy=1,
	colorRamp=colorRampPalette(c("#FCF9EE", "#FFF272", "#FFAC28", "#EC5931", "#EC354D", "#ECA6B1")),
	barColor="#CCCCCC",
	barSides=ifelse(nrow(x) < 100, TRUE, FALSE),
	horizontal=TRUE,
	labels=abbreviate(rownames(x), 9),
	width=0.7,
	...) {
	
	d <- dim(x)
	if (length(labels) != d[1])
		stop("labels does not match the dimensions of x.")
	if (!is.logical(horizontal) && length(horizontal)==1)
		stop("horizontal must be a single logical.")
	if (!is.logical(barSides) && length(barSides)==1)
		stop("barSides must be a single logical.")
	if (!is.numeric(width) || length(width) != 1)
		stop("width must be a single number.")
	if (width <= 0 || width >= 1)
		stop("width must be between zero and one.")
	# set the width of each sequence bar
	half_width <- width/2 # between 0.0 and 0.5
	c <- lapply(diag(x), cumsum)
	if (is.character(colorBy)) {
		COLORBY <- c("neighbor", "frequency", "none")
		colorBy <- pmatch(colorBy[1], COLORBY)
		if (is.na(colorBy))
			stop("Invalid colorBy.")
		if (colorBy == -1)
			stop("Ambiguous colorBy.")
		colorBy <- -colorBy
	} else if (is.numeric(colorBy)) {
		if (length(colorBy) > 1)
			stop("colorBy must be a single number.")
		if (colorBy <= 0)
			stop("colorBy must be greater than zero.")
		if (floor(colorBy) != colorBy)
			stop("colorBy must be a whole number.")
		if (colorBy > d[1])
			stop("colorBy is out of range in x.")
	} else {
		stop("colorBy must be a character string or index.")
	}
	if (typeof(colorRamp) != "closure")
		stop("colorRamp must be a function.")
	
	dev.hold()
	on.exit(dev.flush(),
		add=TRUE)
	
	if (colorBy > 0) {
		dev_width <- dev.size("px")[1]
		c0 <- c(0L, c[[colorBy]])
		gradient <- colorRamp(dev_width + 1)
		endpoint <- seq(1, c0[length(c0)], length.out=dev_width + 1)
	} else if (colorBy == -2) {
		dev_width <- dev.size("px")[1]
		gradient <- colorRamp(dev_width + 1)
	}
	
	if (horizontal) {
		plot(NA,
			ylim=c(d[1] + 0.5, 0.5),
			xlim=c(0, max(unlist(c)) + 1), # allow space
			yaxt='n',
			yaxs='i',
			xaxs='i',
			ylab="",
			xlab="Cumulative Position (nucleotides)",
			...)
		usr <- par()$usr
		cex.labels <- min(1,
				0.7/(strheight(labels, "figure")[1]*d[1]))
		text(y=seq_len(d[1]),
			x=usr[2]*-0.015,
			labels=labels,
			adj=c(1, 0.5),
			cex=cex.labels,
			xpd=TRUE)
	} else {
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
		cex.labels <- min(1,
				0.7/(strheight(labels, "figure")[1]*d[1]))
		text(x=seq_len(d[1]),
			y=usr[3] - 0.01*(usr[4] - usr[3]),
			labels=labels,
			srt=45,
			adj=c(1, 1),
			cex=cex.labels,
			xpd=TRUE)
	}
	
	for (i in seq_len(d[1])) {
		c1 <- c(0L, c[[i]])
		
		if (colorBy == -1) { # color by neighbor
			if (i < d[1]) {
				if (i==1) {
					x0 <- i - half_width
					x1 <- i + half_width
					y0 <- 0
					y1 <- c1
					if (horizontal) {
						rect(y0, x0, y1, x1, col=barColor, border=NA)
					} else {
						rect(x0, y0, x1, y1, col=barColor, border=NA)
					}
				}
				
				x0 <- i + 1 - half_width
				x1 <- i + 1 + half_width
				y0 <- 0
				y1 <- tail(c[[i + 1]], n=1)
				if (horizontal) {
					rect(y0, x0, y1, x1, col=barColor, border=NA)
				} else {
					rect(x0, y0, x1, y1, col=barColor, border=NA)
				}
				
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
					
					cols <- colorRamp(length(s1))
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
					
					cols <- colorRamp(length(s1))
					o <- order(s1)
					cols <- cols[o]
					
					x0 <- j - half_width
					x1 <- j - (1 - half_width)
					y0 <- s2
					y1 <- s1
					if (horizontal) {
						segments(y0, x0, y1, x1, cols)
					} else {
						segments(x0, y0, x1, y1, cols)
					}
					
					y0 <- e2
					y1 <- e1
					if (horizontal) {
						segments(y0, x0, y1, x1, cols)
					} else {
						segments(x0, y0, x1, y1, cols)
					}
				} else {
					s1 <- integer()
				}
			}
			
			if (length(s1) > 0 &&
				(i < d[1] || d[1] > 2)) {
				x0 <- i + half_width
				x1 <- i + (1 - half_width)
				y0 <- s1
				y1 <- s2
				if (horizontal) {
					segments(y0, x0, y1, x1, cols)
				} else {
					segments(x0, y0, x1, y1, cols)
				}
				
				y0 <- e1
				y1 <- e2
				if (horizontal) {
					segments(y0, x0, y1, x1, cols)
				} else {
					segments(x0, y0, x1, y1, cols)
				}
				
				x0 <- i
				x1 <- i + half_width
				y0 <- s1
				y1 <- e1
				if (horizontal) {
					rect(y0, x0, y1, x1, col=cols, border=NA)
				} else {
					rect(x0, y0, x1, y1, col=cols, border=NA)
				}
				
				x0 <- j
				x1 <- j - half_width
				y0 <- s2
				y1 <- e2
				if (horizontal) {
					rect(y0, x0, y1, x1, col=cols, border=NA)
				} else {
					rect(x0, y0, x1, y1, col=cols, border=NA)
				}
			}
		} else if (colorBy > 0) { # color by position in reference
			endpoints <- seq(1, c1[length(c1)], length.out=dev_width + 1)
			if (i==colorBy) {
				cols <- gradient
			} else {
				cols <- rep(barColor, dev_width)
				
				if (i > colorBy) {
					s <- x[i, colorBy][[1]]
					if (dim(s)[1] > 0) {
						s1 <- s[, "start2"] + c1[s[, "index2"]]
						s0 <- s[, "start1"] + c0[s[, "index1"]]
						e1 <- s[, "end2"] + c1[s[, "index2"]]
						e0 <- s[, "end1"] + c0[s[, "index1"]]
					}
				} else { # i < colorBy
					s <- x[colorBy, i][[1]]
					if (dim(s)[1] > 0) {
						s1 <- s[, "start1"] + c1[s[, "index1"]]
						s0 <- s[, "start2"] + c0[s[, "index2"]]
						e1 <- s[, "end1"] + c1[s[, "index1"]]
						e0 <- s[, "end2"] + c0[s[, "index2"]]
					}
				}
				
				start1 <- sapply(s1, function (x) {
						which.min(abs(endpoints - x))
					})
				end1 <- sapply(e1, function (x) {
						which.min(abs(endpoints - x))
					})
				start0 <- sapply(s0, function (x) {
						which.min(abs(endpoint - x))
					})
				end0 <- sapply(e0, function (x) {
						which.min(abs(endpoint - x))
					})
				
				for (j in seq_len(dim(s)[1])) {
					if ((end0[j] - start0[j])==(end1[j] - start1[j])) {
						if (s[j, "strand"]==1) {
							cols[start1[j]:end1[j]] <- gradient[end0[j]:start0[j]]
						} else {
							cols[start1[j]:end1[j]] <- gradient[start0[j]:end0[j]]
						}
					} else { # interpolate colors
						if (s[j, "strand"]==1) {
							cols[start1[j]:end1[j]] <- colorRampPalette(gradient[end0[j]:start0[j]])(end1[j] - start1[j] + 1L)
						} else {
							cols[start1[j]:end1[j]] <- colorRampPalette(gradient[start0[j]:end0[j]])(end1[j] - start1[j] + 1L)
						}
					}
				}
			}
			
			x0 <- i - half_width
			x1 <- i + half_width
			y0 <- endpoints[-length(endpoints)]
			y1 <- endpoints[-1]
			if (horizontal) {
				rect(y0, x0, y1, x1, col=cols, border=NA)
			} else {
				rect(x0, y0, x1, y1, col=cols, border=NA)
			}
		} else if (colorBy == -2) { # color by frequency
			endpoints <- seq(1, c1[length(c1)], length.out=dev_width)
			freqs <- integer(dev_width)
			
			for (k in seq_len(d[1])) {
				if (k==i) {
					next
				} else {
					freq <- integer(dev_width)
					
					if (i > k) {
						s <- x[i, k][[1]]
						if (dim(s)[1] > 0) {
							s1 <- s[, "start2"] + c1[s[, "index2"]]
							e1 <- s[, "end2"] + c1[s[, "index2"]]
						}
					} else { # i < k
						s <- x[k, i][[1]]
						if (dim(s)[1] > 0) {
							s1 <- s[, "start1"] + c1[s[, "index1"]]
							e1 <- s[, "end1"] + c1[s[, "index1"]]
						}
					}
					
					start1 <- sapply(s1, function (x) {
							which.min(abs(endpoints - x))
						})
					end1 <- sapply(e1, function (x) {
							which.min(abs(endpoints - x))
						})
					
					for (j in seq_len(dim(s)[1]))
						freq[start1[j]:end1[j]] <- 1L
				}
				freqs <- freqs + freq
			}
			
			freqs <- freqs/(d[1] - 1)*(dev_width)
			cols <- rep(barColor, dev_width)
			w <- which(freqs > 0)
			if (length(w) > 0)
				cols[w] <- gradient[freqs[w]]
			
			x0 <- i - half_width
			x1 <- i + half_width
			y0 <- endpoints[-length(endpoints)]
			y1 <- endpoints[-1]
			if (horizontal) {
				rect(y0, x0, y1, x1, col=cols, border=NA)
			} else {
				rect(x0, y0, x1, y1, col=cols, border=NA)
			}
		}
		
		if (barSides) {
			if (horizontal) {
				x0 <- 1
				x1 <- c1[length(c1)]
				y0 <- i - half_width
				y1 <- i + half_width
			} else {
				x0 <- i - half_width
				x1 <- i + half_width
				y0 <- 1
				y1 <- c1[length(c1)]
			}
			if (i==1) {
				x01 <- x0
				x11 <- x1
				y01 <- y0
				y11 <- y1
				if (horizontal) {
					on.exit(segments(x01, c(y01, y11), x11, c(y01, y11)), add=TRUE)
				} else {
					on.exit(segments(c(x01, x11), y01, c(x01, x11), y11), add=TRUE)
				}
			} else {
				if (horizontal) {
					segments(x0, c(y0, y1), x1, c(y0, y1))
				} else {
					segments(c(x0, x1), y0, c(x0, x1), y1)
				}
			}
		}
		
		# delineate sequence ends
		if (horizontal) {
			x0 <- c1
			x1 <- c1 + 1L
			y0 <- i - half_width
			y1 <- i + half_width
		} else {
			x0 <- i - half_width
			x1 <- i + half_width
			y0 <- c1
			y1 <- c1 + 1L
		}
		if (i==1) {
			x02 <- x0
			x12 <- x1
			y02 <- y0
			y12 <- y1
			on.exit(rect(x02, y02, x12, y12, col="black"), add=TRUE)
		} else {
			rect(x0, y0, x1, y1, col="black")
		}
	}
	
	invisible(NULL)
}
