/****************************************************************************
 *                         Creates Distance Marix                           *
 *                           Author: Erik Wright                            *
 ****************************************************************************/

/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

// for math functions
#include <math.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

static double distance(const Chars_holder *P, const Chars_holder *S, int start, int end, int pGapsGaps, int pGapLetters)
{
	double distance;
	int i, j, mismatches, gapGapMatches, gapLetterMatches, count;
	const char *p, *s;
	
	distance = 0;
	gapGapMatches = 0;
	gapLetterMatches = 0;
	mismatches = 0;
	count = 0;
	
	// walk along the sequence from (position start + 1) to (length - end - 1)
	for (i = start, j = start, p = (P->ptr + start), s = (S->ptr + start);
	     (i < (P->length - end)) && (i < (S->length - end));
	     i++, j++, p++, s++)
	{
		count++; // increment the length covered
		if (!((*p) & (*s))) { // sequences are not equal
			if (((*p) & 0x40 && (*s) & 0x10) || ((*p) & 0x10 && (*s) & 0x40)) { // gap-gap match
				if (!pGapsGaps) { // don't penalize gap-gap matches
					gapGapMatches++; // don't include gap-gap matches in length
				} else { // penalize gap-gap matches
					mismatches++; // count gap-gap matches as mis-matches
				}
			} else if ((*p) & 0x10 || (*s) & 0x10 || (*p) & 0x40 || (*s) & 0x40) { // gap-letter match
				if (!pGapLetters) { // don't penalize gap-letter matches
					gapLetterMatches++; // don't include gap-letter matches in length
				} else { // penalize gap-letter matches
					mismatches++; // count gap-letter matches as mis-matches
				}
			} else {
				mismatches++; // mis-match
			}
		} else { // sequences are equal
			if (((*p) & 0x10 && (*s) & 0x10) || ((*p) & 0x40 && (*s) & 0x40)) { // gap-gap match
				if (!pGapsGaps) { // don't penalize gap-gap matches
					gapGapMatches++; // don't include gap-gap matches in length
				} else { // penalize gap-gap matches
					mismatches++; // count gap-gap matches as mis-matches
				}
			}
		}
	}
	
	//Rprintf("start%d end%d",start,end);
	//Rprintf("\nmismatches:%d gapGapMatches:%d gapLetterMatches:%d count:%d",mismatches,gapGapMatches,gapLetterMatches,count);
	
	// calculate distance as the percent mis-matches
	distance = (double)mismatches/((double)count - (double)gapGapMatches - (double)gapLetterMatches);
	return distance;
}

static double distanceAA(const Chars_holder *P, const Chars_holder *S, int start, int end, int pGapsGaps, int pGapLetters)
{
	double distance;
	int i, j, mismatches, gapGapMatches, gapLetterMatches, count;
	const char *p, *s;
	
	distance = 0;
	gapGapMatches = 0;
	gapLetterMatches = 0;
	mismatches = 0;
	count = 0;
	
	// walk along the sequence from (position start + 1) to (length - end - 1)
	for (i = start, j = start, p = (P->ptr + start), s = (S->ptr + start);
	     (i < (P->length - end)) && (i < (S->length - end));
	     i++, j++, p++, s++)
	{
		count++; // increment the length covered
		if ((*p) ^ (*s) && // sequences are not equal
			!(!((*p) ^ 0x58) && !(!((*s) ^ 0x2D) || !((*s) ^ 0x2B) || !((*s) ^ 0x2A))) && !(!((*s) ^ 0x58) && !(!((*p) ^ 0x2D) || !((*p) ^ 0x2B) || !((*p) ^ 0x2A))) && // not (X && !(non-letter))
			!(!((*p) ^ 0x42) && (!((*s) ^ 0x4E) || !((*s) ^ 0x44))) && !(!((*s) ^ 0x42) && (!((*p) ^ 0x4E) || !((*p) ^ 0x44))) && // not (B && (N or D))
			!(!((*p) ^ 0x4A) && (!((*s) ^ 0x49) || !((*s) ^ 0x4C))) && !(!((*s) ^ 0x4A) && (!((*p) ^ 0x49) || !((*p) ^ 0x4C))) && // not (J && (I or L))
			!(!((*p) ^ 0x5A) && (!((*s) ^ 0x51) || !((*s) ^ 0x45))) && !(!((*s) ^ 0x5A) && (!((*p) ^ 0x51) || !((*p) ^ 0x45)))) { // not (Z && (Q or E))
			if ((!((*p) ^ 0x2D) && !((*s) ^ 0x2E)) || (!((*p) ^ 0x2E) && !((*s) ^ 0x2D))) { // gap-gap match
				if (!pGapsGaps) { // don't penalize gap-gap matches
					gapGapMatches++; // don't include gap-gap matches in length
				} else { // penalize gap-gap matches
					mismatches++; // count gap-gap matches as mis-matches
				}
			} else if (!((*p) ^ 0x2D) || !((*s) ^ 0x2D) || !((*p) ^ 0x2E) || !((*s) ^ 0x2E)) { // gap-letter match
				if (!pGapLetters) { // don't penalize gap-letter matches
					gapLetterMatches++; // don't include gap-letter matches in length
				} else { // penalize gap-letter matches
					mismatches++; // count gap-letter matches as mis-matches
				}
			} else {
				mismatches++; // mis-match
			}
		} else { // sequences are equal
			if ((!((*p) ^ 0x2D) && !((*s) ^ 0x2D)) || (!((*p) ^ 0x2E) && !((*s) ^ 0x2E))) { // gap-gap match
				if (!pGapsGaps) { // don't penalize gap-gap matches
					gapGapMatches++; // don't include gap-gap matches in length
				} else { // penalize gap-gap matches
					mismatches++; // count gap-gap matches as mis-matches
				}
			}
		}
	}
	
	//Rprintf("start%d end%d",start,end);
	//Rprintf("\nmismatches:%d gapGapMatches:%d gapLetterMatches:%d count:%d",mismatches,gapGapMatches,gapLetterMatches,count);
	
	// calculate distance as the percent mis-matches
	distance = (double)mismatches/((double)count - (double)gapGapMatches - (double)gapLetterMatches);
	return distance;
}

static int frontTerminalGaps(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
	     i < P->length;
	     i++, p++)
	{
		if ((*p) & 0x10 || (*p) & 0x40) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGaps(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
	     i >= 0;
	     i--, p--)
	{
		if ((*p) & 0x10 || (*p) & 0x40) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int frontTerminalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
	     i < P->length;
	     i++, p++)
	{
		if (!((*p) ^ 0x2D) || !((*p) ^ 0x2E)) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
	     i >= 0;
	     i--, p--)
	{
		if (!((*p) ^ 0x2D) || !((*p) ^ 0x2E)) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

SEXP distMatrix(SEXP x, SEXP t, SEXP terminalGaps, SEXP penalizeGapGaps, SEXP penalizeGapLetters, SEXP fullMatrix, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	XStringSet_holder x_set;
	Chars_holder x_i, x_j;
	int x_length, start, end, i, j, seqLength_i, seqLength_j, last;
	int pGapLetters, pGapsGaps, tGaps, fM = asLogical(fullMatrix);
	int before, v, *rPercentComplete;
	double *rans, soFar;
	int nthreads = asInteger(nThreads);
	SEXP ans, percentComplete, utilsPackage;
	v = asLogical(verbose);
	if (v) { // percent complete variables
		soFar = 0;
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	int gapLengths[x_length][2];
	
	if (x_length < 2) { // there is only one sequence
		PROTECT(ans = NEW_INTEGER(0));
	} else {
		if (fM) {
			last = x_length - 1;
			PROTECT(ans = allocMatrix(REALSXP, x_length, x_length));
		} else {
			last = 1;
			PROTECT(ans = allocMatrix(REALSXP, 1, x_length));
		}
		rans = REAL(ans);
		
		// find the terminal gap lengths
		// always needed to identify no-overlap
		if (asInteger(t)==3) { // AAStringSet
			for (i = 0; i < x_length; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				gapLengths[i][0] = frontTerminalGapsAA(&x_i);
				gapLengths[i][1] = endTerminalGapsAA(&x_i);
				//Rprintf("\nstart:%dend:%d",gapLengths[i][0],gapLengths[i][1]);
			}
		} else { // DNAStringSet or RNAStringSet
			for (i = 0; i < x_length; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				gapLengths[i][0] = frontTerminalGaps(&x_i);
				gapLengths[i][1] = endTerminalGaps(&x_i);
				//Rprintf("\nstart:%dend:%d",gapLengths[i][0],gapLengths[i][1]);
			}
		}
		
		tGaps = asLogical(terminalGaps);
		pGapsGaps = asLogical(penalizeGapGaps);
		pGapLetters = asLogical(penalizeGapLetters);
		for (i = 0; i < last; i++) {
			// extract each ith DNAString from the DNAStringSet
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			seqLength_i = x_i.length;
			
			#pragma omp parallel for private(j,x_j,seqLength_j,start,end) schedule(guided) num_threads(nthreads)
			for (j = (i+1); j < x_length; j++) {
				// extract each jth DNAString from the DNAStringSet
				x_j = get_elt_from_XStringSet_holder(&x_set, j);
				seqLength_j = x_j.length;
				
				// find the distance for each row of the matrix
				if ((seqLength_i - gapLengths[i][1]) <= gapLengths[j][0] ||
					gapLengths[i][0] >= (seqLength_j - gapLengths[j][1])) {
					// no overlap between sequences
					if (tGaps) { // include terminal gaps
						rans[j + x_length*i] = 1;
					} else {
						rans[j + x_length*i] = NA_REAL;
					}
				} else {
					if (!tGaps) { // don't include terminal gaps
						// find the intersection of both string's ranges
						// to shorten the sequence comparison for speed
						if (gapLengths[i][0] >= gapLengths[j][0]) {
							start = gapLengths[i][0];
						} else {
							start = gapLengths[j][0];
						}
						if ((seqLength_i - gapLengths[i][1]) <= (seqLength_j - gapLengths[j][1])) {
							end = gapLengths[i][1];
						} else {
							end = gapLengths[j][1];
						}
						
					} else { // use whole sequence including terminal gaps
						start = 0;
						end = 0;
					}
					if (asInteger(t)==3) { // AAStringSet
						rans[j + x_length*i] = distanceAA(&x_i, &x_j, start, end, pGapsGaps, pGapLetters);
					} else {
						rans[j + x_length*i] = distance(&x_i, &x_j, start, end, pGapsGaps, pGapLetters);
					}
				}
			}
			if (fM) // make the matrix symetrical
				for (j = (i+1); j < x_length; j++)
					rans[i + x_length*j] = rans[j + x_length*i];
			
			// set the matrix diagonal to zero distance
			rans[i*x_length+i] = 0;
			
			if (v) { // print the percent completed so far
				soFar = (2*last - i)*(i + 1);
				*rPercentComplete = floor(100*soFar/(last*(last + 1)));
				if (*rPercentComplete > before) { // when the percent has changed
					// tell the progress bar to update in the R console
					eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
					before = *rPercentComplete;
				}
			} else {
				R_CheckUserInterrupt();
			}
		}
		if (fM) // set the last element of the diagonal to zero
			rans[(x_length - 1)*x_length+(x_length - 1)] = 0;
	}
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;	
}

//ans_start <- .Call("gaps", myDNAStringSet, 1L, PACKAGE="DECIPHER")
SEXP gaps(SEXP x, SEXP t)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i;
	double *rans;
	SEXP ans;
	
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);

	PROTECT(ans = allocMatrix(REALSXP, x_length, 3));
	rans = REAL(ans);
	
	// find the three lengths for each sequence
	if (asInteger(t)==3) { // AAStringSet
		for (i = 0; i < x_length; i++) {
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			rans[i + x_length*0] = frontTerminalGapsAA(&x_i);
			rans[i + x_length*1] = endTerminalGapsAA(&x_i);
			rans[i + x_length*2] = x_i.length - rans[i + x_length*1] - rans[i + x_length*0];
		}
	} else { // DNAStringSet or RNAStringSet
		for (i = 0; i < x_length; i++) {
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			rans[i + x_length*0] = frontTerminalGaps(&x_i);
			rans[i + x_length*1] = endTerminalGaps(&x_i);
			rans[i + x_length*2] = x_i.length - rans[i + x_length*1] - rans[i + x_length*0];
		}
	}
	
	UNPROTECT(1);
	return ans;
}

//ans_start <- .Call("firstSeqsEqual", dna1, dna2, start1, end1, start2, end2, first1, first2, PACKAGE="DECIPHER")
SEXP firstSeqsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP first_x, SEXP first_y)
{	
	int i, j;
	XStringSet_holder x_set;
	XStringSet_holder y_set;
	Chars_holder x_i, y_i;
	int sx = asInteger(start_x);
	int ex = asInteger(end_x);
	int fx = asInteger(first_x);
	int sy = asInteger(start_y);
	int ey = asInteger(end_y);
	int fy = asInteger(first_y);
	
	
	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1));
	int *rans;
	rans = INTEGER(ans);
	*(rans) = 1; // equal
	if ((sx - ex) != (sy - ey)) { // different sequence lengths
		*(rans) = 0; // not equal
	} else {
		x_set = hold_XStringSet(x);
		y_set = hold_XStringSet(y);
		x_i = get_elt_from_XStringSet_holder(&x_set, fx - 1);
		y_i = get_elt_from_XStringSet_holder(&y_set, fy - 1);
		for (i = sx - 1, j = sy - 1;
			 i < ex; // i <= ex - 1 covers j <= ey - 1 because equal length
			 i++, j++) {
			if (x_i.ptr[i] != y_i.ptr[j]) {
				*(rans) = 0; // not equal
				break;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP firstSeqsGapsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP t, SEXP first_x, SEXP first_y)
{	
	int i, j;
	XStringSet_holder x_set;
	XStringSet_holder y_set;
	Chars_holder x_i, y_i;
	int sx = asInteger(start_x);
	int ex = asInteger(end_x);
	int fx = asInteger(first_x);
	int sy = asInteger(start_y);
	int ey = asInteger(end_y);
	int fy = asInteger(first_y);
	
	
	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1));
	int *rans;
	rans = INTEGER(ans);
	*(rans) = 1; // equal
	if ((sx - ex) != (sy - ey)) { // different sequence lengths
		*(rans) = 0; // not equal
	} else {
		x_set = hold_XStringSet(x);
		y_set = hold_XStringSet(y);
		x_i = get_elt_from_XStringSet_holder(&x_set, fx - 1);
		y_i = get_elt_from_XStringSet_holder(&y_set, fy - 1);
		if (asInteger(t)==3) { // AAStringSet
			for (i = sx - 1, j = sy - 1;
				 i < ex; // i <= ex - 1 covers j <= ey - 1 because equal length
				 i++, j++) {
				if ((!((*((char *)x_i.ptr + i)) ^ 0x2D) || !((*((char *)x_i.ptr + i)) ^ 0x2E)) ^
					(!((*((char *)y_i.ptr + j)) ^ 0x2D) || !((*((char *)y_i.ptr + j)) ^ 0x2E))) {
						*(rans) = 0; // not equal
						break;
					}
			}
		} else { // DNAStringSet or RNAStringSet
			for (i = sx - 1, j = sy - 1;
				 i < ex; // i <= ex - 1 covers j <= ey - 1 because equal length
				 i++, j++) {
				if (((*((char *)x_i.ptr + i)) & 0x10 || (*((char *)x_i.ptr + i)) & 0x40) ^
					((*((char *)y_i.ptr + j)) & 0x10 || (*((char *)y_i.ptr + j)) & 0x40)) {
					*(rans) = 0; // not equal
					break;
				}
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP firstSeqsPosEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP t, SEXP first_x, SEXP first_y)
{
	XStringSet_holder x_set;
	XStringSet_holder y_set;
	Chars_holder x_i, y_i;
	int sx = asInteger(start_x);
	int ex = asInteger(end_x);
	int fx = asInteger(first_x);
	int sy = asInteger(start_y);
	int ey = asInteger(end_y);
	int fy = asInteger(first_y);
	int type = asInteger(t);
	
	int sizex = 100, sizey = 100;
	int *nx = Calloc(sizex, int); // number of gaps in x
	int *px = Calloc(sizex, int); // position of gaps in x
	int *ny = Calloc(sizey, int); // number of gaps in y
	int *py = Calloc(sizey, int); // position of gaps in y
	
	int i = sx - 1, j = sy - 1; // position in x or y
	int cx = 0, cy = 0; // number of sites
	int ci = 1, cj = 1; // include site in count
	int gx = -1, gy = -1; // index of current gap
	int itx = 0, ity = 0; // iterate position
	
	x_set = hold_XStringSet(x);
	y_set = hold_XStringSet(y);
	x_i = get_elt_from_XStringSet_holder(&x_set, fx - 1);
	y_i = get_elt_from_XStringSet_holder(&y_set, fy - 1);
	
	while (i < ex && j < ey) {
		if (ci) {
			if (type==3) { // AAStringSet
				if (!(!((*((char *)x_i.ptr + i)) ^ 0x2D) || !((*((char *)x_i.ptr + i)) ^ 0x2E))) {
					cx++; // site
				}
			} else { // DNAStringSet or RNAStringSet
				if (!((*((char *)x_i.ptr + i)) & 0x10 || (*((char *)x_i.ptr + i)) & 0x40)) {
					cx++; // site
				}
			}
		}
		if (cj) {
			if (type==3) { // AAStringSet
				if (!(!((*((char *)y_i.ptr + j)) ^ 0x2D) || !((*((char *)y_i.ptr + j)) ^ 0x2E))) {
					cy++; // site
				}
			} else { // DNAStringSet or RNAStringSet
				if (!((*((char *)y_i.ptr + j)) & 0x10 || (*((char *)y_i.ptr + j)) & 0x40)) {
					cy++; // site
				}
			}
		}
		
		if (cx > cy) {
			j++;
			if (itx) {
				nx[gx]++;
			} else {
				ci = 0;
				itx = 1;
				gx++;
				
				if (gx >= (sizex - 1)) {
					sizex += 100;
					nx = Realloc(nx, sizex, int);
					px = Realloc(px, sizex, int);
				}
				
				nx[gx] = 1;
				px[gx] = i + 1;
			}
		} else if (cx < cy) {
			i++;
			if (ity) {
				ny[gy]++;
			} else {
				cj = 0;
				ity = 1;
				gy++;
				
				if (gy >= (sizey - 1)) {
					sizey += 100;
					ny = Realloc(ny, sizey, int);
					py = Realloc(py, sizey, int);
				}
				
				ny[gy] = 1;
				py[gy] = j + 1;
			}
		} else {
			ci = 1;
			cj = 1;
			i++;
			j++;
			itx = 0;
			ity = 0;
		}
	}
	
	if (i < ex) {
		gy++;
		ny[gy] = ex - i;
		py[gy] = ey + 1;
		
		if (!ci)
			i++;
			
			while (i < ex) {
				if (type==3) { // AAStringSet
					if (!(!((*((char *)x_i.ptr + i)) ^ 0x2D) || !((*((char *)x_i.ptr + i)) ^ 0x2E))) {
						cx++; // site
					}
				} else { // DNAStringSet or RNAStringSet
					if (!((*((char *)x_i.ptr + i)) & 0x10 || (*((char *)x_i.ptr + i)) & 0x40)) {
						cx++; // site
					}
				}
				i++;
			}
	}
	
	if (j < ey) {
		gx++;
		nx[gx] = ey - j;
		px[gx] = ex + 1;
		
		if (!cj)
			j++;
			
			while (j < ey) {
				if (type==3) { // AAStringSet
					if (!(!((*((char *)y_i.ptr + j)) ^ 0x2D) || !((*((char *)y_i.ptr + j)) ^ 0x2E))) {
						cy++; // site
					}
				} else { // DNAStringSet or RNAStringSet
					if (!((*((char *)y_i.ptr + j)) & 0x10 || (*((char *)y_i.ptr + j)) & 0x40)) {
						cy++; // site
					}
				}
				j++;
			}
	}
	
	SEXP ret_list, ans;
	if (cx==cy) { // same number of sites
		PROTECT(ret_list = allocVector(VECSXP, 4));
		int *rans;
		
		PROTECT(ans = allocVector(INTSXP, gx + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gx; i++) {
			*(rans + i) = px[i];
		}
		SET_VECTOR_ELT(ret_list, 0, ans);
		PROTECT(ans = allocVector(INTSXP, gx + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gx; i++) {
			*(rans + i) = nx[i];
		}
		SET_VECTOR_ELT(ret_list, 1, ans);
		PROTECT(ans = allocVector(INTSXP, gy + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gy; i++) {
			*(rans + i) = py[i];
		}
		SET_VECTOR_ELT(ret_list, 2, ans);
		PROTECT(ans = allocVector(INTSXP, gy + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gy; i++) {
			*(rans + i) = ny[i];
		}
		SET_VECTOR_ELT(ret_list, 3, ans);
	}
	
	Free(nx);
	Free(px);
	Free(ny);
	Free(py);
	
	if (cx==cy) { // same number of sites
		UNPROTECT(5);
		
		return ret_list;
	} else {
		return R_NilValue;
	}
}
