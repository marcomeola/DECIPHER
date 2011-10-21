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

static double distance(const cachedCharSeq *P, const cachedCharSeq *S, int start, int end, int pGapsGaps, int pGapLetters)
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
	for (i = start, j = start, p = (P->seq + start), s = (S->seq + start);
	     (i < (P->length - end)) && (i < (S->length - end));
	     i++, j++, p++, s++)
	{
		count++; // increment the length covered
		if (!((*p) & (*s))) { // sequences are not equal
			if ((*p) & 0x10 || (*s) & 0x10) { // gap-letter match
				if (!pGapLetters) { // don't penalize gap-letter matches
					gapLetterMatches++; // don't include gap-letter matches in length
				} else { // penalize gap-letter matches
					mismatches++; // count gap-letter matches as mis-matches
				}
			} else {
				mismatches++; // mis-match
			}
		} else {
			if ((*p) & 0x10 && (*s) & 0x10) { // gap-gap match
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

static int frontTerminalGaps(const cachedCharSeq *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->seq;
	     i < P->length;
	     i++, p++)
	{
		if ((*p) & 0x10) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGaps(const cachedCharSeq *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->seq + P->length - 1);
	     i >= 0;
	     i--, p--)
	{
		if ((*p) & 0x10) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

//ans_start <- .Call("distMatrix", myDNAStringSet, pBar, PACKAGE="DECIPHER")
SEXP distMatrix(SEXP x, SEXP terminalGaps, SEXP penalizeGapGaps, SEXP penalizeGapLetters, SEXP verbose, SEXP pBar)
{
	cachedXStringSet x_set;
	cachedCharSeq x_i, x_j;
	int x_length, start, end, i, j, seqLength_i, seqLength_j;
	int pGapLetters, pGapsGaps, tGaps;
	int soFar, before, v, *rPercentComplete;
	double *rans;
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
	
	x_set = cache_XStringSet(x);
	x_length = get_cachedXStringSet_length(&x_set);
	int gapLengths[x_length][2];
	
	if (x_length < 2) { // there is only one sequence
		PROTECT(ans = NEW_INTEGER(0));
	} else {
		PROTECT(ans = allocMatrix(REALSXP, x_length, x_length));
		rans = REAL(ans);
		
		// find the terminal gap lengths
		// always needed to identify no-overlap
		for (i = 0; i < x_length; i++) {
			x_i = get_cachedXStringSet_elt(&x_set, i);
			gapLengths[i][1] = frontTerminalGaps(&x_i);
			gapLengths[i][2] = endTerminalGaps(&x_i);
			//Rprintf("\nstart:%dend:%d",gapLengths[i][1],gapLengths[i][2]);
		}
		
		tGaps = asLogical(terminalGaps);
		pGapsGaps = asLogical(penalizeGapGaps);
		pGapLetters = asLogical(penalizeGapLetters);
		for (i = 0; i < (x_length - 1); i++) {
			// extract each ith DNAString from the DNAStringSet
			x_i = get_cachedXStringSet_elt(&x_set, i);
			seqLength_i = x_i.length;
			
			#pragma omp parallel for private(j,x_j,seqLength_j,start,end) schedule(guided)
			for (j = (i+1); j < x_length; j++) {
				// extract each jth DNAString from the DNAStringSet
				x_j = get_cachedXStringSet_elt(&x_set, j);
				seqLength_j = x_j.length;
				
				// find the distance for each row of the matrix
				if ((seqLength_i - gapLengths[i][2]) < gapLengths[j][1] ||
					gapLengths[i][1] > (seqLength_j - gapLengths[j][2])) {
					// no overlap between sequences
					rans[i + x_length*j] = NA_REAL;
				} else {
					if (!tGaps) { // don't include terminal gaps
						// find the intersection of both string's ranges
						// to shorten the sequence comparison for speed
						if (gapLengths[i][1] >= gapLengths[j][1]) {
							start = gapLengths[i][1];
						} else {
							start = gapLengths[j][1];
						}
						if ((seqLength_i - gapLengths[i][2]) <= (seqLength_j - gapLengths[j][2])) {
							end = gapLengths[i][2];
						} else {
							end = gapLengths[j][2];
						}
						rans[i + x_length*j] = distance(&x_i, &x_j, start, end, pGapsGaps, pGapLetters);
					} else { // use whole sequence including terminal gaps
						rans[i + x_length*j] = distance(&x_i, &x_j, 0, 0, pGapsGaps, pGapLetters);
					}

				}
				// make the matrix symetrical
				rans[j + x_length*i] = rans[i + x_length*j];
			}
			// set the matrix diagonal to zero distance
			rans[i*x_length+i] = 0;
			
			if (v) { // print the percent completed so far
				soFar = (i + 1)*x_length+(i + 1);
				*rPercentComplete = floor(100*(double)soFar/((x_length - 1)*x_length+(x_length - 1)));
				if (*rPercentComplete > before) { // when the percent has changed
					// tell the progress bar to update in the R console
					eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
					before = *rPercentComplete;
				}
			} else {
				void R_CheckUserInterrupt(void);
			}
		}
		// set the last element of the diagonal to zero
		rans[(x_length - 1)*x_length+(x_length - 1)] = 0;
	}
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;	
}

//ans_start <- .Call("gaps", myDNAStringSet, PACKAGE="DECIPHER")
SEXP gaps(SEXP x)
{
	cachedXStringSet x_set;
	cachedCharSeq x_i;
	int x_length, i;
	double *rans;
	SEXP ans;
	
	x_set = cache_XStringSet(x);
	x_length = get_cachedXStringSet_length(&x_set);

	PROTECT(ans = allocMatrix(REALSXP, x_length, 3));
	rans = REAL(ans);
	
	// find the three lengths for each sequence
	for (i = 0; i < x_length; i++) {
		x_i = get_cachedXStringSet_elt(&x_set, i);
		rans[i + x_length*0] = frontTerminalGaps(&x_i);
		rans[i + x_length*1] = endTerminalGaps(&x_i);
		rans[i + x_length*2] = x_i.length - rans[i + x_length*1] - rans[i + x_length*0];
	}
	
	UNPROTECT(1);
	return ans;
}
