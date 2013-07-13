/****************************************************************************
 *                         Quick Matching Functions                         *
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

// DECIPHER header file
#include "DECIPHER.h"

// first matches of x[z...] == y[1]
SEXP multiMatch(SEXP x, SEXP y, SEXP z)
{	
	int i, size = length(x);
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *u = INTEGER(z);
	int start = -1, stop = -1;
	
	for (i = *u - 1; i < size; i++) {
		if (v[i] == w[0]) {
			start = i;
			stop = i;
			for (i = start + 1; i < size; i++) {
				if (v[i] == w[0]) {
					stop = i;
				} else {
					break;
				}
			}
			break;
		}
	}
	
	SEXP ans;
	if (start != -1) {
		PROTECT(ans = allocVector(INTSXP, stop - start + 1));
		int *rans = INTEGER(ans);
		for (i = start; i <= stop; i++) {
			rans[i - start] = i + 1;
		}
	} else {
		PROTECT(ans = allocVector(INTSXP, 0));
	}
	UNPROTECT(1);
	
	return ans;
}

// first matches of x[z...] >= y[1]
SEXP multiMatchUpper(SEXP x, SEXP y, SEXP z)
{	
	int i, size = length(x);
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *u = INTEGER(z);
	int start = -1, stop = -1;
	
	for (i = *u - 1; i < size; i++) {
		if (v[i] >= w[0]) {
			start = i;
			stop = i;
			for (i = start + 1; i < size; i++) {
				if (v[i] == v[start]) {
					stop = i;
				} else {
					break;
				}
			}
			break;
		}
	}
	
	SEXP ans;
	if (start != -1) {
		PROTECT(ans = allocVector(INTSXP, stop - start + 1));
		int *rans = INTEGER(ans);
		for (i = start; i <= stop; i++) {
			rans[i - start] = i + 1;
		}
	} else {
		PROTECT(ans = allocVector(INTSXP, 0));
	}
	UNPROTECT(1);
	
	return ans;
}

// first matches of x[...z] <= y[1]
SEXP multiMatchLower(SEXP x, SEXP y, SEXP z)
{	
	int i, size = length(x);
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *u = INTEGER(z);
	int start = -1, stop = -1;
	
	for (i = *u - 1; i >= 0; i--) {
		if (v[i] <= w[0]) {
			start = i;
			stop = i;
			for (i = start + 1; i < size; i++) {
				if (v[i] == v[start]) {
					stop = i;
				} else {
					break;
				}
			}
			break;
		}
	}
	
	SEXP ans;
	if (start != -1) {
		PROTECT(ans = allocVector(INTSXP, stop - start + 1));
		int *rans = INTEGER(ans);
		for (i = start; i <= stop; i++) {
			rans[i - start] = i + 1;
		}
	} else {
		PROTECT(ans = allocVector(INTSXP, 0));
	}
	UNPROTECT(1);
	
	return ans;
}

// index of first non-NA elements
SEXP multiMatchCharNotNA(SEXP x)
{	
	int i, size = length(x);
	int stop = 0;
	
	for (i = 0; i < size; i++) {
		if (STRING_ELT(x, i) != NA_STRING) {
			stop = i + 1;
		} else {
			break;
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, stop));
	int *rans = INTEGER(ans);
	for (i = 0; i < stop; i++) {
		rans[i] = i + 1;
	}
	
	UNPROTECT(1);
	
	return ans;
}

// same as x %in% y for integer vectors
SEXP intMatch(SEXP x, SEXP y)
{	
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int i, j;
	int size_x = length(x);
	int size_y = length(y);
	
	SEXP ans;
	PROTECT(ans = allocVector(LGLSXP, size_x));
	int *rans = INTEGER(ans);
	
	#pragma omp parallel for private(i, j) schedule(guided)
	for (i = 0; i < size_x; i++) {
		rans[i] = 0;
		for (j = 0; j < size_y; j++) {
			if (v[i] == w[j]) {
				rans[i] = 1;
				break;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// first match in y >= x[...]
SEXP firstMatchUpper(SEXP x, SEXP y)
{	
	int i, j, size_x = length(x), size_y = length(y);
	double *v = REAL(x);
	double *w = REAL(y);
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, size_x));
	int *rans = INTEGER(ans);
	
	#pragma omp parallel for private(i, j) schedule(guided)
	for (i = 0; i < size_x; i++) {
		rans[i] = NA_INTEGER;
		for (j = 0; j < size_y; j++) {
			if (w[j] >= v[i]) {
				rans[i] = j + 1;
				break;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// matrix of d[i, j] = 1 - length x[i] %in% x[j] / min(length)
// requires a list of ordered integers
SEXP matchLists(SEXP x, SEXP verbose, SEXP pBar)
{	
	int i, j, size_x = length(x), before, v, *rPercentComplete;
	int o, p, start, lx, ly, *X, *Y, count;
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size_x, size_x));
	double *rans = REAL(ans);
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	
	if (v) { // initialize progress variables
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	for (i = 0; i < size_x; i++)
		*(rans + i*size_x + i) = 0;
	
	for (i = 0; i < size_x; i++) {
		#pragma omp parallel for private(j, o, p, start, count, X, Y, lx, ly) schedule(guided)
		for (j = i + 1; j < size_x; j++) {
			X = INTEGER(VECTOR_ELT(x, i));
			Y = INTEGER(VECTOR_ELT(x, j));
			lx = length(VECTOR_ELT(x, i));
			ly = length(VECTOR_ELT(x, j));
			
			if (lx > 0 && ly > 0) {
				int first = -1;
				int last = -1;
				for (o = 0; o < lx; o++) {
					if (X[o] >= Y[0]) {
						first = o;
						break;
					}
				}
				if (first == -1) { // no overlap
					*(rans + i*size_x + j) = NA_REAL;
					*(rans + i + j*size_x) = NA_REAL;
					continue;
				}
				
				for (o = lx - 1; o >= 0; o--) {
					if (X[o] <= Y[ly - 1]) {
						last = o;
						break;
					}
				}
				if (last == -1) { // no overlap
					*(rans + i*size_x + j) = NA_REAL;
					*(rans + i + j*size_x) = NA_REAL;
					continue;
				}
				
				int lz = last - first + 1;
				
				count = 0;
				start = 0;
				for (o = first; o <= last; o++) {
					for (p = start; p < ly; p++) {
						if (X[o] == Y[p]) {
							count++;
							start = p + 1;
							break;
						} else if (Y[p] > X[o]) {
							break;
						}
					}
				}
				
				*(rans + i*size_x + j) = 1 - (double)count/(double)lz;
			} else {
				*(rans + i*size_x + j) = NA_REAL;
			}
			*(rans + i + j*size_x) = *(rans + i*size_x + j);
		}
		
		if (v) {
			// print the percent completed so far
			*rPercentComplete = floor(100*(double)((i + 1)*size_x+(i + 1))/((size_x - 1)*size_x+(size_x - 1)));
			
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}
