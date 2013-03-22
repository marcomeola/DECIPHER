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
