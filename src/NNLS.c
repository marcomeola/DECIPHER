/****************************************************************************
 *                         Solves SCA NNLS PROBLEM                          *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

// for math functions
#include <math.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("NNLS", i, j, x, nrows, ncols, b, tol, verbose, pBar, PACKAGE="DECIPHER")
SEXP NNLS(SEXP row, SEXP col, SEXP value, SEXP nrows, SEXP ncols, SEXP b, SEXP tol, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	int i, j, k, before, v, *rPercentComplete;
	int *rows = INTEGER(row);
	int *cols = INTEGER(col);
	double *values = REAL(value);
	unsigned long int n = asInteger(ncols);
	unsigned long int r = asInteger(nrows);
	int reps = length(b)/r;
	R_len_t l = length(row);
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	int nthreads = asInteger(nThreads);
	
	float *H = Calloc(n*n, float); // initialized to zero
	
	if (v) { // initialize progress variables
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	// calculate H = t(A) %*% A
	
	int start = 0, stop = 0;
	for (i = 0; i < l; i++) {
		if (i == stop) {
			start = i;
			stop = l;
			for (j = start + 1; j < l; j++) {
				if (*(rows + i) != *(rows + j)) {
					stop = j;
					break;
				}
			}
			
			if (v) {
				// print the percent completed so far
				*rPercentComplete = floor(100*((double)i + 1)/(double)l);
				
				if (*rPercentComplete > before) { // when the percent has changed
					// tell the progress bar to update in the R console
					eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
					before = *rPercentComplete;
				}
			} else {
				R_CheckUserInterrupt();
			}
		}
		
		#pragma omp parallel for private(j) schedule(guided) num_threads(nthreads)
		for (j = start; j < stop; j++)
			H[(*(cols + i) - 1)*n + *(cols + j) - 1] += *(values + i) * *(values + j);
	}
	
	if (v) { // ensure 100% completion
		// print the percent completed so far
		*rPercentComplete = 100;
		
		if (*rPercentComplete > before) { // when the percent has changed
			// tell the progress bar to update in the R console
			eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
			before = *rPercentComplete;
		}
	} else {
		R_CheckUserInterrupt();
	}
	
	double *B = REAL(b);
	SEXP ans;
	//PROTECT(ans = allocVector(REALSXP, n));
	PROTECT(ans = allocMatrix(REALSXP, n, reps));
	double *rans = REAL(ans);
	
	for (j = 0; j < reps; j++) {
		int startIndex1 = j*n;
		int startIndex2 = j*r;
		double *mu = (double *) R_alloc(n, sizeof(double));
		
		// calculate mu = -t(A) %*% b
		for (i = 0; i < n; i++)
			*(mu + i) = 0;
		for (i = 0; i < l; i++)
			*(mu + (*(cols + i) - 1)) -= *(values + i) * *(B + (*(rows + i) - 1) + startIndex2);
		
		// Solve NNLS
		
		double temp;
		double ep = asReal(tol);
		double error = ep;
		
		for (i = 0; i < n; i++)
			*(rans + i + startIndex1) = 0;
		
		while (error >= ep) {
			error = 0;
			for (i = 0; i < n; i++) {
				temp = *(rans + i + startIndex1);
				*(rans + i + startIndex1) -= *(mu + i) / H[i*n + i];
				if (*(rans + i + startIndex1) < 0)
					*(rans + i + startIndex1) = 0;
				if (*(rans + i + startIndex1) != temp) {
					if (fabs(*(rans + i + startIndex1) - temp) > error)
						error = fabs(*(rans + i + startIndex1) - temp);
					for (k = 0; k < n; k++)
						*(mu + k) += (*(rans + i + startIndex1) - temp) * H[i*n + k];
				}
			}
			R_CheckUserInterrupt();
		}
	}
	
	Free(H);
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}

//ans_start <- .Call("sparseMult", i, j, x, nrows, ncols, b, PACKAGE="DECIPHER")
SEXP sparseMult(SEXP row, SEXP col, SEXP value, SEXP nrows, SEXP ncols, SEXP b)
{
	int i, j;
	int *rows = INTEGER(row);
	int *cols = INTEGER(col);
	double *values = REAL(value);
	R_len_t l = length(row);
	unsigned long int r = asInteger(nrows);
	unsigned long int n = asInteger(ncols);
	int reps = length(b)/n;
	double *B = REAL(b);
	
	double *rans;
	SEXP ans;
	//PROTECT(ans = allocVector(REALSXP, r));
	PROTECT(ans = allocMatrix(REALSXP, r, reps));
	rans = REAL(ans);
	
	// zero the result
	for (i = 0; i < r*reps; i++)
		*(rans + i) = 0;
	
	// calculate A %*% x
	for (j = 0; j < reps; j++) {
		int startIndex1 = j*n;
		int startIndex2 = j*r;
		for (i = 0; i < l; i++)
			*(rans + (*(rows + i) - 1) + startIndex2) += *(values + i) * *(B + (*(cols + i) - 1) + startIndex1);
	}
	
	UNPROTECT(1);
	
	return ans;
}
