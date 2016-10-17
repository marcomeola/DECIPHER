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
SEXP intMatch(SEXP x, SEXP y, SEXP nThreads)
{	
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int i, j;
	int size_x = length(x);
	int size_y = length(y);
	int nthreads = asInteger(nThreads);
	
	SEXP ans;
	PROTECT(ans = allocVector(LGLSXP, size_x));
	int *rans = INTEGER(ans);
	
	#pragma omp parallel for private(i, j) schedule(guided) num_threads(nthreads)
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
SEXP firstMatchUpper(SEXP x, SEXP y, SEXP nThreads)
{	
	int i, j, size_x = length(x), size_y = length(y);
	double *v = REAL(x);
	double *w = REAL(y);
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, size_x));
	int *rans = INTEGER(ans);
	int nthreads = asInteger(nThreads);
	
	#pragma omp parallel for private(i, j) schedule(guided) num_threads(nthreads)
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
SEXP matchLists(SEXP x, SEXP verbose, SEXP pBar, SEXP nThreads)
{	
	int i, j, size_x = length(x), before, v, *rPercentComplete;
	int o, p, start, lx, ly, *X, *Y, count;
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size_x, size_x));
	double *rans = REAL(ans);
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	int nthreads = asInteger(nThreads);
	
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
		#pragma omp parallel for private(j, o, p, start, count, X, Y, lx, ly) schedule(guided) num_threads(nthreads)
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
				
				if (lx < ly) {
					*(rans + i*size_x + j) = 1 - (double)count/(double)lx;
				} else {
					*(rans + i*size_x + j) = 1 - (double)count/(double)ly;
				}
			} else {
				*(rans + i*size_x + j) = NA_REAL;
			}
			*(rans + i + j*size_x) = *(rans + i*size_x + j);
		}
		
		if (v) {
			// print the percent completed so far
			//*rPercentComplete = floor(100*(double)((i + 1)*size_x+(i + 1))/((size_x - 1)*size_x+(size_x - 1)));
			*rPercentComplete = floor(100*(double)(2*size_x - 2 - i)*(i + 1)/((size_x - 1)*size_x));
			
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

// matrix of d[i, j] = length x[i] %in% y[j] / min(length)
// requires a list of ordered integers
SEXP matchListsDual(SEXP x, SEXP y, SEXP verbose, SEXP pBar, SEXP nThreads)
{	
	int i, j, size_x = length(x), size_y = length(y), before, v, *rPercentComplete;
	int o, p, start, lx, ly, *X, *Y, count;
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size_x, size_y));
	double *rans = REAL(ans);
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	int nthreads = asInteger(nThreads);
	
	if (v) { // initialize progress variables
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	for (i = 0; i < size_x; i++) {
		#pragma omp parallel for private(j, o, p, start, count, X, Y, lx, ly) schedule(guided) num_threads(nthreads)
		for (j = 0; j < size_y; j++) {
			X = INTEGER(VECTOR_ELT(x, i));
			Y = INTEGER(VECTOR_ELT(y, j));
			lx = length(VECTOR_ELT(x, i));
			ly = length(VECTOR_ELT(y, j));
			
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
					*(rans + j*size_x + i) = NA_REAL;
					continue;
				}
				
				for (o = lx - 1; o >= 0; o--) {
					if (X[o] <= Y[ly - 1]) {
						last = o;
						break;
					}
				}
				if (last == -1) { // no overlap
					*(rans + j*size_x + i) = NA_REAL;
					continue;
				}
				
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
				
				if (lx > ly) {
					*(rans + j*size_x + i) = (double)count/(double)ly;
				} else {
					*(rans + j*size_x + i) = (double)count/(double)lx;
				}
			} else {
				*(rans + j*size_x + i) = NA_REAL;
			}
		}
		
		if (v) {
			// print the percent completed so far
			*rPercentComplete = floor(100*(double)(i + 1)/size_x);
			
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

// matrix of d[i, j] = ordered matches (x[i], x[j]) / min(length)
// requires a list of unsorted integers retaining the order of x
SEXP matchOrder(SEXP x, SEXP verbose, SEXP pBar, SEXP nThreads)
{	
	int i, j, size_x = length(x), before, v, *rPercentComplete;
	int lx, ly, *X, *Y;
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size_x, size_x));
	double *rans = REAL(ans);
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	int nthreads = asInteger(nThreads);
	
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
		#pragma omp parallel for private(j, X, Y, lx, ly) schedule(guided) num_threads(nthreads)
		for (j = i + 1; j < size_x; j++) {
			X = INTEGER(VECTOR_ELT(x, i));
			Y = INTEGER(VECTOR_ELT(x, j));
			lx = length(VECTOR_ELT(x, i));
			ly = length(VECTOR_ELT(x, j));
			
			int link_x = -1, link_y = -1; // last established link between X and Y
			int matches = 0; // running number of matches
			int pos_x, pos_y, off; // current position
			int offset = 1; // distance offset from link
			int delta, forward;
			while (offset + link_x < lx && offset + link_y < ly) { // within sequence
				pos_y = link_y + 1;
				pos_x = link_x + offset;
				
				if (matches) {
					if (forward) {
						for (off = 1; off <= offset; off++, pos_x--, pos_y++) {
							if (X[pos_x]==Y[pos_y]) {
								link_x = pos_x;
								link_y = pos_y;
								offset = 0;
								matches++;
							}
						}
					} else { // backward
						for (off = 1; off <= offset; off++, pos_x--, pos_y++) {
							if (X[lx - pos_x - 1]==Y[ly - pos_y - 1]) {
								link_x = pos_x;
								link_y = pos_y;
								offset = 0;
								matches++;
							}
						}
					}
				} else {
					for (off = 1; off <= offset; off += delta, pos_x -= delta, pos_y += delta) {
						if (X[pos_x]==Y[pos_y]) {
							link_x = pos_x;
							link_y = pos_y;
							offset = 0;
							matches++;
							forward = 1;
							break;
						} else if (X[lx - pos_x - 1]==Y[ly - pos_y - 1]) {
							link_x = pos_x;
							link_y = pos_y;
							offset = 0;
							matches++;
							forward = 0; // backward
							break;
						}
						delta = (off < 10) ? 1 : (off/5);
					}
				}
				
				offset++;
			}
			
			if (lx < ly) {
				*(rans + i*size_x + j) = 1 - (double)matches/(double)lx;
			} else {
				*(rans + i*size_x + j) = 1 - (double)matches/(double)ly;
			}
			
			*(rans + i + j*size_x) = *(rans + i*size_x + j);
		}
		
		if (v) {
			// print the percent completed so far
			//*rPercentComplete = floor(100*((double)i/((double)size_x - 1)));
			//*rPercentComplete = floor(100*(double)((i + 1)*size_x+(i + 1))/((size_x - 1)*size_x+(size_x - 1)));
			*rPercentComplete = floor(100*(double)(2*size_x - 2 - i)*(i + 1)/((size_x - 1)*size_x));
			
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

// matrix of d[i, j] = ordered matches (x[i], y[j]) / min(length bound by consecutive matches)
// requires dual lists of unsorted integers retaining the order of x (typically length(x) < length(y))
SEXP matchOrderDual(SEXP x, SEXP y, SEXP nThreads)
{	
	int i, j, size_x = length(x), size_y = length(y);
	int lx, ly, *X, *Y;
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size_x, size_y));
	double *rans = REAL(ans);
	int nthreads = asInteger(nThreads);
	
	for (i = 0; i < size_x; i++) {
		#pragma omp parallel for private(j, X, Y, lx, ly) schedule(guided) num_threads(nthreads)
		for (j = 0; j < size_y; j++) {
			X = INTEGER(VECTOR_ELT(x, i));
			Y = INTEGER(VECTOR_ELT(y, j));
			lx = length(VECTOR_ELT(x, i));
			ly = length(VECTOR_ELT(y, j));
			
			int link_x = -1, link_y = -1; // last established link between X and Y
			int matches = 0; // running number of matches
			int pos_x, pos_y, off; // current position
			int offset = 1; // distance offset from link
			while (offset <= (lx - link_x + ly - link_y - 2)) { // within sequence
				for (off = 1; off <= offset; off++) {
					pos_y = link_y + off;
					pos_x = link_x + offset - off + 1;
					if (pos_y < ly &&
						pos_x < lx &&
						X[pos_x]==Y[pos_y]) {
						link_x = pos_x;
						link_y = pos_y;
						offset = 0;
						matches++;
					}
				}
				offset++;
			}
			
			if (lx < ly) {
				*(rans + i + j*size_x) = 1 - (double)matches/(double)(lx);
			} else {
				*(rans + i + j*size_x) = 1 - (double)matches/(double)(ly);
			}
			
			link_x = -1, link_y = -1; // last established link between X and Y
			matches = 0; // running number of matches
			offset = 1; // distance offset from link
			while (offset <= (lx - link_x + ly - link_y - 2)) { // within sequence
				for (off = 1; off <= offset; off++) {
					pos_y = link_y + off;
					pos_x = link_x + offset - off + 1;
					if (pos_y < ly &&
						pos_x < lx &&
						X[lx - pos_x - 1]==Y[ly - pos_y - 1]) {
						link_x = pos_x;
						link_y = pos_y;
						offset = 0;
						matches++;
					}
				}
				offset++;
			}
			
			double temp;
			if (lx < ly) {
				temp = 1 - (double)matches/(double)(lx);
			} else {
				temp = 1 - (double)matches/(double)(ly);
			}
			if (temp < *(rans + i + j*size_x))
				*(rans + i + j*size_x) = temp;
		}
		
		R_CheckUserInterrupt();
	}
	
	UNPROTECT(1);
	
	return ans;
}

// returns shared ranges between pairs in two unordered lists
SEXP matchRanges(SEXP x, SEXP y, SEXP wordSize, SEXP maxLength, SEXP threshold)
{	
	int i, j, size_x = length(x), size_y = length(y), size;
	int lx, ly, *X, *Y, *X_pos, *Y_pos, wS, *rans;
	int l = asInteger(maxLength);
	double thresh = asReal(threshold);
	wS = asInteger(wordSize);
	int *bits = Calloc(l*2, int); // initialized to zero
	
	if (size_x > size_y) {
		size = size_x;
	} else {
		size = size_y;
	}
	
	for (i = 0; i < size; i += 2) {
		j = i;
		X = INTEGER(VECTOR_ELT(x, i));
		Y = INTEGER(VECTOR_ELT(y, j));
		X_pos = INTEGER(VECTOR_ELT(x, i + 1));
		Y_pos = INTEGER(VECTOR_ELT(y, j + 1));
		lx = length(VECTOR_ELT(x, i));
		ly = length(VECTOR_ELT(y, j));
		
		int link_x = -1, link_y = -1; // last established link between X and Y
		int pos_x, pos_y; // current position
		
		int offset = 1; // distance offset from link
		while (offset <= (lx - link_x + ly - link_y - 2)) { // within sequence
			for (int off = 1; off <= offset; off++) {
				pos_y = link_y + off;
				pos_x = link_x + offset - off + 1;
				if (pos_y < ly &&
					pos_x < lx &&
					X[pos_x]==Y[pos_y] &&
					X[pos_x] != NA_INTEGER) {
					if (pos_x==(link_x + 1) &&
						pos_y==(link_y + 1)) {
						if (*(bits + l + X_pos[pos_x] - 1)==0) { // new anchor
							*(bits + l + X_pos[pos_x] - 1) = Y_pos[pos_y]; // anchor position
							*(bits + X_pos[pos_x] - 1) += 1; // increment anchoring
						} else if (*(bits + l + X_pos[pos_x] - 1)==Y_pos[pos_y]) { // previous anchor
							*(bits + X_pos[pos_x] - 1) += 1; // increment anchoring
						} else { // does not match previous anchoring
							*(bits + X_pos[pos_x] - 1) = 0; // reset inconsistent anchoring
							*(bits + l + X_pos[pos_x] - 1) = 0;
						}
					}
					link_x = pos_x;
					link_y = pos_y;
					offset = 0;
				}
			}
			offset++;
		}
		
		R_CheckUserInterrupt();
	}
	
	size /= 2;
	// calculate ranges of anchors
	int *temp = Calloc(l, int); // initialized to zero
	int match = 0, count = -1, last_end_x = -1, last_end_y = -1;
	for (i = 0; i < l; i++) {
		//Rprintf("\n%d start_x=%d start_y=%d end_x=%d end_y=%d last_end_x=%d last_end_y=%d", i, i - wS + 2, *(bits + i + l) - wS + 1, i+1, *(bits + i + l), last_end_x, last_end_y);
		if ((double)*(bits + i)/(double)size >= thresh) {
			if (match==0) { // start of anchor range
				if (i - wS + 2 > last_end_x + 100 &&
					*(bits + i + l) - wS + 1 > last_end_y + 100) {
					match = 1;
					count++;
					*(temp + count*4) = i + 1;//i - wS + 2;
					*(temp + count*4 + 2) = *(bits + i + l);//*(bits + i + l) - wS + 1;
					*(temp + count*4 + 1) = i + 1;
					*(temp + count*4 + 3) = *(bits + i + l);
				}
			} else if (i - wS + 2 > last_end_x + 1000 &&
				*(bits + i + l) - wS + 1 > last_end_y + 1000) {
				count++; // start new range
				last_end_x = i - 100 - wS;
				last_end_y = *(bits + i + l - 1) - 100 - wS;
				*(temp + count*4) = i + 1;//i - wS + 2;
				*(temp + count*4 + 2) = *(bits + i + l);//*(bits + i + l) - wS + 1;
				*(temp + count*4 + 1) = i + 1;
				*(temp + count*4 + 3) = *(bits + i + l);
			} else if (i + 1 > last_end_x + 100 &&
				*(bits + i + l) > last_end_y + 100 &&
				*(bits + i + l) > *(temp + count*4 + 3)) { // extend end
				*(temp + count*4 + 1) = i + 1;
				*(temp + count*4 + 3) = *(bits + i + l);
			}
		} else if (match==1) { // set to last match
			match = 0;
			last_end_x = *(temp + count*4 + 1);
			last_end_y = *(temp + count*4 + 3);
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, 4, count + 1));
	rans = INTEGER(ans);
	for (i = 0; i <= count; i++) {
		*(rans + i*4) = *(temp + i*4);
		*(rans + i*4 + 1) = *(temp + i*4 + 1);
		*(rans + i*4 + 2) = *(temp + i*4 + 2);
		*(rans + i*4 + 3) = *(temp + i*4 + 3);
	}
	
	UNPROTECT(1);
	Free(bits);
	Free(temp);
	
	return ans;
}

// which (bl <= x <= bu)
// x must be in ascending order
SEXP boundedMatches(SEXP x, SEXP bl, SEXP bu)
{	
	int i, mid, start = 0, end = length(x), count = 0, size_x = length(x);
	int lowBound = asInteger(bl);
	int upBound = asInteger(bu);
	int *v = INTEGER(x);
	int *buffer = (int *) R_alloc(size_x, sizeof(int));
	
	while (start < end) {
		mid = floor(start + (end - start)/2);
		if (v[mid] >= lowBound) {
			end = mid;
		} else if (start == mid) {
			break;
		} else {
			start = mid;
		}
	}
	
	for (i = end; i < size_x; i++) {
		if (v[i] >= lowBound && v[i] <= upBound) {
			buffer[count] = i + 1;
			count++;
		} else {
			break;
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, count));
	memcpy(INTEGER(ans), buffer, sizeof(int) * count);
	
	UNPROTECT(1);
	
	return ans;
}

// first unmatched occurrence of x in y for ascending order integer vectors
// requires NAs to be first in the ordering (treated as the most negative)
// similar to match(x, y, incomparables=NA) without repetition in the result
// performs the following reorderings of the inputs and output:
// ans[o1] <- o2[intMatchOnce(x[o1], y[o2])] + 1
// where o1 and o2 are indexed starting at zero, and ans starts at 1
SEXP intMatchOnce(SEXP x, SEXP y, SEXP o1, SEXP o2)
{	
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *p = INTEGER(o1);
	int *q = INTEGER(o2);
	
	int i, j, k, temp, start = 0;
	int size_x = length(x);
	int size_y = length(y);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, size_x));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < size_x; i++) {
		rans[p[i]] = NA_INTEGER;
		if (v[p[i]] != NA_INTEGER)
			break;
	}
	
	for (; i < size_x; i++) { // i = i
		temp = NA_INTEGER;
		for (j = start; j < size_y; j++) {
			if (v[p[i]] < w[q[j]]) {
				start = j;
				break;
			} else if (v[p[i]] == w[q[j]]) {
				k = j + 1;
				if (k < size_y &&
					w[q[j]] == w[q[k]]) {
					start = k; // prevent repeated matching
					temp = j;
				} else {
					start = j; // allow repeated matching
					temp = j;
				}
				break;
			}
		}
		
		if (temp == NA_INTEGER) {
			rans[p[i]] = NA_INTEGER;
		} else {
			rans[p[i]] = q[temp] + 1;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}
