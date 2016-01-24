/****************************************************************************
 *                       Obtain Ordering of a Vector                        *
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

// DECIPHER header file
#include "DECIPHER.h"

SEXP radixOrder(SEXP x, SEXP startAt)
{	
	int i, k, o1;
	int l = length(x);
	int *v = INTEGER(x);
	int s = asInteger(startAt); // starting index
	
	int size = sizeof(int)*CHAR_BIT; // unit size of memory
	int R = size/2; // size of radix key
	int o2 = size - R; // right shift
	int count = 1 << R; // 2^R
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l));
	int *rans = INTEGER(ans);
	for (i = 0; i < l; i++)
		rans[i] = i;
	
	unsigned int shifted;
	
	// least significant digit first
	for (k = 1; k <= (size/R); k++) {
		// subset bits in Radix k
		o1 = size - k*R;
		int *counts = Calloc(count, int); // initialized to zero
		for (i = 0; i < l; i++) {
			shifted = v[rans[i]] << o1; // clear left bits
			shifted >>= o2; // clear right bits
			counts[shifted]++;
		}
		
		// cumulative sum
		for (i = 1; i < count; i++)
			counts[i] = counts[i - 1] + counts[i];
		// shift values right by one
		for (i = count - 1; i > 0; i--)
			counts[i] = counts[i - 1];
		counts[0] = 0;
		
		// move orders
		int *temp = Calloc(l, int);
		for (i = 0; i < l; i++) {
			shifted = v[rans[i]] << o1; // clear left bits
			shifted >>= o2; // clear right bits
			temp[counts[shifted]++] = rans[i];
		}
		
		// replace orders
		for (i = 0; i < l; i++)
			rans[i] = temp[i];
		
		Free(counts);
		Free(temp);
	}
	
	// reorder on sign bit
	for (i = 0; i < l; i++) {
		if (v[rans[i]] < 0)
			break;
	}
	if (i < l) {
		int *temp = Calloc(l, int);
		
		// move negatives before positives
		for (k = i, count = 0; k < l; k++, count++)
			temp[count] = rans[k];
		for (k = 0; k < i; k++, count++)
			temp[count] = rans[k];
		
		if (s != 0)
			for (i = 0; i < l; i++)
				rans[i] = temp[i] + s;
		
		Free(temp);
	} else if (s != 0) {
		for (i = 0; i < l; i++)
			rans[i]++;
	}
	
	UNPROTECT(1);
	
	return ans;
}
