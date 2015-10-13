/****************************************************************************
 *         Helper Functions for Assigning Insertions and Deletions          *
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

// DECIPHER header file
#include "DECIPHER.h"

// clears "ins" attribute in-place within a nested list
SEXP clearIns(SEXP x)
{	
	if (length(x) > 1) { // not a leaf
		clearIns(VECTOR_ELT(x, 0));
		clearIns(VECTOR_ELT(x, 1));
	}
	
	// clear "ins" attribute
	setAttrib(x, install("ins"), R_NilValue);
	
	return R_NilValue;
}

// returns TRUE when all are true and NA if all are NA
SEXP all(SEXP x)
{	
	int i, l = length(x), *v = INTEGER(x), count = 0;
	for (i = 0; i < l; i++) {
		if (v[i]!=NA_LOGICAL) {
			if (count==0)
				count = 1;
			if (v[i]==0)
				break;
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(LGLSXP, 1));
	int *rans = INTEGER(ans);
	
	if (count==0) {
		rans[0] = NA_LOGICAL;
	} else if (i < l) {
		rans[0] = 0;
	} else {
		rans[0] = 1;
	}
	
	UNPROTECT(1);
	
	return ans;
}

// returns TRUE when any are true and NA if all are NA
SEXP any(SEXP x)
{	
	int i, l = length(x), *v = INTEGER(x), count = 0, en = 0;
	for (i = 0; i < l; i++) {
		if (v[i]!=NA_LOGICAL) {
			if (count==0)
				count = 1;
			if (v[i]==1) {
				en = 1;
				break;
			}
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(LGLSXP, 1));
	int *rans = INTEGER(ans);
	
	if (count==0) {
		rans[0] = NA_LOGICAL;
	} else if (en) {
		rans[0] = 1;
	} else {
		rans[0] = 0;
	}
	
	UNPROTECT(1);
	
	return ans;
}
