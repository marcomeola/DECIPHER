/****************************************************************************
 *            Difference Between Consecutive Elements in a Vector           *
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

SEXP intDiff(SEXP x)
{	
	int l = length(x);
	int *v = INTEGER(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l - 1));
	int *rans = INTEGER(ans);
	
	l -= 1;
	for (int i = 0; i < l; i++)
		rans[i] = v[i + 1] - v[i];
	
	UNPROTECT(1);
	
	return ans;
}
