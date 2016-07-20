/****************************************************************************
 *               Provides a Character Vector of Pool Addresses              *
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
SEXP getPools(SEXP x)
{
	int l = length(x);
	char buf[40];
	SEXP s, ans;
	
	PROTECT(ans = allocVector(STRSXP, l));
	
	for (int i = 0; i < l; i++) {
		s = VECTOR_ELT(x, i);
		snprintf(buf, sizeof(buf), "%p", s);
		SET_STRING_ELT(ans, i, mkChar(buf));
	}
	
	UNPROTECT(1);
	
	return ans;
}
