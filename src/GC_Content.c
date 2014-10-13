/****************************************************************************
 *                     Calculates GC Content In Windows                     *
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("gcContent", myDNAString, begins, ends, PACKAGE="DECIPHER")
SEXP gcContent(SEXP x, SEXP begins, SEXP ends)
{
	int i, j, l = length(begins);
	int *b = INTEGER(begins), *e = INTEGER(ends);
	double bits, tot;
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, l));
	double *rans = REAL(ans);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	
	for (i = 0; i < l; i++) {
		tot = 0;
		bits = 0;
		for (j = *(b + i) - 1; j < *(e + i); j++) {
			tot++;
			switch (x_i.seq[j]) {
				case 2: // C
					bits++;
					break;
				case 3: // M
					bits += 0.5;
					break;
				case 4: // G
					bits++;
					break;
				case 5: // R
					bits += 0.5;
					break;
				case 6: // S
					bits++;
					break;
				case 7: // V
					bits += 2/3;
					break;
				case 10: // Y
					bits += 0.5;
					break;
				case 11: // H
					bits += 1/3;
					break;
				case 12: // K
					bits += 0.5;
					break;
				case 13: // D
					bits += 1/3;
					break;
				case 14: // B
					bits += 2/3;
					break;
				case 15: // N
					bits += 0.5;
					break;
				case 16: // -
					tot--;
					break;
				case 32: // +
					tot--;
					break;
				case 64: // .
					tot--;
					break;
				default: // any A/T
					break;
			}
		}
		*(rans + i) = bits/tot;
	}
	
	UNPROTECT(1);
	
	return ans;
}
