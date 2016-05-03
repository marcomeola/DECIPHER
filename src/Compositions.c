/****************************************************************************
 *                   Calculates Compositions of Sequences                   *
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
#include "XVector_interface.h"

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
			switch (x_i.ptr[j]) {
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

//ans_start <- .Call("composition", x, PACKAGE="DECIPHER")
SEXP composition(SEXP x)
{
	int i, j, l;
	int n = length(x);
	int cD = 0, cR = 0, cA = 0, cB = 0;
	const char *seq;
	
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		
		for (j = 0; j < l; j++) {
			switch (seq[j]) {
				case '.':
				case '-':
				case '+':
					break;
				case 'A':
				case 'a':
				case 'C':
				case 'c':
				case 'G':
				case 'g':
					cD++;
					cR++;
					cA++;
					cB++;
					break;
				case 'T':
				case 't':
					cD++;
					cA++;
					cB++;
					break;
				case 'U':
				case 'u':
					cR++;
					cB++;
					break;
				case 'R':
				case 'r':
				case 'N':
				case 'n':
				case 'D':
				case 'd':
				case 'Q':
				case 'q':
				case 'E':
				case 'e':
				case 'H':
				case 'h':
				case 'I':
				case 'i':
				case 'L':
				case 'l':
				case 'K':
				case 'k':
				case 'M':
				case 'm':
				case 'F':
				case 'f':
				case 'P':
				case 'p':
				case 'S':
				case 's':
				case 'W':
				case 'w':
				case 'Y':
				case 'y':
				case 'V':
				case 'v':
				case 'X':
				case 'x':
				case '*':
					cA++;
					cB++;
					break;
				default:
					cB++;
					break;
			}
		}
	}
	
	double *rans;
	SEXP ans;
	
	PROTECT(ans = allocVector(REALSXP, 3));
	rans = REAL(ans);
	
	if (cB==0)
		cB = 1; // prevent divide by zero
	
	rans[0] = (double)cD/(double)cB;
	rans[1] = (double)cR/(double)cB;
	rans[2] = (double)cA/(double)cB;
	
	UNPROTECT(1);
	
	return ans;
}

SEXP positionWeightMatrix(SEXP x, SEXP begins, SEXP ends, SEXP width)
{
	int i, j, p, l = length(begins);
	int *b = INTEGER(begins), *e = INTEGER(ends);
	int w = asInteger(width);
	
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, 5, w));
	int *rans = INTEGER(ans);
	for (i = 0; i < w*5; i++)
		*(rans + i) = 0;
	
	Chars_holder xstring;
	xstring = hold_XRaw(x);
	
	for (i = 0; i < l; i++) {
		for (j = *(b + i) - 1, p = 0; j < *(e + i); j++, p++) {
			switch (xstring.ptr[j]) {
				case 1: // A
					*(rans + p*5) = *(rans + p*5) + 1;
					break;
				case 2: // C
					*(rans + p*5 + 1) = *(rans + p*5 + 1) + 1;
					break;
				case 4: // G
					*(rans + p*5 + 2) = *(rans + p*5 + 2) + 1;
					break;
				case 8: // T
					*(rans + p*5 + 3) = *(rans + p*5 + 3) + 1;
					break;
				default: // other
					*(rans + p*5 + 4) = *(rans + p*5 + 4) + 1;
					break;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}
