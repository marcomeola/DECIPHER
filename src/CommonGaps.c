/****************************************************************************
 *                           Removes Common Gaps                            *
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

// for math functions
#include <math.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("commonGaps", sequences, PACKAGE="DECIPHER")
SEXP commonGaps(SEXP x)
{
	int i, j, l, p, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	// initialize vector of positions
	int bits[longest];
	for (i = 0; i < longest; i++)
		bits[i] = 0;
	
	// determine which positions are not gaps
	p = 0; // number of non-gap positions
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		for (j = 0; j < l; j++) {
			if (bits[j]==0) {
				if (seq[j]!='-') {
					bits[j] = 1;
					p++;
				}
			}
		}
	}
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char s[p + 1]; // each sequence
	
	// write new character vector
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		count = 0;
		for (j = 0; j < l; j++) {
			if (bits[j]==1) {
				s[count] = seq[j];
				count ++;
			}
		}
		s[count] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(s));
	}
	
	UNPROTECT(1);
	
	return seqs;	
}
