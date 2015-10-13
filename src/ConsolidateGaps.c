/****************************************************************************
 *                    Consolidates Gaps in XStringSets                      *
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"

SEXP consolidateGaps(SEXP x, SEXP type)
{
	int i, j, k, l, mp, x_length, count, c;
	int t = asInteger(type);
	char p;
	
	// determine the length of the XStringSet
	XStringSet_holder x_set;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	Chars_holder x_i;
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	l = x_i.length; // assume sequences are aligned
	
	// pos keeps track of last (non-gap) position
	int *pos = (int *) R_alloc(x_length, sizeof(int));
	for (i = 0; i < x_length; i++) {
		pos[i] = NA_INTEGER;
	}
	
	// w keeps track of which sequences are non-gap
	int *w = (int *) R_alloc(x_length, sizeof(int));
	
	for (i = 0; i < l; i++) {
		count = 0;
		mp = 0; // max pos
		for (j = 0; j < x_length; j++) {
			x_i = get_elt_from_XStringSet_holder(&x_set, j);
			
			if ((t==3 && x_i.ptr[i] ^ 0x2D && x_i.ptr[i] ^ 0x2E) ||
				(t!=3 && !(x_i.ptr[i] & 0x10 || x_i.ptr[i] & 0x40))) {
				// not a gap ("-") or unknown (".") in this position
				w[count] = j;
				count++;
				
				if (pos[j]==NA_INTEGER) // pos is NA
					pos[j] = i;
				
				if (pos[j] > mp)
					mp = pos[j];
			}
		}
		
		if (count==0) // column is all gaps
			continue;
		
		if ((i - mp) > 1) { // swap columns
			c = 0;
			for (j = 0; j < x_length; j++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, j);
				
				if (c < count && j==w[c]) { // shift to left
					c++;
					p = x_i.ptr[mp + 1];
					*((char *)x_i.ptr + mp + 1) = *((char *)x_i.ptr + i);
					*((char *)x_i.ptr + i) = p;
					pos[j] = mp + 1;
				} else { // shift to right
					p = x_i.ptr[i];
					for (k = i; k >= (mp + 2); k--) {
						*((char *)x_i.ptr + k) = *((char *)x_i.ptr + k - 1);
					}
					*((char *)x_i.ptr + mp + 1) = p;
					
					if (pos[j] >= mp + 1)
						pos[j]++;
				}
			}
		} else {
			for (j = 0; j < count; j++)
				pos[w[j]] = i;
		}
	}
	
	return R_NilValue;
}
