/****************************************************************************
 *               Translate DNA/RNAStringSet into AAStringSet                *
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

// strcpy
#include <string.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

SEXP basicTranslate(SEXP x, SEXP code, SEXP starts)
{
	int i, j, l, pos, index;
	int *s = INTEGER(starts);
	SEXP ans_width, ans;
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_XStringSet_xsbaseclassname(code);
	
	// determine the length of the XStringSet
	XStringSet_holder x_set, code_set, ans_holder;
	x_set = hold_XStringSet(x);
	l = get_length_from_XStringSet_holder(&x_set);
	code_set = hold_XStringSet(code);
	
	// determine the widths of the XStringSet
	PROTECT(ans_width = NEW_INTEGER(l));
	int *width = INTEGER(ans_width);
	Chars_holder x_s;
	for (i = 0; i < l; i++) {
		x_s = get_elt_from_XStringSet_holder(&x_set, i);
		width[i] = (x_s.length - s[i] + 1)/3;
	}
	
	// set the class of the XStringSet
	char ans_classname[40];
	strcpy(ans_classname, "AAStringSet");
	
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	Chars_holder ans_elt_holder;
	
	Chars_holder code_s = get_elt_from_XStringSet_holder(&code_set, 0);
	for (i = 0; i < l; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_s = get_elt_from_XStringSet_holder(&x_set, i);
		
		pos = s[i] - 1;
		for (j = 0; j < width[i]; j++) {
			// first codon position
			if (*((char *)x_s.ptr + pos) & 0x1) {
				index = 0;
			} else if (*((char *)x_s.ptr + pos) & 0x2) {
				index = 16;
			} else if (*((char *)x_s.ptr + pos) & 0x4) {
				index = 32;
			} else {
				index = 48;
			}
			
			// second codon position
			pos++;
			if (*((char *)x_s.ptr + pos) & 0x2) {
				index += 4;
			} else if (*((char *)x_s.ptr + pos) & 0x4) {
				index += 8;
			} else if (*((char *)x_s.ptr + pos) & 0x8) {
				index += 12;
			}
			
			// third codon position
			pos++;
			if (*((char *)x_s.ptr + pos) & 0x2) {
				index += 1;
			} else if (*((char *)x_s.ptr + pos) & 0x4) {
				index += 2;
			} else if (*((char *)x_s.ptr + pos) & 0x8) {
				index += 3;
			}
			
			pos++;
			*((char *)ans_elt_holder.ptr + j) = *((char *)code_s.ptr + index);
		}
	}
	
	UNPROTECT(2);
	return ans;
}
