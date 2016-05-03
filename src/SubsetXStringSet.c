/****************************************************************************
 *           Returns the Specified Subset of an XStringSet Object           *
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

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// DECIPHER header file
#include "DECIPHER.h"

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"

// strcpy
#include <string.h>

SEXP subsetXStringSet(SEXP x, SEXP subset, SEXP type, SEXP nThreads)
{
	int i, *width;
	SEXP ans_width, ans;
	int t = asInteger(type);
	int *s = INTEGER(subset);
	int nthreads = asInteger(nThreads);
	int l = length(subset);
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_XStringSet_xsbaseclassname(x);
	
	// determine the length of the XStringSet
	XStringSet_holder x_set, ans_holder;
	x_set = hold_XStringSet(x);
	
	// determine the widths of the XStringSet
	PROTECT(ans_width = NEW_INTEGER(l));
	Chars_holder x_s;
	for (i = 0, width = INTEGER(ans_width); i < l; i++, width++) {
		x_s = get_elt_from_XStringSet_holder(&x_set, *(s + i) - 1);
		*width = x_s.length;
	}
	
	// set the class of the XStringSet
	char ans_classname[40];
	if (t==1) {
		strcpy(ans_classname, "DNAStringSet");
	} else if (t==2) {
		strcpy(ans_classname, "RNAStringSet");
	} else { // t==3
		strcpy(ans_classname, "AAStringSet");
	}
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	Chars_holder ans_elt_holder;
	
	#pragma omp parallel for private(i,ans_elt_holder,x_s) schedule(guided) num_threads(nthreads)
	for (i = 0; i < l; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_s = get_elt_from_XStringSet_holder(&x_set, *(s + i) - 1);
		memcpy((char *) ans_elt_holder.ptr, x_s.ptr, ans_elt_holder.length * sizeof(char));
	}
	
	UNPROTECT(2);
	return ans;
}
