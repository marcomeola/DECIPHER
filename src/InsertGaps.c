/****************************************************************************
 *                  Inserts Gaps into Aligned XStringSets                   *
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

//ans_start <- .Call("insertGaps", sequences, positions, lengths, processors, PACKAGE="DECIPHER")
SEXP insertGaps(SEXP x, SEXP positions, SEXP lengths, SEXP type, SEXP nThreads)
{
	int i, j, x_length, *width, sum, start;
	SEXP ans_width, ans;
	int t = asInteger(type);
	int *p = INTEGER(positions);
	int *l = INTEGER(lengths);
	int n = length(positions);
	int nthreads = asInteger(nThreads);
	
	// determine the cumulative width of insertions
	sum = 0;
	for (i = 0; i < n; i++)
		sum += l[i];
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_XStringSet_xsbaseclassname(x);
	
	// determine the length of the XStringSet
	XStringSet_holder x_set, ans_holder;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	// determine the widths of the aligned (equal width) XStringSet
	PROTECT(ans_width = NEW_INTEGER(x_length));
	Chars_holder x_s;
	x_s = get_elt_from_XStringSet_holder(&x_set, 0);
	for (i = 0, width = INTEGER(ans_width); i < x_length; i++, width++)
		*width = x_s.length + sum;
	
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
	
	#pragma omp parallel for private(i,j,ans_elt_holder,x_s,sum,start) schedule(guided) num_threads(nthreads)
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		//ans_elt_holder.length = 0;
		x_s = get_elt_from_XStringSet_holder(&x_set, i);
		sum = 0; // position in ans_elt_holder.ptr
		start = 0; // position in x_s.ptr
		for (j = 0; j < n; j++) {
			if ((p[j] - 1) > start) { // copy over sequence
				memcpy((char *) ans_elt_holder.ptr + sum, x_s.ptr + start, (p[j] - 1 - start) * sizeof(char));
				sum += (p[j] - 1 - start);
				start += (p[j] - 1 - start);
			}
			if (l[j] > 0) { // insert gaps
				if (t==3) { // AAStringSet
					memset((char *) ans_elt_holder.ptr + sum, 45, l[j] * sizeof(char));
				} else { // DNAStringSet or RNAStringSet
					memset((char *) ans_elt_holder.ptr + sum, 16, l[j] * sizeof(char));
				}
				sum += l[j];
			}
		}
		if (sum < ans_elt_holder.length) {
			memcpy((char *) ans_elt_holder.ptr + sum, x_s.ptr + start, (ans_elt_holder.length - sum) * sizeof(char));
		}
	}
	
	UNPROTECT(2);
	return ans;
}
