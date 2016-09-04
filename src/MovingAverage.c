/****************************************************************************
 *             Calculates Bounds of an Exponential Moving Average           *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

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

SEXP movAvg(SEXP x, SEXP type, SEXP alpha, SEXP thresh, SEXP maxAvg, SEXP start, SEXP end)
{
	int i, j, m;
	
	double a = asReal(alpha);
	double b = 1 - a;
	double t = asReal(thresh);
	t *= 2; // moving average is doubled
	int k = asInteger(type);
	double mA = asReal(maxAvg);
	
	SEXP lefts, rights;
	PROTECT(lefts = duplicate(start));
	int *left = INTEGER(lefts);
	PROTECT(rights = duplicate(end));
	int *right = INTEGER(rights);
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	int l = get_length_from_XStringSet_holder(&x_set);
	
	for (i = 0; i < l; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		m = x_i.length;
		if (m==0)
			continue;
		
		// initialize array of error probabilities
		double *p = Calloc(m, double); // initialized to zero
		
		// initialize arrays of weighted averages
		double *s1 = Calloc(m, double); // initialized to zero
		double *s2 = Calloc(m, double); // initialized to zero
		
		if (k==1) { // Phred
			for (j = 0; j < m; j++)
				p[j] = pow(10, ((double)x_i.ptr[j] - 33)/-10);
		} else if (k==2) { // Solexa
			for (j = 0; j < m; j++)
				p[j] = 1 - 1/(1 + pow(10, ((double)x_i.ptr[j] - 64)/-10));
		} else { // Illumina
			for (j = 0; j < m; j++)
				p[j] = pow(10, ((double)x_i.ptr[j] - 64)/-10);
		}
		
		// trailing moving average
		s1[0] = p[0];
		for (j = 1; j < m; j++)
			s1[j] = a*p[j] + b*s1[j - 1];
		
		// leading moving average
		s2[m - 1] = p[m - 1];
		for (j = m - 2; j >= 0; j--)
			s2[j] = a*p[j] + b*s2[j + 1];
		
		// combined moving average
		for (j = 0; j < m; j++)
			s1[j] += s2[j];
		
		Free(s2);
		
		// find the longest region below threshold
		int longest = 0;
		int temp = 0;
		int lastStart = left[i];
		int bestEnd = -2;
		for (j = left[i] - 1; j < right[i]; j++) {
			if (s1[j] <= t) {
				if (temp==0)
					lastStart = j + 1;
				temp++;
				if (temp > longest) {
					longest = temp;
					left[i] = lastStart;
					bestEnd = j;
				}
			} else {
				temp = 0;
			}
		}
		
		Free(s1);
		
		right[i] = bestEnd + 1;
		if (longest==0) {
			left[i] = 0;
		} else {
			double avgError = 0;
			for (j = left[i] - 1; j <= bestEnd; j++)
				avgError += p[j];
			avgError /= bestEnd - left[i] + 2;
			if (avgError > mA) {
				left[i] = 0;
				right[i] = -1;
			}
		}
		
		Free(p);
	}
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret_list, 0, lefts);
	SET_VECTOR_ELT(ret_list, 1, rights);
	
	UNPROTECT(3);
	
	return ret_list;
}
