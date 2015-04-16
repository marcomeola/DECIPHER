/****************************************************************************
 *                      Expands Sequence Ambiguities                        *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

SEXP expandAmbiguities(SEXP x, SEXP c)
{
	int i, j, k, p;
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length;
	SEXP ans, ret_list;
	const char *C = CHAR(STRING_ELT(c, 0));
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	PROTECT(ret_list = allocVector(VECSXP, x_length));
	
	// loop through each sequence in the DNAStringSet
	for (i = 0; i < x_length; i++) {
		// extract each ith DNAString from the DNAStringSet
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// make an array of position, ambiguity
		int pos[x_i.length][2];
		int count = 0;
		int maxCombos = 1;
		
		// determine the non-gap pos and ambiguities
		for (j = 0; j < x_i.length; j++) {
			if (x_i.ptr[j] < 16) {
				pos[count][0] = j;
				switch (x_i.ptr[j]) {
					case 3:
					case 5:
					case 6:
					case 9:
					case 10:
					case 12:
						pos[count][1] = 2;
						break;
					case 7:
					case 11:
					case 13:
					case 14:
						pos[count][1] = 3;
						break;
					case 15:
						pos[count][1] = 4;
						break;
					default:
						pos[count][1] = 1;
						break;
				}
				
				maxCombos *= pos[count][1];
				count++;
			}
		}
		
		int *pers = Calloc(maxCombos*count, int); // initialized to zero
		int alt; // whether to alternate
		int combos = 1; // number of ambiguities
		for (j = 0; j < count; j++) {
			if ((pos[j][1] > 1) && // ambiguity AND
				(!((combos & 1) ^ // XOR
				   (pos[j][1] & 1)) || // both odd or even OR
				 ((combos % pos[j][1]) == 0))) { // remainder is zero
				// alternate order every combo
				alt = 1;
			} else {
				alt = 0;
			}
			
			switch (x_i.ptr[pos[j][0]]) {
				case 1: // A
					for (k = 0; k < maxCombos; k++) {
						pers[k*count + j] = 0;
					}
					break;
				case 2: // C
					for (k = 0; k < maxCombos; k++) {
						pers[k*count + j] = 1;
					}
					break;
				case 3: // M = A or C
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = (k + (int)(k/combos)) & 1;
						} else {
							pers[k*count + j] = k & 1;
						}
					}
					break;
				case 4: // G
					for (k = 0; k < maxCombos; k++) {
						pers[k*count + j] = 2;
					}
					break;
				case 5: // R = A or G
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = ((k + (int)(k/combos)) & 1)*2;
						} else {
							pers[k*count + j] = (k & 1)*2;
						}
					}
					break;
				case 6: // S = C or G
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = ((k + (int)(k/combos)) & 1) + 1;
						} else {
							pers[k*count + j] = (k & 1) + 1;
						}
					}
					break;
				case 7: // V = A or C or G
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = (k + (int)(k/combos)) % 3;
						} else {
							pers[k*count + j] = k % 3;
						}
					}
					break;
				case 8: // T
					for (k = 0; k < maxCombos; k++) {
						pers[k*count + j] = 3;
					}
					break;
				case 9: // W = A or T
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = ((k + (int)(k/combos)) & 1)*3;
						} else {
							pers[k*count + j] = (k & 1)*3;
						}
					}
					break;
				case 10: // Y = C or T
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = ((k + (int)(k/combos)) & 1)*2 + 1;
						} else {
							pers[k*count + j] = (k & 1)*2 + 1;
						}
					}
					break;
				case 11: // H = A or C or T
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = (k + (int)(k/combos)) % 3;
						} else {
							pers[k*count + j] = k % 3;
						}
						if (pers[k*count + j]==2)
							pers[k*count + j]++;
					}
					break;
				case 12: // K = G or T
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = ((k + (int)(k/combos)) & 1) + 2;
						} else {
							pers[k*count + j] = (k & 1) + 2;
						}
					}
					break;
				case 13: // D = A or G or T
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = (k + (int)(k/combos)) % 3;
						} else {
							pers[k*count + j] = k % 3;
						}
						if (pers[k*count + j] > 0)
							pers[k*count + j]++;
					}
					break;
				case 14: // B = C or G or T
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = ((k + (int)(k/combos)) % 3) + 1;
						} else {
							pers[k*count + j] = (k % 3) + 1;
						}
					}
					break;
				case 15: // N = A or C or G or T
					for (k = 0; k < maxCombos; k++) {
						if (alt) {
							pers[k*count + j] = (k + (int)(k/combos)) & 3;
						} else {
							pers[k*count + j] = k & 3;
						}
					}
					break;
			}
			combos *= pos[j][1];
		}
		
		PROTECT(ans = allocVector(STRSXP, combos));
		for (j = 0; j < maxCombos; j++) {
			char perm[x_i.length + 1]; // last position is for null terminating
			p = 0;
			for (k = 0; k < x_i.length; k++) {
				if (x_i.ptr[k] < 16) {
					switch (pers[j*count + p]) {
						case 0:
							perm[k] = 'A';
							break;
						case 1:
							perm[k] = 'C';
							break;
						case 2:
							perm[k] = 'G';
							break;
						case 3:
							perm[k] = C[0];
							break;
					}
					p++;
				} else {
					switch (x_i.ptr[k]) {
						case 16:
							perm[k] = '-';
							break;
						case 32:
							perm[k] = '+';
							break;
						case 64:
							perm[k] = '.';
							break;
						default:
							error("not DNA!");
							break;
					}
				}
			}
			perm[k] = '\0'; // end (null terminate) the string
			SET_STRING_ELT(ans, j, mkChar(perm));
		}
		
		Free(pers);
		
		SET_VECTOR_ELT(ret_list, i, ans);
		UNPROTECT(1);
	}
	
	UNPROTECT(1);
	
	return ret_list;
}
