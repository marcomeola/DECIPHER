/****************************************************************************
 *               Predicts 3 State Secondary Protein Structure               *
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

SEXP predictHEC(SEXP x, SEXP windowSize, SEXP background, SEXP HEC_MI1, SEXP HEC_MI2, SEXP output)
{
	int i, j, l, k1, k2, p1, p2;
	int wS = asInteger(windowSize);
	double f1 = (2*(double)wS - 1)/(2*(double)wS + 1); // fraction of each single
	double f2 = 2/(2*(double)wS + 1); // fraction of each double
	double *bg = REAL(background);
	double *MI1 = REAL(HEC_MI1); // [AA, pos, HEC]
	int o = asInteger(output);
	int total = length(HEC_MI1)/60; // maximum window
	int center = (total - 1)/2; // center of window
	double *MI2 = REAL(HEC_MI2); // [AA1, AA2, pos1, pos2, HEC]
	double H, E, C, sum, *rans;
	char *states;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	int x_length = get_length_from_XStringSet_holder(&x_set);
	
	SEXP ret, ans;
	if (o==1) { // return a character vector
		PROTECT(ret = allocVector(STRSXP, x_length));
	} else { // return a list of matrices (x_i.length x 3)
		PROTECT(ret = allocVector(VECSXP, x_length));
	}
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		int *residues = Calloc(x_i.length, int);
		l = 0;
		for (j = 0; j < x_i.length; j++) {
			switch (x_i.ptr[j]) {
				case 65: // A
					residues[l] = 0;
					l++;
					break;
				case 82: // R
					residues[l] = 1;
					l++;
					break;
				case 78: // N
					residues[l] = 2;
					l++;
					break;
				case 68: // D
					residues[l] = 3;
					l++;
					break;
				case 67: // C
					residues[l] = 4;
					l++;
					break;
				case 81: // Q
					residues[l] = 5;
					l++;
					break;
				case 69: // E
					residues[l] = 6;
					l++;
					break;
				case 71: // G
					residues[l] = 7;
					l++;
					break;
				case 72: // H
					residues[l] = 8;
					l++;
					break;
				case 73: // I
					residues[l] = 9;
					l++;
					break;
				case 76: // L
					residues[l] = 10;
					l++;
					break;
				case 75: // K
					residues[l] = 11;
					l++;
					break;
				case 77: // M
					residues[l] = 12;
					l++;
					break;
				case 70: // F
					residues[l] = 13;
					l++;
					break;
				case 80: // P
					residues[l] = 14;
					l++;
					break;
				case 83: // S
					residues[l] = 15;
					l++;
					break;
				case 84: // T
					residues[l] = 16;
					l++;
					break;
				case 87: // W
					residues[l] = 17;
					l++;
					break;
				case 89: // Y
					residues[l] = 18;
					l++;
					break;
				case 86: // V
					residues[l] = 19;
					l++;
					break;
				case 85: // U
					residues[l] = 20; // unknown
					l++;
					break;
				case 79: // O
					residues[l] = 20; // unknown
					l++;
					break;
				case 66: // B = N or D
					residues[l] = 2; // treat as N
					l++;
					break;
				case 90: // Z = Q or E
					residues[l] = 5; // treat as Q
					l++;
					break;
				case 74: // J = I or L
					residues[l] = 9; // treat as I
					l++;
					break;
				case 88: // X = any letter
					residues[l] = 20; // unknown
					l++;
					break;
				case 42: // * (stop)
					residues[l] = 20; // unknown
					l++;
					break;
				case 45: // -
					break; // does not participate
				case 43: // +
					break; // does not participate
				case 46: // . treated as -
					break; // does not participate
				default:
					error("not AA!");
					break;
			}
		}
		
		if (o == 1) {
			states = Calloc(l + 1, char); // last position is for null terminating
		} else  {
			PROTECT(ans = allocMatrix(REALSXP, 3, l)); // [state][pos]
			rans = REAL(ans);
		}
		
		for (j = 0; j < l; j++) {
			H = *(bg);
			E = *(bg + 1);
			C = *(bg + 2);
			
			for (k1 = -1*wS; k1 <= wS; k1++) {
				p1 = j + k1;
				if (p1 < 0)
					continue;
				if (p1 >= l)
					break;
				
				if (residues[p1] < 20) {
					// add mutual information from single residues
					H -= f1 * *(MI1 + residues[p1] + 20*(center + k1));
					E -= f1 * *(MI1 + residues[p1] + 20*(center + k1) + 20*total);
					C -= f1 * *(MI1 + residues[p1] + 20*(center + k1) + 40*total);
				} else {
					continue;
				}
				
				for (k2 = wS; k2 > k1; k2--) {
					p2 = j + k2;
					if (p2 < 0)
						break;
					if (p2 >= l)
						continue;
					
					if (residues[p2] < 20) {
						// add mutual information from pairs of residues
						H += f2 * *(MI2 + residues[p1] + 20*residues[p2] + 400*(center + k1) + 400*total*(center + k2));
						E += f2 * *(MI2 + residues[p1] + 20*residues[p2] + 400*(center + k1) + 400*total*(center + k2) + 400*total*total);
						C += f2 * *(MI2 + residues[p1] + 20*residues[p2] + 400*(center + k1) + 400*total*(center + k2) + 800*total*total);
					}
				}
			}
			
			if (o == 1) { // states
				if (H > E && H > C) {
					states[j] = 'H';
				} else if (E > C) {
					states[j] = 'E';
				} else {
					states[j] = 'C';
				}
			} else if (o == 2) { // log-odds scores
				*(rans + 3*j) = H;
				*(rans + 3*j + 1) = E;
				*(rans + 3*j + 2) = C;
			} else { // normalized probabilities
				H = exp(H)/3;
				E = exp(E)/3;
				C = exp(C)/3;
				sum = H + E + C;
				*(rans + 3*j) = H/sum;
				*(rans + 3*j + 1) = E/sum;
				*(rans + 3*j + 2) = C/sum;
			}
		}
		
		if (o == 1) {
			states[l] = '\0'; // end (null terminate) the string
			SET_STRING_ELT(ret, i, mkChar(states));
			Free(states);
		} else {
			SET_VECTOR_ELT(ret, i, ans);
			UNPROTECT(1); // ans
		}
		
		Free(residues);
	}
	
	UNPROTECT(1);
	
	return ret;
}
