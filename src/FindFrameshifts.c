/****************************************************************************
 *                      Finds Protein Reading Frames                        *
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

static void assignNumsAA(const Chars_holder *P, int *residues)
{
	int j;
	const char *p;
	
	for (j = 0, p = (P->ptr);
		 j < P->length;
		 j++, p++)
	{
		switch (*p) {
			case 65: // A
				*(residues + j) = 0;
				break;
			case 82: // R
				*(residues + j) = 1;
				break;
			case 78: // N
				*(residues + j) = 2;
				break;
			case 68: // D
				*(residues + j) = 3;
				break;
			case 67: // C
				*(residues + j) = 4;
				break;
			case 81: // Q
				*(residues + j) = 5;
				break;
			case 69: // E
				*(residues + j) = 6;
				break;
			case 71: // G
				*(residues + j) = 7;
				break;
			case 72: // H
				*(residues + j) = 8;
				break;
			case 73: // I
				*(residues + j) = 9;
				break;
			case 76: // L
				*(residues + j) = 10;
				break;
			case 75: // K
				*(residues + j) = 11;
				break;
			case 77: // M
				*(residues + j) = 12;
				break;
			case 70: // F
				*(residues + j) = 13;
				break;
			case 80: // P
				*(residues + j) = 14;
				break;
			case 83: // S
				*(residues + j) = 15;
				break;
			case 84: // T
				*(residues + j) = 16;
				break;
			case 87: // W
				*(residues + j) = 17;
				break;
			case 89: // Y
				*(residues + j) = 18;
				break;
			case 86: // V
				*(residues + j) = 19;
				break;
			case 85: // U
				*(residues + j) = 20; // unknown
				break;
			case 79: // O
				*(residues + j) = 20; // unknown
				break;
			case 66: // B = N or D
				*(residues + j) = 2; // treat as N
				break;
			case 90: // Z = Q or E
				*(residues + j) = 5; // treat as Q
				break;
			case 74: // J = I or L
				*(residues + j) = 9; // treat as I
				break;
			case 88: // X = any letter
				*(residues + j) = 20; // unknown
				break;
			case 42: // * (stop)
				*(residues + j) = 20; // unknown
				break;
			default:
				error("not AA!");
				break;
		}
	}
}

SEXP findFrameshifts(SEXP t, SEXP l, SEXP f, SEXP index, SEXP maxComp, SEXP go, SEXP ge, SEXP fs, SEXP minD, SEXP maxD, SEXP subMatrix, SEXP verbose, SEXP pBar)
{
	int s, o, i, j, k, I, J, K, n, m, w, r, c, rc;
	int nIns, nDels, pm, d, pos, newL;
	double m0, m1, m2, m3, match, gapCol, gapRow;
	double bestD, tempD;
	double GO = asReal(go);
	double GE = asReal(ge);
	double FS = asReal(fs);
	int mC = asInteger(maxComp);
	double aD = asReal(minD);
	double rD = asReal(maxD);
	int *ind = INTEGER(index);
	double *subM = REAL(subMatrix);
	int *lengths = INTEGER(l);
	
	XStringSet_holder t_set, f_set;
	t_set = hold_XStringSet(t);
	f_set = hold_XStringSet(f);
	int f_length = get_length_from_XStringSet_holder(&f_set);
	f_length /= 3;
	Chars_holder x_s;
	
	SEXP distance, refNum, insertions, deletions, list_comp, ret_list;
	PROTECT(ret_list = allocVector(VECSXP, f_length));
	
	int before, v, *rPercentComplete;
	v = asLogical(verbose);
	SEXP percentComplete, utilsPackage;
	if (v) { // percent complete variables
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	for (s = 0; s < f_length; s++) { // sequence
		x_s = get_elt_from_XStringSet_holder(&f_set, s); // 1st reading frame
		m = x_s.length;
		int *f1 = Calloc(m, int); // initialized to zero
		assignNumsAA(&x_s, f1);
		x_s = get_elt_from_XStringSet_holder(&f_set, s + f_length); // 2nd reading frame
		if (x_s.length < m)
			m = x_s.length;
		int *f2 = Calloc(x_s.length, int); // initialized to zero
		assignNumsAA(&x_s, f2);
		x_s = get_elt_from_XStringSet_holder(&f_set, s + 2*f_length); // 3rd reading frame
		if (x_s.length < m)
			m = x_s.length;
		int *f3 = Calloc(x_s.length, int); // initialized to zero
		assignNumsAA(&x_s, f3);
		
		bestD = 1000; // initialize to large distance
		
		for (o = 0; o < mC; o++) { // order of reference translations
			x_s = get_elt_from_XStringSet_holder(&t_set, *(ind + o*f_length + s) - 1);
			n = x_s.length;
			int *ref = Calloc(n, int); // initialized to zero
			assignNumsAA(&x_s, ref);
			
			r = m + 1; // number of rows
			c = n + 1; // number of columns
			rc = r*c; // length of each matrix
			
			// initialize an array of scores
			double *A = Calloc(rc*3, double); // initialized to zero
			// initialize an array of arrows
			int *B = Calloc(rc*3, int); // initialized to zero
			// initialize an array of shifts
			int *C = Calloc(rc*3, int); // initialized to zero
			memset(C, -1, rc*3 * sizeof(int)); // initialize to -1
			
			// indexing:  A[k*rc + j*r + i]==A[i, j, k]
			
			for (i = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					// find the best scoring frame
					m1 = A[j*r + i];
					m2 = A[rc + j*r + i];
					m3 = A[2*rc + j*r + i];
					if (m1 > m2 && m1 > m3) {
						w = 0;
					} else if (m2 > m3) {
						w = 1;
					} else {
						w = 2;
					}
					
					for (k = 0; k < 3; k++) {
						// look for frameshifts from previous i and j
						if (k != w &&
							(A[w*rc + j*r + i] + FS) > A[k*rc + j*r + i]) {
							A[k*rc + j*r + i] = A[w*rc + j*r + i] + FS;
							B[k*rc + j*r + i] = 0;
							C[k*rc + j*r + i] = w;
						}
						
						// score matching
						match = A[k*rc + j*r + i];
						if (k==0) {
							match += *(subM + 21*ref[j] + f1[i]);
						} else if (k==1) {
							match += *(subM + 21*ref[j] + f2[i]);
						} else {
							match += *(subM + 21*ref[j] + f3[i]);
						}
						
						// score gap in ref
						if (B[k*rc + j*r + i + 1] > 0) {
							gapCol = GE;
						} else {
							gapCol = GO;
						}
						gapCol += A[k*rc + (j + 1)*r + i];
						
						// score gap in query
						if (B[k*rc + (j + 1)*r + i] < 0) {
							gapRow = GE;
						} else {
							gapRow = GO;
						}
						gapRow += A[k*rc + j*r + i + 1];
						
						if (match >= gapRow &&
							match >= gapCol) {
							A[k*rc + (j + 1)*r + i + 1] = match;
							B[k*rc + (j + 1)*r + i + 1] = 0;
						} else if (gapRow > gapCol) {
							A[k*rc + (j + 1)*r + i + 1] = gapRow;
							if (B[k*rc + j*r + i + 1] > 0) {
								B[k*rc + (j + 1)*r + i + 1] = B[k*rc + j*r + i + 1] + 1;
							} else {
								B[k*rc + (j + 1)*r + i + 1] = 1;
							}
						} else {
							A[k*rc + (j + 1)*r + i + 1] = gapCol;
							if (B[k*rc + (j + 1)*r + i] < 0) {
								B[k*rc + (j + 1)*r + i + 1] = B[k*rc + (j + 1)*r + i] - 1;
							} else {
								B[k*rc + (j + 1)*r + i + 1] = -1;
							}
						}
					}
				}
			}
			
			// find the max score
			m0 = -1e53; // max value
			j = n;
			for (i = 0; i <= m; i++) {
				for (k = 0; k < 3; k++) {
					if (A[k*rc + j*r + i] > m0) {
						I = i;
						J = j;
						K = k;
						m0 = A[k*rc + j*r + i];
					}
				}
			}
			i = m;
			for (j = 0; j < n; j++) { // corner already visited
				for (k = 0; k < 3; k++) {
					if (A[k*rc + j*r + i] > m0) {
						I = i;
						J = j;
						K = k;
						m0 = A[k*rc + j*r + i];
					}
				}
			}
			
			i = I;
			j = J;
			k = K;
			
			// traceback
			int *ins = Calloc(1000, int); // initialized to zero
			int *dels = Calloc(1000, int); // initialized to zero
			nIns = 0;
			nDels = 0;
			pm = 0;
			
			while (i > 0 && j > 0) {
				if (C[k*rc + j*r + i] >= 0) {
					d = k - C[k*rc + j*r + i];
					pos = i*3 + k + 1;
					if (d==1 || d==-2) {
						if (B[C[k*rc + j*r + i]*rc + j*r + i] > 0) {
							// missing nucleotides
							if (nDels > 998) {
								error("Out of memory");
							}
							dels[nDels++] = pos;
							dels[nDels++] = pos;
						} else {
							// extra nucleotide
							if (nIns > 999) {
								error("Out of memory");
							}
							ins[nIns++] = pos;
						}
					} else { // d==-1 || d==2
						if (B[C[k*rc + j*r + i]*rc + j*r + i] >= 0) {
							// missing nucleotide
							if (nDels > 999) {
								error("Out of memory");
							}
							dels[nDels++] = pos;
						} else {
							// extra nucleotides
							if (nIns > 998) {
								error("Out of memory");
							}
							ins[nIns++] = pos;
							ins[nIns++] = pos;
						}
					}
					k = C[k*rc + j*r + i];
				}
				
				if (B[k*rc + j*r + i]==0) {
					i--;
					j--;
					
					if (k==0) {
						if (ref[j]==f1[i])
							pm++;
					} else if (k==1) {
						if (ref[j]==f2[i])
							pm++;
					} else { // k==2
						if (ref[j]==f3[i])
							pm++;
					}
				} else if (B[k*rc + j*r + i] < 0) {
					i += B[k*rc + j*r + i];
				} else {
					j -= B[k*rc + j*r + i];
				}
			}
			
			if (k==1) {
				if (nIns > 999) {
					error("Out of memory");
				}
				ins[nIns++] = 1;
			} else if (k==2) {
				if (nIns > 998) {
					error("Out of memory");
				}
				ins[nIns++] = 1;
				ins[nIns++] = 1;
			}
			
			tempD = 1 - (double)pm/m;
			if (tempD >= bestD) { // no improvement
				Free(ref);
				Free(A);
				Free(B);
				Free(C);
				Free(ins);
				Free(dels);
				
				continue;
			}
			bestD = tempD;
			
			PROTECT(distance = allocVector(REALSXP, 1));
			double *dist = REAL(distance);
			*(dist) = tempD;
			
			PROTECT(refNum = allocVector(INTSXP, 1));
			
			if (*(dist) > rD) {
				nIns = 0;
				nDels = 0;
				INTEGER(refNum)[0] = 0;
			} else {
				INTEGER(refNum)[0] = *(ind + o*f_length + s);
			}
			
			newL = *(lengths + s) - nIns + nDels;
			if ((newL % 3)==1) {
				if (nIns > 999) {
					error("Out of memory");
				}
				ins[nIns++] = *(lengths + s);
			} else if ((newL % 3)==2) {
				if (nIns > 998) {
					error("Out of memory");
				}
				pos = *(lengths + s) - 1;
				ins[nIns++] = pos;
				ins[nIns++] = pos;
			}
			
			PROTECT(insertions = allocVector(INTSXP, nIns));
			int *indels1 = INTEGER(insertions);
			for (i = 0; i < nIns; i++)
				*(indels1 + i) = ins[i];
			PROTECT(deletions = allocVector(INTSXP, nDels));
			int *indels2 = INTEGER(deletions);
			for (i = 0; i < nDels; i++)
				*(indels2 + i) = dels[i];
			
			PROTECT(list_comp = allocVector(VECSXP, 4));
			SET_VECTOR_ELT(list_comp, 0, insertions);
			SET_VECTOR_ELT(list_comp, 1, deletions);
			SET_VECTOR_ELT(list_comp, 2, distance);
			SET_VECTOR_ELT(list_comp, 3, refNum);
			SET_VECTOR_ELT(ret_list, s, list_comp);
			UNPROTECT(5);
			
			Free(ref);
			Free(A);
			Free(B);
			Free(C);
			Free(ins);
			Free(dels);
			
			if (*(dist) <= aD)
				break;
		}
		
		Free(f1);
		Free(f2);
		Free(f3);
		
		if (v) {
			*rPercentComplete = floor(100*(double)(s + 1)/f_length);
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ret_list;
}
