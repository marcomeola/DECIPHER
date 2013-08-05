/****************************************************************************
 *                       Aligns Two Sequence Profiles                       *
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

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// DECIPHER header file
#include "DECIPHER.h"

// for div_t div
#include <stdlib.h>

//ans_start <- .Call("alignProfiles", p.profile, s.profile, pm, mm, go, ge, endGapPenalty, PACKAGE="DECIPHER")//, codon, subMatrix
SEXP alignProfiles(SEXP p, SEXP s, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP nThreads)//, SEXP codonBonus, SEXP subMatrix
{
	int i, j, k, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, S, M, GP, GS;
	SEXP ans;
	
	pprofile = REAL(p);
	sprofile = REAL(s);
	double PM = asReal(pm);
	double MM = asReal(mm);
	double GO = asReal(go);
	double GE = asReal(ge);
	double egpL = asReal(endGapPenaltyLeft);
	double egpR = asReal(endGapPenaltyRight);
	//double codon = asReal(codonBonus);
	//double *subM = REAL(subMatrix);
	int nthreads = asInteger(nThreads);
	
	R_len_t lp = length(p)/7;
	R_len_t ls = length(s)/7;
	
	PROTECT(ans = allocVector(INTSXP, lp + ls + 4));
	rans = INTEGER(ans);
	
	float *m = Calloc((lp+1)*(ls+1), float); // initialized to zero
	int *o = Calloc(lp*ls, int); // initialized to zero
	
	if (egpL != 0) {
		*(m + 1) = egpL*(1 - sprofile[4]);
		for (i = 2; i <= ls; i++)
			*(m + i) = *(m + i - 1) + egpL*(1 - sprofile[4 + 7*(i - 1)]);
		*(m + ls + 1) = egpL*(1 - pprofile[4]);
		for (i = 2; i <= lp; i++)
			*(m + (ls + 1)*i) = *(m + (ls + 1)*(i - 1)) + egpL*(1 - pprofile[4 + 7*(i - 1)]);
	}
	
	for (k = 0; k < lp + ls; k++) {
		if (k > ls) {
			start = k - ls;
		} else {
			start = 0;
		}
		if (k >= lp) {
			end = lp - 1;
		} else {
			end = k - 1;
		}
		
		#pragma omp parallel for private(i,j,gp,gs,S,M,GP,GS) schedule(guided) num_threads(nthreads)
		for (i = start; i <= end; i++) {
			// determine column index
			if (k >= lp) {
				if (k >= ls) {
					j = ls - i + start - 1;
				} else {
					j = end - i + 1 - lp + k - 1;
				}
			} else {
				j = end - i;
			}
			
			// determine insertion penalty (gap extension or gap opening)
			gp = 2*GO*(1 - sprofile[5 + 7*j]);// + GE*(1 - sprofile[4 + 7*j]);
			gs = 2*GO*(1 - pprofile[5 + 7*i]);// + GE*(1 - pprofile[4 + 7*i]);
			if (j > 0) {
				if (*(o + i*ls + j - 1) > 0) {
					if (*(o + i*ls + j - 1) == 1) {
						gp = GE*(1 - sprofile[4 + 7*j]) - GO*(1 - sprofile[5 + 7*(j - 1)]) + GO*(1 - sprofile[6 + 7*j]);
					} else {
						gp = GE*(1 - sprofile[4 + 7*j]) + GO*(1 - sprofile[6 + 7*j]) - GO*(1 - sprofile[6 + 7*(j - 1)]);//*sprofile[4 + 7*(j - 1)] + GO*(1 - sprofile[4 + 7*(j - 1)]);
					}
				}
			}
			if (i > 0) {
				if (*(o + (i - 1)*ls + j) < 0) {
					if (*(o + (i - 1)*ls + j) == -1) {
						gs = GE*(1 - pprofile[4 + 7*i]) - GO*(1 - pprofile[5 + 7*(i - 1)]) + GO*(1 - pprofile[6 + 7*i]);
					} else {
						gs = GE*(1 - pprofile[4 + 7*i]) + GO*(1 - pprofile[6 + 7*i]) - GO*(1 - pprofile[6 + 7*(i - 1)]);//*pprofile[4 + 7*(i - 1)] + GO*(1 - pprofile[4 + 7*(i - 1)]);
					}
				}
			}
			// score aligning
			S = sqrt(pprofile[7*i]*sprofile[7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[4 + 7*j]);
			//M = *(m + i*(ls + 1) + j) + (sqrt(pprofile[0 + 7*i]*sprofile[1 + 7*j]) * *(subM + 1*4 + 0) + sqrt(pprofile[0 + 7*i]*sprofile[2 + 7*j]) * *(subM + 2*4 + 0) + sqrt(pprofile[0 + 7*i]*sprofile[3 + 7*j]) * *(subM + 3*4 + 0) + sqrt(pprofile[1 + 7*i]*sprofile[0 + 7*j]) * *(subM + 0*4 + 1) + sqrt(pprofile[1 + 7*i]*sprofile[2 + 7*j]) * *(subM + 2*4 + 1) + sqrt(pprofile[1 + 7*i]*sprofile[3 + 7*j]) * *(subM + 3*4 + 1) + sqrt(pprofile[2 + 7*i]*sprofile[0 + 7*j]) * *(subM + 0*4 + 2) + sqrt(pprofile[2 + 7*i]*sprofile[1 + 7*j]) * *(subM + 1*4 + 2) + sqrt(pprofile[2 + 7*i]*sprofile[3 + 7*j]) * *(subM + 3*4 + 2) + sqrt(pprofile[3 + 7*i]*sprofile[0 + 7*j]) * *(subM + 0*4 + 3) + sqrt(pprofile[3 + 7*i]*sprofile[1 + 7*j]) * *(subM + 1*4 + 3) + sqrt(pprofile[3 + 7*i]*sprofile[2 + 7*j]) * *(subM + 2*4 + 3) + sqrt(pprofile[4 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[4 + 7*j]))*MM + S*PM;
			M = *(m + i*(ls + 1) + j) + (sqrt(pprofile[0 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[4 + 7*j]))*MM + S*PM;
			GP = *(m + (i + 1)*(ls + 1) + j) + gp + PM*sprofile[4 + 7*j] - PM*pprofile[4 + 7*i];// + MM*(1 - pprofile[4 + 7*i]);
			GS = *(m + i*(ls + 1) + j + 1) + gs + PM*pprofile[4 + 7*i] - PM*sprofile[4 + 7*j];// + MM*(1 - sprofile[4 + 7*j]);
			
			/*
			if (codon != 0 && i > 0) {
				div_t s_div = div(*(o + (i - 1)*ls + j), 3);
				if (s_div.rem == -2) {
					GS -= (s_div.quot - 1)*codon;
				} else if (s_div.quot < 0 && s_div.rem == 0) {
					GS += s_div.quot*codon;
				}
			}
			
			if (codon != 0 && j > 0) {
				div_t p_div = div(*(o + i*ls + j - 1), 3);
				if (p_div.rem == 2) {
					GP += (p_div.quot + 1)*codon;
				} else if (p_div.quot > 0 && p_div.rem == 0) {
					GP -= p_div.quot*codon;
				}
			}
			*/
			
			if ((GS >= M) && (GS >= GP)) {
				if (i > 0) {
					if (*(o + (i - 1)*ls + j) < 0) {
						*(o + i*ls + j) = *(o + (i - 1)*ls + j) - 1;
					} else {
						*(o + i*ls + j) = -1;
					}
				} else {
					*(o + i*ls + j) = -1;
				}
				
				*(m + (i + 1)*(ls + 1) + j + 1) = GS;
			} else if (GP >= M) {
				if (j > 0) {
					if (*(o + i*ls + j - 1) > 0) {
						*(o + i*ls + j) = *(o + i*ls + j - 1) + 1;
					} else {
						*(o + i*ls + j) = 1;
					}
				} else {
					*(o + i*ls + j) = 1;
				}
				
				*(m + (i + 1)*(ls + 1) + j + 1) = GP;
			} else {
				*(m + (i + 1)*(ls + 1) + j + 1) = M;
				*(o + i*ls + j) = 0;
			}
		}
	}
	
	if (egpR != 0) {
		GP = egpR*(1 - sprofile[4 + 7*(ls - 1)]);
		*(m + lp*(ls + 1) + ls - 1) = *(m + lp*(ls + 1) + ls - 1) + GP;
		for (i = ls - 2; i > -1; i--) {
			GP += egpR*(1 - sprofile[4 + 7*i]);
			*(m + lp*(ls + 1) + i) = *(m + lp*(ls + 1) + i) + GP;
		}
		GP = egpR*(1 - pprofile[4 + 7*(lp - 1)]);
		*(m + (lp - 1)*(ls + 1) + ls) = *(m + (lp - 1)*(ls + 1) + ls) + GP;
		for (i = lp - 2; i > -1; i--) {
			GP += egpR*(1 - pprofile[4 + 7*i]);
			*(m + i*(ls + 1) + ls) = *(m + i*(ls + 1) + ls) + GP;
		}
	}
	
	// find the max scoring alignment
	int maxp = 0;
	int maxs = 0;
	for (i = 1; i <= ls; i++)
		if (*(m + lp*(ls + 1) + i) >= *(m + lp*(ls + 1) + maxs))
			maxs = i;
	for (i = 1; i <= lp; i++)
		if (*(m + i*(ls + 1) + ls) >= *(m + maxp*(ls + 1) + ls))
			maxp = i;
	
	if (*(m + lp*(ls + 1) + maxs) > *(m + maxp*(ls + 1) + ls)) {
		*(rans) = lp;
		*(rans + 1) = maxs;
	} else {
		*(rans) = maxp;
		*(rans + 1) = ls;
	}
	
	i = (int)*(rans) - 1;
	j = (int)*(rans + 1) - 1;
	count = lp + ls + 3;
	while ((i > -1) && (j > -1)) {
		*(rans + count) = *(o + i*ls + j);
		z = *(o + i*ls + j);
		//*(o + i*ls + j) = 0;
		if (z <= 0)
			i--;
		if (z >= 0)
			j--;
		count--;
	}
	*(rans + 2) = i + 2;
	*(rans + 3) = j + 2;
	*(rans + 4) = count;
	
	/*
	for (i = 0; i < (lp+1)*(ls+1); i++) {
		if (i % (ls + 1) == 0) {
			Rprintf("\n");
		}
		Rprintf("%d ", (int)m[i]);
	}
	Rprintf("\n");
	for (i = 0; i < lp*ls; i++) {
		if (i % ls == 0) {
			Rprintf("\n");
		}
		Rprintf("%d ", (int)o[i]);
	}
	Rprintf("\n");
	*/
	
	Free(m);
	Free(o);
	
	UNPROTECT(1);
	
	return ans; // vector = [end.p end.s start.p start.s count t]
}
