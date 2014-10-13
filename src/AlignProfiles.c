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
//#include <stdlib.h>

//ans_start <- .Call("alignProfiles", p.profile, s.profile, subMatrix, pm, mm, go, ge, endGapPenalty, boundary, processors, PACKAGE="DECIPHER")//, codon
SEXP alignProfiles(SEXP p, SEXP s, SEXP subMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)//, SEXP codonBonus
{
	int i, j, k, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, S, M, GP, GS;
	double max, tot, PMnorm, MMnorm, GOnorm, GEnorm;
	SEXP ans;
	
	pprofile = REAL(p);
	sprofile = REAL(s);
	double PM = asReal(pm);
	double MM = asReal(mm);
	double GO = asReal(go);
	double GE = asReal(ge);
	double bound = asReal(boundary);
	double egpL = asReal(endGapPenaltyLeft);
	double egpR = asReal(endGapPenaltyRight);
	//double codon = asReal(codonBonus);
	double *subM;
	int do_subM;
	if (length(subMatrix)==0) {
		do_subM = 0;
	} else {
		do_subM = 1;
		subM = REAL(subMatrix);
	}
	int NTHREADS, nthreads = asInteger(nThreads);
	
	R_len_t lp = length(p)/8;
	R_len_t ls = length(s)/8;
	
	PROTECT(ans = allocVector(INTSXP, lp + ls + 4));
	rans = INTEGER(ans);
	
	float *m = Calloc((lp+1)*(ls+1), float); // initialized to zero
	int *o = Calloc(lp*ls, int); // initialized to zero
	
	// normalize profile positions to sqrt(coverge)
	for (i = 0; i < lp; i++)
		pprofile[7 + 8*i] = sqrt(pprofile[7 + 8*i]);
	for (i = 0; i < ls; i++)
		sprofile[7 + 8*i] = sqrt(sprofile[7 + 8*i]);
	
	if (egpL != 0) {
		*(m + 1) = egpL*(1 - sprofile[4])*sprofile[7];
		for (i = 2; i <= ls; i++)
			*(m + i) = *(m + i - 1) + egpL*(1 - sprofile[4 + 8*(i - 1)])*sprofile[7 + 8*(i - 1)];
		*(m + ls + 1) = egpL*(1 - pprofile[4])*pprofile[7];
		for (i = 2; i <= lp; i++)
			*(m + (ls + 1)*i) = *(m + (ls + 1)*(i - 1)) + egpL*(1 - pprofile[4 + 8*(i - 1)])*pprofile[7 + 8*(i - 1)];
	}
	
	int left = 0, top = 0; // boundaries of DP
	int START, END; // constrained start and end points
	max = -1e53;
	for (k = 1; k < lp + ls; k++) {
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
		
		max += bound; // threshold for constraint boundary
		
		if (top > start) {
			START = top;
		} else {
			START = start;
		}
		for (i = START + 1; i <= end; i++) {
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
			//Rprintf("\ni%d j%d top%d START%d start%d value%d", i, j, top, START, start, (int)*(m + i*(ls + 1) + j));
			if (*(m + i*(ls + 1) + j) < max) {
				if (j < ls && i > 1) {
					//Rprintf("\n!\ntop %d i %d j %d ls %d lp %d\n!\n", top, i, j, ls, lp);
					*(o + i*ls + j) = 1;
					*(m + i*(ls + 1) + ls) = -1e53; // block restricted area from traceback
					top = i;
					START = i;
				}
			} else {
				break;
			}
		}
		//Rprintf("\nk %d START %d END %d start %d end %d top %d left %d", k, START, END, start, end, top, left);
		END = end;
		for (i = END; i > START; i--) {
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
			//Rprintf("\ni%d j%d left%d END%d value%d", i, j, left, END, (int)*(m + i*(ls + 1) + j));
			if (j < left) {
				i -= left - j;
				continue;
			} else if (*(m + i*(ls + 1) + j) < max && i != (lp - 1)) {
				left = j;
				*(o + i*ls + j) = -1;
				*(m + lp*(ls + 1) + left + 1) = -1e53; // block restricted area from traceback
			} else {
				if (i < END) {
					END = i + 1;
				} else {
					END = i;
				}
				left = j - 1;
				break;
			}
		}
		if (END < end)
			*(m + (END + 1)*(ls + 1) + left) = max + 2*GO; // set edge boundary
		if (START > start) {
			if (k >= lp) {
				if (k >= ls) {
					j = ls - START + start - 1;
				} else {
					j = end - START + 1 - lp + k - 1;
				}
			} else {
				j = end - START;
			}
			*(m + START*(ls + 1) + j + 1) = max + 2*GO; // set edge boundary
		}
		//Rprintf("\nk %d START %d END %d start %d end %d top %d left %d", k, START, END, start, end, top, left);
		max = -1e53;
		if (END - START > 2000) {
			NTHREADS = nthreads;
		} else {
			NTHREADS = 1;
		}
		#pragma omp parallel for private(i,j,gp,gs,S,M,GP,GS) num_threads(NTHREADS)
		for (i = START; i <= END; i++) {
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
			
			tot = sqrt(sprofile[7 + 8*j]*pprofile[7 + 8*i]);
			PMnorm = PM*tot;
			MMnorm = MM*tot;
			GOnorm = GO*tot;
			GEnorm = GE*tot;
			
			// determine insertion penalty (gap extension or gap opening)
			gp = 2*GOnorm*(1 - sprofile[5 + 8*j]);
			gs = 2*GOnorm*(1 - pprofile[5 + 8*i]);
			
			if (j > 0) {
				if (*(o + i*ls + j - 1) > 0) {
					if (*(o + i*ls + j - 1) == 1) {
						gp = GEnorm*(1 - sprofile[4 + 8*j]) - GOnorm*(1 - sprofile[5 + 8*(j - 1)]) + GOnorm*(1 - sprofile[6 + 8*j]);
					} else {
						gp = GEnorm*(1 - sprofile[4 + 8*j]) + GOnorm*(1 - sprofile[6 + 8*j]) - GOnorm*(1 - sprofile[6 + 8*(j - 1)]);
					}
				}
			}
			if (i > 0) {
				if (*(o + (i - 1)*ls + j) < 0) {
					if (*(o + (i - 1)*ls + j) == -1) {
						gs = GEnorm*(1 - pprofile[4 + 8*i]) - GOnorm*(1 - pprofile[5 + 8*(i - 1)]) + GOnorm*(1 - pprofile[6 + 8*i]);
					} else {
						gs = GEnorm*(1 - pprofile[4 + 8*i]) + GOnorm*(1 - pprofile[6 + 8*i]) - GOnorm*(1 - pprofile[6 + 8*(i - 1)]);
					}
				}
			}
			
			// apply end-gap penalties at corners of matrix
			if (i==0 && j==0) {
				gp = egpL*(1 - sprofile[5 + 8*j])*sprofile[7 + 8*j];
				gs = egpL*(1 - pprofile[5 + 8*i])*pprofile[7 + 8*i];
			} else if (i==(lp - 1) && j==(ls - 1)) {
				gp = egpR*(1 - sprofile[5 + 8*j])*sprofile[7 + 8*j];
				gs = egpR*(1 - pprofile[5 + 8*i])*pprofile[7 + 8*i];
			}
			
			// score aligning
			if (do_subM) {
				S = 0;
				M = 0;
				for (int n = 0; n < 5; n++) {
					if (pprofile[n + 8*i]==0)
						continue;
					if (n < 4) {
						S += pprofile[n + 8*i]*sprofile[n + 8*j] * *(subM + n*4 + n)*tot;
					} else {
						S += pprofile[n + 8*i]*sprofile[n + 8*j]*PMnorm;
					}
					for (int p = 0; p < 5; p++) {
						if (p==n || sprofile[p + 8*j]==0)
							continue;
						if (p < 4 && n < 4) {
							M += pprofile[n + 8*i]*sprofile[p + 8*j] * *(subM + n*4 + p)*tot;
						} else {
							M += pprofile[n + 8*i]*sprofile[p + 8*j]*MMnorm;
						}
					}
				}
				M = *(m + i*(ls + 1) + j) + S + M;
			} else {
				S = 0;
				M = 0;
				for (int n = 0; n < 5; n++) {
					if (pprofile[n + 8*i]==0)
						continue;
					S += pprofile[n + 8*i]*sprofile[n + 8*j];
				}
				M = *(m + i*(ls + 1) + j) + S*PMnorm + (1 - S)*MMnorm; // % mismatched = 100% - % matched
			}
			GP = *(m + (i + 1)*(ls + 1) + j) + gp + PMnorm*sprofile[4 + 8*j] - PMnorm*pprofile[4 + 8*i] + MMnorm*(1 - sprofile[4 + 8*j]);
			GS = *(m + i*(ls + 1) + j + 1) + gs + PMnorm*pprofile[4 + 8*i] - PMnorm*sprofile[4 + 8*j] + MMnorm*(1 - pprofile[4 + 8*i]);
			
			if ((GS >= M) && (GS >= GP) && sprofile[4 + 8*j] != 1) {
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
			} else if (GP >= M && pprofile[4 + 8*i] != 1) {
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
			
			//#pragma omp critical // not necessary because will only underestimate max
			if (*(m + (i + 1)*(ls + 1) + j + 1) > max)
				max = *(m + (i + 1)*(ls + 1) + j + 1);
		}
	}
	
	if (egpR != 0) {
		GP = egpR*(1 - sprofile[4 + 8*(ls - 1)])*sprofile[7 + 8*(ls - 1)];
		*(m + lp*(ls + 1) + ls - 1) = *(m + lp*(ls + 1) + ls - 1) + GP;
		for (i = ls - 2; i > -1; i--) {
			GP += egpR*(1 - sprofile[4 + 8*i])*sprofile[7 + 8*i];
			*(m + lp*(ls + 1) + i) = *(m + lp*(ls + 1) + i) + GP;
		}
		GP = egpR*(1 - pprofile[4 + 8*(lp - 1)])*pprofile[7 + 8*(lp - 1)];
		*(m + (lp - 1)*(ls + 1) + ls) = *(m + (lp - 1)*(ls + 1) + ls) + GP;
		for (i = lp - 2; i > -1; i--) {
			GP += egpR*(1 - pprofile[4 + 8*i])*pprofile[7 + 8*i];
			*(m + i*(ls + 1) + ls) = *(m + i*(ls + 1) + ls) + GP;
		}
	}
	
	// find the max scoring alignment
	int maxp = 0;
	int maxs = 0;
	for (i = 1; i <= ls; i++)
		if (*(m + lp*(ls + 1) + i) >= *(m + lp*(ls + 1) + maxs))
			maxs = i;
	for (i = top + 1; i <= lp; i++)
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

//ans_start <- .Call("alignProfilesAA", p.profile, s.profile, subMatrix, pm, mm, go, ge, endGapPenalty, boundary, processors, PACKAGE="DECIPHER")
SEXP alignProfilesAA(SEXP p, SEXP s, SEXP subMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)
{
	int i, j, k, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, M, GP, GS;
	double max, tot, PMnorm, MMnorm, GOnorm, GEnorm;
	SEXP ans;
	
	pprofile = REAL(p);
	sprofile = REAL(s);
	double PM = asReal(pm);
	double MM = asReal(mm);
	double GO = asReal(go);
	double GE = asReal(ge);
	double bound = asReal(boundary);
	double egpL = asReal(endGapPenaltyLeft);
	double egpR = asReal(endGapPenaltyRight);
	double *subM = REAL(subMatrix);
	int NTHREADS, nthreads = asInteger(nThreads);
	
	R_len_t lp = length(p)/27;
	R_len_t ls = length(s)/27;
	
	PROTECT(ans = allocVector(INTSXP, lp + ls + 4));
	rans = INTEGER(ans);
	
	float *m = Calloc((lp+1)*(ls+1), float); // initialized to zero
	int *o = Calloc(lp*ls, int); // initialized to zero
	
	// normalize profile positions to sqrt(coverge)
	for (i = 0; i < lp; i++)
		pprofile[26 + 27*i] = sqrt(pprofile[26 + 27*i]);
	for (i = 0; i < ls; i++)
		sprofile[26 + 27*i] = sqrt(sprofile[26 + 27*i]);
	
	if (egpL != 0) {
		*(m + 1) = egpL*(1 - sprofile[23])*sprofile[26];
		for (i = 2; i <= ls; i++)
			*(m + i) = *(m + i - 1) + egpL*(1 - sprofile[23 + 27*(i - 1)])*sprofile[23 + 26*(i - 1)];
		*(m + ls + 1) = egpL*(1 - pprofile[23])*pprofile[26];
		for (i = 2; i <= lp; i++)
			*(m + (ls + 1)*i) = *(m + (ls + 1)*(i - 1)) + egpL*(1 - pprofile[23 + 27*(i - 1)])*pprofile[26 + 27*(i - 1)];
	}
	
	int left = 0, top = 0; // boundaries of DP
	int START, END; // constrained start and end points
	max = -1e53;
	for (k = 1; k < lp + ls; k++) {
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
		
		max += bound; // threshold for constraint boundary
		
		if (top > start) {
			START = top;
		} else {
			START = start;
		}
		for (i = START + 1; i <= end; i++) {
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
			if (*(m + i*(ls + 1) + j) < max) {
				if (j < ls && i > 1) {
					*(o + i*ls + j) = 1;
					*(m + i*(ls + 1) + ls) = -1e53; // block restricted area from traceback
					top = i;
					START = i;
				}
			} else {
				break;
			}
		}
		
		END = end;
		for (i = END; i > START; i--) {
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
			if (j < left) {
				i -= left - j;
				continue;
			} else if (*(m + i*(ls + 1) + j) < max && i != (lp - 1)) {
				left = j;
				*(o + i*ls + j) = -1;
				*(m + lp*(ls + 1) + left + 1) = -1e53; // block restricted area from traceback
			} else {
				if (i < END) {
					END = i + 1;
				} else {
					END = i;
				}
				left = j - 1;
				break;
			}
		}
		if (END < end)
			*(m + (END + 1)*(ls + 1) + left) = max + 2*GO; // set edge boundary
		if (START > start) {
			if (k >= lp) {
				if (k >= ls) {
					j = ls - START + start - 1;
				} else {
					j = end - START + 1 - lp + k - 1;
				}
			} else {
				j = end - START;
			}
			*(m + START*(ls + 1) + j + 1) = max + 2*GO; // set edge boundary
		}
		
		max = -1e53;
		if (END - START > 2000) {
			NTHREADS = nthreads;
		} else {
			NTHREADS = 1;
		}
		#pragma omp parallel for private(i,j,gp,gs,M,GP,GS) num_threads(NTHREADS)
		for (i = START; i <= END; i++) {
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
			
			tot = sqrt(sprofile[26 + 27*j]*pprofile[26 + 27*i]); // normalization factor
			PMnorm = PM*tot;
			MMnorm = MM*tot;
			GOnorm = GO*tot;
			GEnorm = GE*tot;
			
			// determine insertion penalty (gap extension or gap opening)
			gp = 2*GOnorm*(1 - sprofile[24 + 27*j]);
			gs = 2*GOnorm*(1 - pprofile[24 + 27*i]);
			if (j > 0) {
				if (*(o + i*ls + j - 1) > 0) {
					if (*(o + i*ls + j - 1) == 1) {
						gp = GEnorm*(1 - sprofile[23 + 27*j]) - GOnorm*(1 - sprofile[24 + 27*(j - 1)]) + GOnorm*(1 - sprofile[25 + 27*j]);
					} else {
						gp = GEnorm*(1 - sprofile[23 + 27*j]) + GOnorm*(1 - sprofile[25 + 27*j]) - GOnorm*(1 - sprofile[25 + 27*(j - 1)]);
					}
				}
			}
			if (i > 0) {
				if (*(o + (i - 1)*ls + j) < 0) {
					if (*(o + (i - 1)*ls + j) == -1) {
						gs = GEnorm*(1 - pprofile[23 + 27*i]) - GOnorm*(1 - pprofile[24 + 27*(i - 1)]) + GOnorm*(1 - pprofile[25 + 27*i]);
					} else {
						gs = GEnorm*(1 - pprofile[23 + 27*i]) + GOnorm*(1 - pprofile[25 + 27*i]) - GOnorm*(1 - pprofile[25 + 27*(i - 1)]);
					}
				}
			}
			
			// apply end-gap penalties at corners of matrix
			if (i==0 && j==0) {
				gp = egpL*(1 - sprofile[23 + 27*j])*sprofile[26 + 27*j];
				gs = egpL*(1 - pprofile[23 + 27*i])*pprofile[26 + 27*i];
			} else if (i==(lp - 1) && j==(ls - 1)) {
				gp = egpR*(1 - sprofile[23 + 27*j])*sprofile[26 + 27*j];
				gs = egpR*(1 - pprofile[23 + 27*i])*pprofile[26 + 27*i];
			}
			
			// score aligning letters
			M = *(m + i*(ls + 1) + j);
			for (int n = 0; n < 20; n++) { // omit U in 20, O in 21
				if (pprofile[n + 27*i]!=0) {
					for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
						if (sprofile[p + 27*j]!=0)
							M += pprofile[n + 27*i]*sprofile[p + 27*j] * *(subM + n*21 + p)*tot;
					}
				}
			}
			
			// score aligning stops
			if (pprofile[22 + 27*i]!=0) { // stop (*)
				for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
					if (sprofile[p + 27*j]!=0)
						M += pprofile[22 + 27*i]*sprofile[p + 27*j] * *(subM + 20*21 + p)*tot;
				}
			}
			if (sprofile[22 + 27*j]!=0) { // stop (*)
				for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
					if (pprofile[p + 27*i]!=0)
						M += pprofile[p + 27*i]*sprofile[22 + 27*j] * *(subM + 20*21 + p)*tot;
				}
			}
			if (pprofile[22 + 27*i]!=0 && sprofile[22 + 27*j]!=0) // both stops (*)
				M += pprofile[22 + 27*i]*sprofile[22 + 27*j] * *(subM + 20*21 + 20)*tot;
			
			// score aligning gaps
			M += PMnorm * pprofile[23 + 27*i]*sprofile[23 + 27*j]; // gap-gap match
			if (pprofile[23 + 27*i] > 0) // gap
				M += MMnorm * pprofile[23 + 27*i]*(1 - sprofile[23 + 27*j]); // gap-letter match
			if (sprofile[22 + 25*j] > 0) // gap
				M += MMnorm * (1 - pprofile[23 + 27*i])*sprofile[23 + 27*j]; // gap-letter match
			
			GP = *(m + (i + 1)*(ls + 1) + j) + gp + PMnorm*sprofile[23 + 27*j] - PMnorm*pprofile[23 + 27*i] + MMnorm*(1 - sprofile[23 + 27*j]);
			GS = *(m + i*(ls + 1) + j + 1) + gs + PMnorm*pprofile[23 + 27*i] - PMnorm*sprofile[23 + 27*j] + MMnorm*(1 - pprofile[23 + 27*i]);
			
			if ((GS >= M) && (GS >= GP) && sprofile[23 + 27*j] != 1) {
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
			} else if (GP >= M && pprofile[23 + 27*i] != 1) {
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
			
			//#pragma omp critical // not necessary because will only underestimate max
			if (*(m + (i + 1)*(ls + 1) + j + 1) > max)
				max = *(m + (i + 1)*(ls + 1) + j + 1);
		}
	}
	
	if (egpR != 0) {
		GP = egpR*(1 - sprofile[23 + 27*(ls - 1)])*sprofile[26 + 27*(ls - 1)];
		*(m + lp*(ls + 1) + ls - 1) = *(m + lp*(ls + 1) + ls - 1) + GP;
		for (i = ls - 2; i > -1; i--) {
			GP += egpR*(1 - sprofile[23 + 27*i])*sprofile[26 + 27*i];
			*(m + lp*(ls + 1) + i) = *(m + lp*(ls + 1) + i) + GP;
		}
		GP = egpR*(1 - pprofile[23 + 27*(lp - 1)])*pprofile[26 + 27*(lp - 1)];
		*(m + (lp - 1)*(ls + 1) + ls) = *(m + (lp - 1)*(ls + 1) + ls) + GP;
		for (i = lp - 2; i > -1; i--) {
			GP += egpR*(1 - pprofile[23 + 27*i])*pprofile[26 + 27*i];
			*(m + i*(ls + 1) + ls) = *(m + i*(ls + 1) + ls) + GP;
		}
	}
	
	// find the max scoring alignment
	int maxp = 0;
	int maxs = 0;
	for (i = 1; i <= ls; i++)
		if (*(m + lp*(ls + 1) + i) >= *(m + lp*(ls + 1) + maxs))
			maxs = i;
	for (i = top + 1; i <= lp; i++)
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
