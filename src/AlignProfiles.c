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

//ans_start <- .Call("alignProfiles", p.profile, s.profile, pm, mm, go, ge, endGapPenalty, boundary, processors, PACKAGE="DECIPHER")//, codon, subMatrix
SEXP alignProfiles(SEXP p, SEXP s, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)//, SEXP codonBonus, SEXP subMatrix
{
	int i, j, k, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, S, M, GP, GS;
	double max;
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
	//double *subM = REAL(subMatrix);
	int NTHREADS, nthreads = asInteger(nThreads);
	
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
		
		START = start;
		for (i = top; i <= end; i++) {
			if (i < start) {
				top = start;
				continue;
			}
			
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
				top = i;
				//*(m + (i + 1)*(ls + 1) + j + 1) = *(m + i*(ls + 1) + j) + GE;
				*(o + i*ls + j) = 1;
			} else {
				if (i==0) {
					START = i;
					top = i;
				} else { // i > 0
					if (i - 1 < start) {
						START = start;
						top = start;
					} else {
						START = i - 1;
						top = START;
					}
				}
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
			//Rprintf("\ni%d j%d left%d END%d value%d", i, j, left, END, (int)*(m + i*(ls + 1) + j));
			if (j < left) {
				i -= left - j;
				continue;
			} else if (*(m + i*(ls + 1) + j) < max) {
				left = j;
				//*(m + (i + 1)*(ls + 1) + j + 1) = *(m + i*(ls + 1) + j) + GE;
				*(o + i*ls + j) = -1;
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
		//Rprintf("\nk%d max%d start%d end%d START%d END %d left %d", k, (int)max, start, end, START, END, left);
		
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
			S = 0;
			M = 0;
			if (MM==0) {
				for (int n = 0; n < 5; n++) {
					if (pprofile[n + 7*i]==0)
						continue;
					S += sqrt(pprofile[n + 7*i]*sprofile[n + 7*j]);
				}
				M = *(m + i*(ls + 1) + j) + S*PM;
			} else {
				for (int n = 0; n < 5; n++) {
					if (pprofile[n + 7*i]==0)
						continue;
					S += sqrt(pprofile[n + 7*i]*sprofile[n + 7*j]);
					for (int p = 0; p < 5; p++) {
						if (p==n || sprofile[p + 7*j]==0)
							continue;
						M += sqrt(pprofile[n + 7*i]*sprofile[p + 7*j]);
					}
				}
				M = *(m + i*(ls + 1) + j) + S*PM + M*MM;
			}
			
			//S = sqrt(pprofile[7*i]*sprofile[7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[4 + 7*j]);
			//M = *(m + i*(ls + 1) + j) + (sqrt(pprofile[0 + 7*i]*sprofile[1 + 7*j]) * *(subM + 1*4 + 0) + sqrt(pprofile[0 + 7*i]*sprofile[2 + 7*j]) * *(subM + 2*4 + 0) + sqrt(pprofile[0 + 7*i]*sprofile[3 + 7*j]) * *(subM + 3*4 + 0) + sqrt(pprofile[1 + 7*i]*sprofile[0 + 7*j]) * *(subM + 0*4 + 1) + sqrt(pprofile[1 + 7*i]*sprofile[2 + 7*j]) * *(subM + 2*4 + 1) + sqrt(pprofile[1 + 7*i]*sprofile[3 + 7*j]) * *(subM + 3*4 + 1) + sqrt(pprofile[2 + 7*i]*sprofile[0 + 7*j]) * *(subM + 0*4 + 2) + sqrt(pprofile[2 + 7*i]*sprofile[1 + 7*j]) * *(subM + 1*4 + 2) + sqrt(pprofile[2 + 7*i]*sprofile[3 + 7*j]) * *(subM + 3*4 + 2) + sqrt(pprofile[3 + 7*i]*sprofile[0 + 7*j]) * *(subM + 0*4 + 3) + sqrt(pprofile[3 + 7*i]*sprofile[1 + 7*j]) * *(subM + 1*4 + 3) + sqrt(pprofile[3 + 7*i]*sprofile[2 + 7*j]) * *(subM + 2*4 + 3) + sqrt(pprofile[4 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[4 + 7*j]))*MM + S*PM;
			//M = *(m + i*(ls + 1) + j) + (sqrt(pprofile[0 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[0 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[1 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[2 + 7*j]) + sqrt(pprofile[4 + 7*i]*sprofile[3 + 7*j]) + sqrt(pprofile[0 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[1 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[2 + 7*i]*sprofile[4 + 7*j]) + sqrt(pprofile[3 + 7*i]*sprofile[4 + 7*j]))*MM + S*PM;
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
			
			if ((GS >= M) && (GS >= GP) && sprofile[4 + 7*j] != 1) {
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
			} else if (GP >= M && pprofile[4 + 7*i] != 1) {
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

//ans_start <- .Call("alignProfilesAA", p.profile, s.profile, subMatrix, pm, mm, go, ge, endGapPenalty, boundary, processors, PACKAGE="DECIPHER")
SEXP alignProfilesAA(SEXP p, SEXP s, SEXP subMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)
{
	int i, j, k, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, S, M, GP, GS;
	double max;
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
	int *subM = INTEGER(subMatrix);
	int NTHREADS, nthreads = asInteger(nThreads);
	
	R_len_t lp = length(p)/26;
	R_len_t ls = length(s)/26;
	
	PROTECT(ans = allocVector(INTSXP, lp + ls + 4));
	rans = INTEGER(ans);
	
	float *m = Calloc((lp+1)*(ls+1), float); // initialized to zero
	int *o = Calloc(lp*ls, int); // initialized to zero
	
	if (egpL != 0) {
		*(m + 1) = egpL*(1 - sprofile[23]);
		for (i = 2; i <= ls; i++)
			*(m + i) = *(m + i - 1) + egpL*(1 - sprofile[23 + 26*(i - 1)]);
		*(m + ls + 1) = egpL*(1 - pprofile[23]);
		for (i = 2; i <= lp; i++)
			*(m + (ls + 1)*i) = *(m + (ls + 1)*(i - 1)) + egpL*(1 - pprofile[23 + 26*(i - 1)]);
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
		
		START = start;
		for (i = top; i <= end; i++) {
			if (i < start) {
				top = start;
				continue;
			}
			
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
				top = i;
				//*(m + (i + 1)*(ls + 1) + j + 1) = *(m + i*(ls + 1) + j) + GE;
				*(o + i*ls + j) = 1;
			} else {
				if (i==0) {
					START = i;
					top = i;
				} else { // i > 0
					if (i - 1 < start) {
						START = start;
						top = start;
					} else {
						START = i - 1;
						top = START;
					}
				}
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
			} else if (*(m + i*(ls + 1) + j) < max) {
				left = j;
				//*(m + (i + 1)*(ls + 1) + j + 1) = *(m + i*(ls + 1) + j) + GE;
				*(o + i*ls + j) = -1;
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
			
			// determine insertion penalty (gap extension or gap opening)
			gp = 2*GO*(1 - sprofile[24 + 26*j]);
			gs = 2*GO*(1 - pprofile[24 + 26*i]);
			if (j > 0) {
				if (*(o + i*ls + j - 1) > 0) {
					if (*(o + i*ls + j - 1) == 1) {
						gp = GE*(1 - sprofile[23 + 26*j]) - GO*(1 - sprofile[24 + 26*(j - 1)]) + GO*(1 - sprofile[25 + 26*j]);
					} else {
						gp = GE*(1 - sprofile[23 + 26*j]) + GO*(1 - sprofile[25 + 26*j]) - GO*(1 - sprofile[25 + 26*(j - 1)]);
					}
				}
			}
			if (i > 0) {
				if (*(o + (i - 1)*ls + j) < 0) {
					if (*(o + (i - 1)*ls + j) == -1) {
						gs = GE*(1 - pprofile[23 + 26*i]) - GO*(1 - pprofile[24 + 26*(i - 1)]) + GO*(1 - pprofile[25 + 26*i]);
					} else {
						gs = GE*(1 - pprofile[23 + 26*i]) + GO*(1 - pprofile[25 + 26*i]) - GO*(1 - pprofile[25 + 26*(i - 1)]);
					}
				}
			}
			
			// score aligning letters
			M = *(m + i*(ls + 1) + j);
			for (int n = 0; n < 20; n++) { // omit U in 20, O in 21
				if (pprofile[n + 26*i]!=0) {
					for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
						if (sprofile[p + 26*j]!=0)
							M += sqrt(pprofile[n + 26*i]*sprofile[p + 26*j]) * *(subM + n*21 + p);
					}
				}
			}
			
			// score aligning stops
			if (pprofile[22 + 26*i]!=0) { // stop (*)
				for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
					if (sprofile[p + 26*j]!=0)
						M += sqrt(pprofile[22 + 26*i]*sprofile[p + 26*j]) * *(subM + 20*21 + p);
				}
			}
			if (sprofile[22 + 26*j]!=0) { // stop (*)
				for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
					if (pprofile[p + 26*i]!=0)
						M += sqrt(pprofile[p + 26*i]*sprofile[22 + 26*j]) * *(subM + 20*21 + p);
				}
			}
			if (pprofile[22 + 26*i]!=0 && sprofile[22 + 26*j]!=0) // both stops (*)
				M += sqrt(pprofile[22 + 26*i]*sprofile[22 + 26*j]) * *(subM + 20*21 + 20);
			
			// score aligning gaps
			M += PM * sqrt(pprofile[23 + 26*i]*sprofile[23 + 26*j]); // gap-gap match
			if (pprofile[23 + 26*i] > 0) // gap
				M += MM * sqrt(pprofile[23 + 26*i]*(1 - sprofile[23 + 26*j])); // gap-letter match
			if (sprofile[22 + 25*j] > 0) // gap
				M += MM * sqrt((1 - pprofile[23 + 26*i])*sprofile[23 + 26*j]); // gap-letter match
			
			/*
			// multiply match score by fractional information content
			double R = 8.643856; // max entropy = 2*log2(20)
			for (int n = 0; n < 20; n++) // omit U in 20
				if (pprofile[n + 25*i]!=0)
					R += pprofile[n + 25*i]*log2(pprofile[n + 25*i]);
			for (int n = 0; n < 20; n++) // omit U in 20
				if (sprofile[n + 25*j]!=0)
					R += sprofile[n + 25*j]*log2(sprofile[n + 25*j]);
			M *= cbrt(R/8.643856);
			*/
			
			GP = *(m + (i + 1)*(ls + 1) + j) + gp + PM*sprofile[23 + 26*j] - PM*pprofile[23 + 26*i];
			GS = *(m + i*(ls + 1) + j + 1) + gs + PM*pprofile[23 + 26*i] - PM*sprofile[23 + 26*j];
			
			if ((GS >= M) && (GS >= GP) && sprofile[23 + 26*j] != 1) {
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
			} else if (GP >= M && pprofile[23 + 26*i] != 1) {
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
		GP = egpR*(1 - sprofile[23 + 26*(ls - 1)]);
		*(m + lp*(ls + 1) + ls - 1) = *(m + lp*(ls + 1) + ls - 1) + GP;
		for (i = ls - 2; i > -1; i--) {
			GP += egpR*(1 - sprofile[23 + 26*i]);
			*(m + lp*(ls + 1) + i) = *(m + lp*(ls + 1) + i) + GP;
		}
		GP = egpR*(1 - pprofile[23 + 26*(lp - 1)]);
		*(m + (lp - 1)*(ls + 1) + ls) = *(m + (lp - 1)*(ls + 1) + ls) + GP;
		for (i = lp - 2; i > -1; i--) {
			GP += egpR*(1 - pprofile[23 + 26*i]);
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
