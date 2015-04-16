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

SEXP alignProfiles(SEXP p, SEXP s, SEXP subMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)
{
	int i, j, k, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, S, M, GP, GS;
	double max, tot;
	SEXP ans;
	
	pprofile = REAL(p);
	sprofile = REAL(s);
	double PM = asReal(pm);
	double MM = asReal(mm);
	double GO = asReal(go);
	double GE = asReal(ge);
	double EX = asReal(exp);
	double POW = asReal(power);
	double bound = asReal(boundary);
	double egpL = asReal(endGapPenaltyLeft);
	double egpR = asReal(endGapPenaltyRight);
	
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
	
	// zipfian distribution for gap cost
	if (lp > ls) {
		j = lp;
	} else {
		j = ls;
	}
	double *zip = Calloc(j + 1, double); // initialized to zero
	for (i = 1; i <= j; i++)
		zip[i] = pow((double)i, EX);
	
	// deterine the fraction of sequences starting or ending
	double *pstarts = Calloc(lp, double); // initialized to zero
	double *pstops = Calloc(lp, double); // initialized to zero
	double *sstarts = Calloc(ls, double); // initialized to zero
	double *sstops = Calloc(ls, double); // initialized to zero
	double past, current;
	past = 0;
	for (i = 0; i < lp; i++) {
		current = pprofile[7 + 8*i];
		if (current > past) {
			pstarts[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i > 0)
			pstarts[i] += pstarts[i - 1];
	}
	past = 0;
	for (i = lp - 1; i >= 0; i--) {
		current = pprofile[7 + 8*i];
		if (current > past) {
			pstops[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i < (lp - 1))
			pstops[i] += pstops[i + 1];
	}
	past = 0;
	for (i = 0; i < ls; i++) {
		current = sprofile[7 + 8*i];
		if (current > past) {
			sstarts[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i > 0)
			sstarts[i] += sstarts[i - 1];
	}
	past = 0;
	for (i = ls - 1; i >= 0; i--) {
		current = sprofile[7 + 8*i];
		if (current > past) {
			sstops[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i < (ls - 1))
			sstops[i] += sstops[i + 1];
	}
	
	// normalize profile positions to coverge^POW
	for (i = 0; i < lp; i++)
		if (pprofile[7 + 8*i] > 0)
			pprofile[7 + 8*i] = pow(pprofile[7 + 8*i], POW);
	for (i = 0; i < ls; i++)
		if (sprofile[7 + 8*i] > 0)
			sprofile[7 + 8*i] = pow(sprofile[7 + 8*i], POW);
	
	if (egpL != 0) {
		*(m + 1) = egpL*sstarts[0];
		for (i = 2; i <= ls; i++)
			*(m + i) = egpL*sstarts[i - 1];
		*(m + ls + 1) = egpL*pstarts[0];
		for (i = 2; i <= lp; i++)
			*(m + (ls + 1)*i) = egpL*pstarts[i - 1];
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
		#pragma omp parallel for private(i,j,gp,gs,S,M,GP,GS,tot) num_threads(NTHREADS)
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
			
			// apply end-gap penalties at corners of matrix
			if (i==0 && j==0) {
				gp = egpL*sstarts[0];
				gs = egpL*pstarts[0];
			} else if (i==(lp - 1) && j==(ls - 1)) {
				gp = egpR*sstops[ls - 1];
				gs = egpR*pstops[lp - 1];
			} else {
				// determine insertion penalty (gap extension or gap opening)
				if (j > 0 && *(o + i*ls + j - 1) > 0) {
					gp = GE*(1 - 0.3*sprofile[4 + 8*j])*zip[*(o + i*ls + j - 1)] + 0.5*GO*(sprofile[6 + 8*(j - 1)] - sprofile[6 + 8*j]);
				} else {
					gp = GO*(2 - 0.5*sprofile[5 + 8*j] - 0.5*sprofile[6 + 8*j]);
				}
				if (i > 0 && *(o + (i - 1)*ls + j) < 0) {
					gs = GE*(1 - 0.3*pprofile[4 + 8*i])*zip[*(o + (i - 1)*ls + j)*-1] + 0.5*GO*(pprofile[6 + 8*(i - 1)] - pprofile[6 + 8*i]);
				} else {
					gs = GO*(2 - 0.5*pprofile[5 + 8*i] - 0.5*pprofile[6 + 8*i]);
				}
			}
			
			tot = (1 - pprofile[4 + 8*i])*(1 - sprofile[4 + 8*j]); // max of dot-product
			GS = gs*tot;
			GP = gp*tot;
			
			if (sprofile[7 + 8*j] > 0 && pprofile[7 + 8*i] > 0) {
				tot = sqrt(sprofile[7 + 8*j]*pprofile[7 + 8*i]); // normalization factor
			} else {
				tot = 0;
			}
			
			// score aligning
			if (do_subM) {
				S = 0;
				M = 0;
				for (int n = 0; n < 4; n++) {
					if (pprofile[n + 8*i]==0)
						continue;
					S += pprofile[n + 8*i]*sprofile[n + 8*j] * *(subM + n*4 + n);
					for (int p = 0; p < 4; p++) {
						if (p==n || sprofile[p + 8*j]==0)
							continue;
						M += pprofile[n + 8*i]*sprofile[p + 8*j] * *(subM + n*4 + p);
					}
				}
				M = S + M;
			} else { // no substitution matrix
				S = 0;
				M = 0;
				for (int n = 0; n < 4; n++) {
					if (pprofile[n + 8*i]==0)
						continue;
					S += pprofile[n + 8*i]*sprofile[n + 8*j];
				}
				M = S*PM + ((1 - pprofile[4 + 8*i])*(1 - sprofile[4 + 8*j]) - S)*MM; // % mismatched = (% ungapped) - (% matched)
			}
			
			if (((GS + *(m + i*(ls + 1) + j + 1)) > (M + *(m + i*(ls + 1) + j))) && ((GS + *(m + i*(ls + 1) + j + 1)) > (GP + *(m + (i + 1)*(ls + 1) + j))) && sprofile[4 + 8*j] != 1) {
				if (i > 0) {
					if (*(o + (i - 1)*ls + j) < 0) {
						*(o + i*ls + j) = *(o + (i - 1)*ls + j) - 1;
					} else {
						*(o + i*ls + j) = -1;
					}
				} else {
					*(o + i*ls + j) = -1;
				}
				
				*(m + (i + 1)*(ls + 1) + j + 1) = GS*tot + *(m + i*(ls + 1) + j + 1);
			} else if ((GP + *(m + (i + 1)*(ls + 1) + j)) > (M + *(m + i*(ls + 1) + j)) && pprofile[4 + 8*i] != 1) {
				if (j > 0) {
					if (*(o + i*ls + j - 1) > 0) {
						*(o + i*ls + j) = *(o + i*ls + j - 1) + 1;
					} else {
						*(o + i*ls + j) = 1;
					}
				} else {
					*(o + i*ls + j) = 1;
				}
				
				*(m + (i + 1)*(ls + 1) + j + 1) = GP*tot + *(m + (i + 1)*(ls + 1) + j);
			} else {
				*(m + (i + 1)*(ls + 1) + j + 1) = M*tot + *(m + i*(ls + 1) + j);
				*(o + i*ls + j) = 0;
			}
			
			//#pragma omp critical // not necessary because will only underestimate max
			if (*(m + (i + 1)*(ls + 1) + j + 1) > max)
				max = *(m + (i + 1)*(ls + 1) + j + 1);
		}
	}
	
	if (egpR != 0) {
		GP = egpR*sstops[ls - 1];
		*(m + lp*(ls + 1) + ls - 1) = *(m + lp*(ls + 1) + ls - 1) + GP;
		for (i = ls - 2; i > -1; i--) {
			GP = egpR*sstops[i];
			*(m + lp*(ls + 1) + i) = *(m + lp*(ls + 1) + i) + GP;
		}
		GP = egpR*pstops[lp - 1];
		*(m + (lp - 1)*(ls + 1) + ls) = *(m + (lp - 1)*(ls + 1) + ls) + GP;
		for (i = lp - 2; i > -1; i--) {
			GP = egpR*pstops[i];
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
	Free(zip);
	Free(pstarts);
	Free(pstops);
	Free(sstarts);
	Free(sstops);
	
	UNPROTECT(1);
	
	return ans; // vector = [end.p end.s start.p start.s count t]
}

SEXP alignProfilesAA(SEXP p, SEXP s, SEXP subMatrix, SEXP hecMatrix, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)
{
	int i, j, k, pos, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, M, GP, GS;
	double max, tot, freq;
	SEXP ans, dims;
	
	pprofile = REAL(p);
	sprofile = REAL(s);
	double GO = asReal(go);
	double GE = asReal(ge);
	double EX = asReal(exp);
	double POW = asReal(power);
	double bound = asReal(boundary);
	double egpL = asReal(endGapPenaltyLeft);
	double egpR = asReal(endGapPenaltyRight);
	double *subM = REAL(subMatrix);
	int NTHREADS, nthreads = asInteger(nThreads);
	
	double *hecM = REAL(hecMatrix);
	int do_HEC, d, size = 27;
	if (length(hecMatrix) > 0) {
		do_HEC = 1;
		PROTECT(dims = GET_DIM(hecMatrix));
		d = INTEGER(dims)[0];
		UNPROTECT(1);
		size += d;
	} else {
		do_HEC = 0;
		d = 0;
	}
	
	R_len_t lp = length(p)/size;
	R_len_t ls = length(s)/size;
	
	int N21[20] = {0, 21, 42, 63, 84, 105, 126, 147, 168, 189, 210, 231, 252, 273, 294, 315, 336, 357, 378, 399};
	
	// initialize arrays of residue orders
	int *Op = Calloc(lp*20, int); // initialized to zero
	int *Os = Calloc(ls*20, int); // initialized to zero
	int SIZE20, ISIZE;
	
	for (i = 0; i < lp; i++) {
		k = 0;
		z = -1;
		SIZE20 = 20*i;
		ISIZE = size*i;
		for (j = 0; j < 20; j++) {
			if (pprofile[j + ISIZE] > 0) {
				Op[k + SIZE20] = j;
				k++;
			} else if (z==-1) {
				z = j;
			}
		}
		if (z!=-1)
			Op[k + SIZE20] = z;
	}
	for (i = 0; i < ls; i++) {
		k = 0;
		z = -1;
		SIZE20 = 20*i;
		ISIZE = size*i;
		for (j = 0; j < 20; j++) {
			if (sprofile[j + ISIZE] > 0) {
				Os[k + SIZE20] = j;
				k++;
			} else if (z==-1) {
				z = j;
			}
		}
		if (z!=-1)
			Os[k + SIZE20] = z;
	}
	
	// gap modulation factor based on position(s) across from the gap
	// added once for each gap added (including opening/closing)
	// A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V:
	double GC[20] = {0.101,-0.453,0.446,0.178,-0.48,0.658,0.245,0.702,0.221,-1.586,-0.979,0.239,-1.574,-1.097,1.188,1.364,0.417,-1.379,-1.486,-0.788};
	
	// gap modulation factor based on +/- positions of gap
	// half of log-odds because added twice (gap open and close)
	// units in third-bits - may not scale with all substitution matrices
	double GM[320] = { // gap modulate position -8 to +8
		0.111,0.109,0.088,0.097,0.11,0.117,0.037,-0.67,0.116,0.003,-0.031,-0.05,0.016,-0.017,0.027,-0.002, // A
		0.067,0.155,0.025,0.001,0.019,0.017,-0.08,0.063,-0.166,-0.235,-0.268,-0.15,-0.207,-0.148,-0.084,-0.073, // R
		-0.123,-0.126,-0.153,-0.207,-0.032,0.048,0.19,0.082,0.564,0.37,0.187,0.079,-0.096,-0.133,-0.191,-0.201, // N
		-0.121,-0.181,-0.164,-0.133,-0.092,0.121,0.195,0.526,0.658,0.579,0.371,0.246,0.15,0.127,0.074,0.071, // D
		0.102,-0.029,-0.085,-0.129,-0.089,-0.363,-0.493,-0.243,-0.728,-0.55,-0.395,-0.193,-0.151,0.049,0.095,0.029, // C
		0.119,0.114,0.08,0.151,0.131,0.273,0.264,0.098,0.246,0.165,0.106,0.091,0.008,0.002,-0.037,-0.048, // Q
		0.076,0.144,0.168,0.134,0.23,0.362,0.361,0.447,0.437,0.368,0.311,0.247,0.153,0.138,0.087,0.059, // E
		-0.219,-0.187,-0.198,-0.231,-0.157,0.031,0.291,0.55,0.605,0.449,0.293,0.086,-0.119,-0.211,-0.194,-0.131, // G
		-0.032,-0.041,0.027,-0.121,-0.035,-0.112,-0.024,-0.302,-0.243,-0.12,-0.132,-0.254,-0.099,-0.187,-0.08,-0.156, // H
		0.038,-0.01,-0.163,-0.11,-0.342,-0.481,-0.693,-0.646,-0.993,-0.802,-0.583,-0.341,-0.142,-0.103,-0.059,0.061, // I
		0.038,-0.013,0.056,0.066,-0.006,-0.228,-0.362,-0.216,-0.895,-0.679,-0.539,-0.416,-0.25,-0.118,-0.155,-0.119, // L
		0.006,0.089,0.111,0.13,0.251,0.302,0.417,0.373,0.276,0.249,0.218,0.138,0.064,-0.046,-0.05,-0.046, // K
		-0.184,-0.252,-0.211,-0.301,-0.481,-0.46,-0.673,-0.61,-1.24,-0.907,-0.91,-0.797,-0.567,-0.564,-0.477,-0.478, // M
		0,-0.094,-0.134,-0.029,-0.159,-0.291,-0.491,-0.229,-1.165,-0.837,-0.533,-0.409,-0.202,-0.137,-0.061,-0.045, // F
		0.003,0.044,0.103,0.128,0.158,0.267,0.441,0.683,0.361,0.57,0.614,0.54,0.494,0.414,0.32,0.205, // P
		0.002,0.02,0.069,0.075,0.17,0.176,0.298,-0.1,0.585,0.412,0.368,0.342,0.221,0.215,0.2,0.188, // S
		-0.046,-0.082,0.005,-0.028,-0.002,0.012,0.006,-0.272,0.203,0.173,0.229,0.212,0.175,0.179,0.122,0.126, // T
		0.089,0.156,-0.07,0.02,-0.012,-0.451,-0.389,-0.109,-1.388,-0.953,-0.617,-0.288,-0.18,-0.145,-0.035,0.043, // W
		-0.039,0.037,0.025,-0.071,-0.15,-0.245,-0.458,-0.224,-1.076,-0.8,-0.561,-0.397,-0.227,-0.167,-0.189,-0.144, // Y
		0.049,0.005,0.04,0.054,-0.131,-0.31,-0.51,-0.505,-0.545,-0.359,-0.2,-0.053,0.101,0.141,0.198,0.168, // V
	};
	
	double *GCp = Calloc(lp, double); // initialized to zero
	double *GCs = Calloc(ls, double); // initialized to zero
	double *GOp = Calloc(lp, double); // initialized to zero
	double *GOs = Calloc(ls, double); // initialized to zero
	
	for (i = 0; i < lp; i++) {
		GOp[i] = GO;
		
		tot = 0;
		for (j = 0; j < 20; j++)
			tot += pprofile[j + size*i];
		
		if (tot > 0) {
			tot = sqrt(tot); // normalize by sqrt of occupancy
			for (j = 0; j < 20; j++) {
				freq = pprofile[j + size*i]/tot;
				if (freq == 0)
					continue;
				
				GCp[i] += freq*GC[j];
				
				for (k = 4; k <= 14; k++) { // -4 to +7
					pos = i + k - 7;
					if (pos < 0 || pos >= lp)
						continue;
					
					GOp[pos] += freq*GM[j*16 + 15 - k];
				}
			}
		}
	}
	for (i = 0; i < ls; i++) {
		GOs[i] = GO;
		
		tot = 0;
		for (j = 0; j < 20; j++)
			tot += sprofile[j + size*i];
		
		if (tot > 0) {
			tot = sqrt(tot); // normalize by sqrt of occupancy
			for (j = 0; j < 20; j++) {
				freq = sprofile[j + size*i]/tot;
				if (freq == 0)
					continue;
				
				GCs[i] += freq*GC[j];
				
				for (k = 4; k <= 14; k++) { // -4 to +7
					pos = i + k - 7;
					if (pos < 0 || pos >= ls)
						continue;
					
					GOs[pos] += freq*GM[j*16 + 15 - k];
				}
			}
		}
	}
	
	PROTECT(ans = allocVector(INTSXP, lp + ls + 4));
	rans = INTEGER(ans);
	
	float *m = Calloc((lp+1)*(ls+1), float); // initialized to zero
	int *o = Calloc(lp*ls, int); // initialized to zero
	
	// zipfian distribution for gap cost
	if (lp > ls) {
		j = lp;
	} else {
		j = ls;
	}
	double *zip = Calloc(j + 1, double); // initialized to zero
	for (i = 1; i <= j; i++)
		zip[i] = pow((double)i, EX);
	
	// deterine the fraction of sequences starting or ending
	double *pstarts = Calloc(lp, double); // initialized to zero
	double *pstops = Calloc(lp, double); // initialized to zero
	double *sstarts = Calloc(ls, double); // initialized to zero
	double *sstops = Calloc(ls, double); // initialized to zero
	double past, current;
	past = 0;
	for (i = 0; i < lp; i++) {
		current = pprofile[26 + size*i];
		if (current > past) {
			pstarts[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i > 0)
			pstarts[i] += pstarts[i - 1];
	}
	past = 0;
	for (i = lp - 1; i >= 0; i--) {
		current = pprofile[26 + size*i];
		if (current > past) {
			pstops[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i < (lp - 1))
			pstops[i] += pstops[i + 1];
	}
	past = 0;
	for (i = 0; i < ls; i++) {
		current = sprofile[26 + size*i];
		if (current > past) {
			sstarts[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i > 0)
			sstarts[i] += sstarts[i - 1];
	}
	past = 0;
	for (i = ls - 1; i >= 0; i--) {
		current = sprofile[26 + size*i];
		if (current > past) {
			sstops[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i < (ls - 1))
			sstops[i] += sstops[i + 1];
	}
	
	// normalize profile positions to coverge^POW
	for (i = 0; i < lp; i++)
		if (pprofile[26 + size*i] > 0)
			pprofile[26 + size*i] = pow(pprofile[26 + size*i], POW);
	for (i = 0; i < ls; i++)
		if (sprofile[26 + size*i] > 0)
			sprofile[26 + size*i] = pow(sprofile[26 + size*i], POW);
	
	if (egpL != 0) {
		*(m + 1) = egpL*sstarts[0];
		for (i = 2; i <= ls; i++)
			*(m + i) = egpL*sstarts[i - 1];
		*(m + ls + 1) = egpL*pstarts[0];
		for (i = 2; i <= lp; i++)
			*(m + (ls + 1)*i) = egpL*pstarts[i - 1];
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
		#pragma omp parallel for private(i,j,gp,gs,M,GP,GS,tot) num_threads(NTHREADS)
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
			
			int SIZEI = size*i;
			int SIZEJ = size*j;
			int SIZE20I = 20*i;
			int SIZE20J = 20*j;
			
			// apply end-gap penalties at corners of matrix
			if (i==0 && j==0) {
				gp = egpL*sstarts[0];
				gs = egpL*pstarts[0];
			} else if (i==(lp - 1) && j==(ls - 1)) {
				gp = egpR*sstops[ls - 1];
				gs = egpR*pstops[lp - 1];
			} else {
				// determine insertion penalty (gap extension or gap opening)
				if (j > 0 && *(o + i*ls + j - 1) > 0) {
					gp = GE*(1 - 0.3*sprofile[23 + SIZEJ])*zip[*(o + i*ls + j - 1)] + 0.5*GOp[i]*(sprofile[25 + size*(j - 1)] - sprofile[25 + SIZEJ]);
				} else {
					gp = GOp[i]*(2 - 0.5*(sprofile[24 + SIZEJ] + sprofile[25 + SIZEJ]));
				}
				if (i > 0 && *(o + (i - 1)*ls + j) < 0) {
					gs = GE*(1 - 0.3*pprofile[23 + SIZEI])*zip[*(o + (i - 1)*ls + j)*-1] + 0.5*GOs[j]*(pprofile[25 + size*(i - 1)] - pprofile[25 + SIZEI]);
				} else {
					gs = GOs[j]*(2 - 0.5*(pprofile[24 + SIZEI] + pprofile[25 + SIZEI]));
				}
				
				// add modulation factor for center gap region
				// amount is based on the opposing sequence
				gp += GCs[j]*(1 - sprofile[23 + SIZEJ]);
				gs += GCp[i]*(1 - pprofile[23 + SIZEI]);
			}
			
			tot = (1 - pprofile[23 + SIZEI])*(1 - sprofile[23 + SIZEJ]); // max of dot-product
			GS = gs*tot;
			GP = gp*tot;
			
			if (sprofile[26 + SIZEJ] > 0 && pprofile[26 + SIZEI] > 0) {
				tot = sqrt(sprofile[26 + SIZEJ]*pprofile[26 + SIZEI]); // normalization factor
			} else {
				tot = 0;
			}
			
			M = 0;
			double pFreq, sFreq;
			// score aligning letters
			for (int n = 0; n < 20; n++) { // omit U in 20, O in 21
				pFreq = pprofile[Op[n + SIZE20I] + SIZEI];
				if (pFreq!=0) {
					for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
						sFreq = sprofile[Os[p + SIZE20J] + SIZEJ];
						if (sFreq!=0) {
							M += pFreq*sFreq * *(subM + N21[Op[n + SIZE20I]] + Os[p + SIZE20J]);
						} else {
							break;
						}
					}
				} else {
					break;
				}
			}
			
			if (do_HEC) {
				for (int n = 0; n < d; n++) {
					for (int p = 0; p < d; p++) {
						M += pprofile[n + 27 + SIZEI]*sprofile[p + 27 + SIZEJ] * *(hecM + n*d + p);
					}
				}
			}
			
			// score aligning stops
			if (pprofile[22 + SIZEI]!=0) { // stop (*)
				for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
					if (sprofile[p + SIZEJ]!=0)
						M += pprofile[22 + SIZEI]*sprofile[p + SIZEJ] * *(subM + 420 + p);
				}
				if (sprofile[22 + SIZEJ]!=0) // both stops (*)
					M += pprofile[22 + SIZEI]*sprofile[22 + SIZEJ] * *(subM + 440);
			}
			if (sprofile[22 + SIZEJ]!=0) { // stop (*)
				for (int p = 0; p < 20; p++) { // omit U in 20, O in 21
					if (pprofile[p + SIZEI]!=0)
						M += pprofile[p + SIZEI]*sprofile[22 + SIZEJ] * *(subM + 420 + p);
				}
			}
			
			if ((GS + *(m + i*(ls + 1) + j + 1) > M + *(m + i*(ls + 1) + j)) && (GS + *(m + i*(ls + 1) + j + 1) > GP + *(m + (i + 1)*(ls + 1) + j)) && sprofile[23 + SIZEJ] != 1) {
				if (i > 0) {
					if (*(o + (i - 1)*ls + j) < 0) {
						*(o + i*ls + j) = *(o + (i - 1)*ls + j) - 1;
					} else {
						*(o + i*ls + j) = -1;
					}
				} else {
					*(o + i*ls + j) = -1;
				}
				
				*(m + (i + 1)*(ls + 1) + j + 1) = GS*tot + *(m + i*(ls + 1) + j + 1);
			} else if (GP + *(m + (i + 1)*(ls + 1) + j) > M + *(m + i*(ls + 1) + j) && pprofile[23 + SIZEI] != 1) {
				if (j > 0) {
					if (*(o + i*ls + j - 1) > 0) {
						*(o + i*ls + j) = *(o + i*ls + j - 1) + 1;
					} else {
						*(o + i*ls + j) = 1;
					}
				} else {
					*(o + i*ls + j) = 1;
				}
				
				*(m + (i + 1)*(ls + 1) + j + 1) = GP*tot + *(m + (i + 1)*(ls + 1) + j);
			} else {
				*(m + (i + 1)*(ls + 1) + j + 1) = M*tot + *(m + i*(ls + 1) + j);
				*(o + i*ls + j) = 0;
			}
			
			//#pragma omp critical // not necessary because will only underestimate max
			if (*(m + (i + 1)*(ls + 1) + j + 1) > max)
				max = *(m + (i + 1)*(ls + 1) + j + 1);
		}
	}
	
	if (egpR != 0) {
		GP = egpR*sstops[ls - 1];
		*(m + lp*(ls + 1) + ls - 1) = *(m + lp*(ls + 1) + ls - 1) + GP;
		for (i = ls - 2; i > -1; i--) {
			GP = egpR*sstops[i];
			*(m + lp*(ls + 1) + i) = *(m + lp*(ls + 1) + i) + GP;
		}
		GP = egpR*pstops[lp - 1];
		*(m + (lp - 1)*(ls + 1) + ls) = *(m + (lp - 1)*(ls + 1) + ls) + GP;
		for (i = lp - 2; i > -1; i--) {
			GP = egpR*pstops[i];
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
		//*(m + (i + 1)*(ls + 1) + j + 1) = 0;
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
	Free(Op);
	Free(Os);
	Free(zip);
	Free(GCp);
	Free(GCs);
	Free(GOp);
	Free(GOs);
	Free(pstarts);
	Free(pstops);
	Free(sstarts);
	Free(sstops);
	
	UNPROTECT(1);
	
	return ans; // vector = [end.p end.s start.p start.s count t]
}
