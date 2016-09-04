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

SEXP alignProfiles(SEXP p, SEXP s, SEXP type, SEXP subMatrix, SEXP dbnMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)
{
	int i, j, k, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, S, M, GP, GS;
	double max, tot;
	SEXP ans1, ans2, ans3, ans4, dims;
	
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
	int RNA = asInteger(type);
	
	double *subM;
	int do_subM;
	if (length(subMatrix)==0) {
		do_subM = 0;
	} else {
		do_subM = 1;
		subM = REAL(subMatrix);
	}
	int NTHREADS, nthreads = asInteger(nThreads);
	
	double *dbnM = REAL(dbnMatrix);
	int do_DBN, d, size = 8;
	if (length(dbnMatrix) > 0) {
		do_DBN = 1;
		PROTECT(dims = GET_DIM(dbnMatrix));
		d = INTEGER(dims)[0];
		UNPROTECT(1);
		size += d;
	} else {
		do_DBN = 0;
		d = 0;
	}
	
	R_len_t lp = length(p)/size;
	R_len_t ls = length(s)/size;
	int l = lp + ls;
	
	// create array of normalization factors (coverge^POW)
	double *pnorm = Calloc(lp, double); // initialized to zero
	double *snorm = Calloc(ls, double); // initialized to zero
	for (i = 0; i < lp; i++)
		if (pprofile[7 + size*i] > 0)
			pnorm[i] = pow(pprofile[7 + size*i], POW);
	for (i = 0; i < ls; i++)
		if (sprofile[7 + size*i] > 0)
			snorm[i] = pow(sprofile[7 + size*i], POW);
	
	// substitution matrix for pairs of nucleotides
	double subD[256] = {
		+0.2,-0.5,-0.4,-0.7,-0.7,-1.1,-0.9,-1.3,-0.4,-1.2,-1.5,-1.4,-0.6,-0.9,-1.4,-0.7, // AA
		-0.5,+0.5,-0.8,-0.4,-0.2,-0.8,-0.2,-1.0,-0.8,-0.2,-1.1,-1.0,-0.5,-0.8,-1.1,-1.4, // AC
		-0.4,-0.8,+0.5,-0.8,-0.9,-1.1,-0.4,-1.2,-0.8,-0.9,-0.2,-0.9,-0.9,-1.1,-1.0,-1.2, // AG
		-0.7,-0.4,-0.8,+0.2,-0.8,-1.6,-0.4,-0.9,-1.3,-1.3,-1.2,-0.1,-0.3,-1.3,-1.0,-0.7, // AU
		-0.7,-0.2,-0.9,-0.8,+0.2,-1.0,-0.1,-1.0,-1.0,-0.6,-1.2,-0.9,-0.3,-0.9,-1.2,-1.3, // CA
		-1.1,-0.8,-1.1,-1.6,-1.0,+0.3,-0.4,-0.6,-1.2,-0.8,-0.5,-1.2,-1.5,-0.4,-1.5,-1.7, // CC
		-0.9,-0.2,-0.4,-0.4,-0.1,-0.4,+0.9,-0.7,-0.7,+0.3,-0.2,+0.1,-0.5,-0.6,-0.2,-0.9, // CG
		-1.3,-1.0,-1.2,-0.9,-1.0,-0.6,-0.7,+0.1,-1.3,-1.1,-0.9,-0.7,-1.3,-1.1,-1.2,-0.7, // CU
		-0.4,-0.8,-0.8,-1.3,-1.0,-1.2,-0.7,-1.3,+0.6,-0.8,-0.2,-0.8,-0.7,-0.8,-1.0,-0.9, // GA
		-1.2,-0.2,-0.9,-1.3,-0.6,-0.8,+0.3,-1.1,-0.8,+0.3,-0.6,-0.1,-0.7,-0.8,-0.5,-1.4, // GC
		-1.5,-1.1,-0.2,-1.2,-1.2,-0.5,-0.2,-0.9,-0.2,-0.6,+0.7,-0.5,-1.1,-0.7,-0.8,-0.9, // GG
		-1.4,-1.0,-0.9,-0.1,-0.9,-1.2,+0.1,-0.7,-0.8,-0.1,-0.5,+0.7,-0.5,-0.9,-0.3,-0.6, // GU
		-0.6,-0.5,-0.9,-0.3,-0.3,-1.5,-0.5,-1.3,-0.7,-0.7,-1.1,-0.5,+0.7,-0.5,-0.2,-0.4, // UA
		-0.9,-0.8,-1.1,-1.3,-0.9,-0.4,-0.6,-1.1,-0.8,-0.8,-0.7,-0.9,-0.5,+0.3,-0.9,-0.6, // UC
		-1.4,-1.1,-1.0,-1.0,-1.2,-1.5,-0.2,-1.2,-1.0,-0.5,-0.8,-0.3,-0.2,-0.9,+0.2,-0.9, // UG
		-0.7,-1.4,-1.2,-0.7,-1.3,-1.7,-0.9,-0.7,-0.9,-1.4,-0.9,-0.6,-0.4,-0.6,-0.9,-0.0, // UU
	};
	int N16[16] = {0, 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240};
	int *Op, *Os;
	double *Pp, *Ps;
	if (RNA==2) { // prepare for scoring pairs
		// initialize arrays of pair frequencies
		Pp = Calloc(lp*16, double); // initialized to zero
		Ps = Calloc(ls*16, double); // initialized to zero
		
		for (i = 1; i < lp; i++) {
			z = 0;
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 4; k++) {
					tot = sqrt(pnorm[i - 1]*pnorm[i]);
					if (tot > 0)
						Pp[z + i*16] = pprofile[j + size*(i - 1)]*pprofile[k + size*i]*tot;
					z++;
				}
			}
		}
		for (i = 1; i < ls; i++) {
			z = 0;
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 4; k++) {
					tot = sqrt(snorm[i - 1]*snorm[i]);
					if (tot > 0)
						Ps[z + i*16] = sprofile[j + size*(i - 1)]*sprofile[k + size*i]*tot;
					z++;
				}
			}
		}
		
		// initialize arrays of nucleotide orders
		Op = Calloc(lp*16, int); // initialized to zero
		Os = Calloc(ls*16, int); // initialized to zero
		
		for (i = 0; i < lp; i++) {
			k = 0;
			z = -1;
			for (j = 0; j < 16; j++) {
				if (Pp[j + i*16] > 0) {
					Op[k + i*16] = j;
					k++;
				} else if (z==-1) {
					z = j;
				}
			}
			if (z!=-1)
				Op[k + i*16] = z;
		}
		for (i = 0; i < ls; i++) {
			k = 0;
			z = -1;
			for (j = 0; j < 16; j++) {
				if (Ps[j + i*16] > 0) {
					Os[k + i*16] = j;
					k++;
				} else if (z==-1) {
					z = j;
				}
			}
			if (z!=-1)
				Os[k + i*16] = z;
		}
	}
	
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
		current = pprofile[7 + size*i];
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
		current = pprofile[7 + size*i];
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
		current = sprofile[7 + size*i];
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
		current = sprofile[7 + size*i];
		if (current > past) {
			sstops[i] = current - past;
			past = current;
		} else if (current < past) {
			past = current;
		}
		if (i < (ls - 1))
			sstops[i] += sstops[i + 1];
	}
	
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
	for (k = 1; k < l; k++) {
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
			
			int SIZEI = size*i;
			int SIZEJ = size*j;
			
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
					gp = GE*(1 - 0.3*sprofile[4 + SIZEJ])*zip[*(o + i*ls + j - 1)] + 0.5*GO*(sprofile[6 + size*(j - 1)] - sprofile[6 + SIZEJ]);
				} else {
					gp = GO*(2 - 0.5*sprofile[5 + SIZEJ] - 0.5*sprofile[6 + SIZEJ]);
				}
				if (i > 0 && *(o + (i - 1)*ls + j) < 0) {
					gs = GE*(1 - 0.3*pprofile[4 + SIZEI])*zip[*(o + (i - 1)*ls + j)*-1] + 0.5*GO*(pprofile[6 + size*(i - 1)] - pprofile[6 + SIZEI]);
				} else {
					gs = GO*(2 - 0.5*pprofile[5 + SIZEI] - 0.5*pprofile[6 + SIZEI]);
				}
			}
			
			tot = (1 - pprofile[4 + SIZEI])*(1 - sprofile[4 + SIZEJ]); // max of dot-product
			GS = gs*tot;
			GP = gp*tot;
			
			if (pnorm[i] > 0 && snorm[j] > 0) {
				tot = sqrt(pnorm[i]*snorm[j]); // normalization factor
			} else {
				tot = 0;
			}
			
			// score aligning
			if (do_subM) {
				S = 0;
				M = 0;
				for (int n = 0; n < 4; n++) {
					if (pprofile[n + SIZEI]==0)
						continue;
					S += pprofile[n + SIZEI]*sprofile[n + SIZEJ] * *(subM + n*4 + n);
					for (int p = 0; p < 4; p++) {
						if (p==n || sprofile[p + SIZEJ]==0)
							continue;
						M += pprofile[n + SIZEI]*sprofile[p + SIZEJ] * *(subM + n*4 + p);
					}
				}
				M = S + M;
			} else { // no substitution matrix
				S = 0;
				M = 0;
				for (int n = 0; n < 4; n++) {
					if (pprofile[n + SIZEI]==0)
						continue;
					S += pprofile[n + SIZEI]*sprofile[n + SIZEJ];
				}
				M = S*PM + ((1 - pprofile[4 + SIZEI])*(1 - sprofile[4 + SIZEJ]) - S)*MM; // % mismatched = (% ungapped) - (% matched)
			}
			
			double pFreq, sFreq;
			if (RNA==2) {
				// score aligning pairs of nucleotides
				for (int n = 0; n < 16; n++) {
					pFreq = Pp[Op[n + i*16] + i*16];
					if (pFreq!=0) {
						for (int p = 0; p < 16; p++) {
							sFreq = Ps[Os[p + j*16] + j*16];
							if (sFreq!=0) {
								M += pFreq*sFreq*subD[N16[Op[n + i*16]] + Os[p + j*16]];
							} else {
								break;
							}
						}
					} else {
						break;
					}
				}
			}
			
			if (do_DBN) {
				for (int n = 0; n < d; n++) {
					for (int p = 0; p < d; p++) {
						M += pprofile[n + 8 + SIZEI]*sprofile[p + 8 + SIZEJ] * *(dbnM + n*d + p);
					}
				}
			}
			
			if (((GS + *(m + i*(ls + 1) + j + 1)) > (M + *(m + i*(ls + 1) + j))) && ((GS + *(m + i*(ls + 1) + j + 1)) > (GP + *(m + (i + 1)*(ls + 1) + j))) && sprofile[4 + SIZEJ] != 1) {
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
			} else if ((GP + *(m + (i + 1)*(ls + 1) + j)) > (M + *(m + i*(ls + 1) + j)) && pprofile[4 + SIZEI] != 1) {
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
	
	Free(pnorm);
	Free(snorm);
	
	if (RNA==2) {
		Free(Pp);
		Free(Ps);
		Free(Op);
		Free(Os);
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
		maxp = lp;
	} else {
		maxs = ls;
	}
	
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
	Free(zip);
	Free(pstarts);
	Free(pstops);
	Free(sstarts);
	Free(sstops);
	
	int *t = Calloc(l, int); // initialized to zero
	
	i = maxp - 1;
	j = maxs - 1;
	count = l - 1;
	while ((i > -1) && (j > -1)) {
		*(t + count) = *(o + i*ls + j);
		z = *(o + i*ls + j);
		//*(o + i*ls + j) = 0;
		//*(m + (i + 1)*(ls + 1) + j + 1) = 0;
		if (z <= 0)
			i--;
		if (z >= 0)
			j--;
		count--;
	}
	int minp = i + 2;
	int mins = j + 2;
	
	int *p_ins = Calloc(l, int); // initialized to zero
	int *s_ins = Calloc(l, int); // initialized to zero
	int *p_ats = Calloc(l, int); // initialized to zero
	int *s_ats = Calloc(l, int); // initialized to zero
	
	int p_count, s_count, value;
	
	if (mins < minp) {
		*(s_ins) = minp - 1;
		*(s_ats) = 1;
		p_count = 0;
		s_count = 1;
	} else if (minp < mins) {
		*(p_ins) = mins - 1;
		*(p_ats) = 1;
		p_count = 1;
		s_count = 0;
	} else {
		p_count = 0;
		s_count = 0;
	}
	
	i = minp;
	j = mins;
	k = count + 1;
	while (k < l) {
		value = *(t + k);
		for (count = k + 1; count < l; count++) {
			if (*(t + count)!=value)
				break;
		}
		z = count - k;
		if (value==0) {
			i += z;
			j += z;
		} else if (value > 0) {
			p_ins[p_count] = z;
			p_ats[p_count] = i;
			j += z;
			p_count++;
		} else { // value < 0
			s_ins[s_count] = z;
			s_ats[s_count] = j;
			i += z;
			s_count++;
		}
		k = count;
	}
	
	if (ls > maxs) {
		p_ins[p_count] = ls - maxs;
		p_ats[p_count] = lp + 1;
		p_count++;
	}
	
	if (lp > maxp) {
		s_ins[s_count] = lp - maxp;
		s_ats[s_count] = ls + 1;
		s_count++;
	}
	
	PROTECT(ans1 = allocVector(INTSXP, p_count));
	rans = INTEGER(ans1);
	for (i = 0; i < p_count; i++) {
		*(rans + i) = p_ats[i];
	}
	PROTECT(ans2 = allocVector(INTSXP, p_count));
	rans = INTEGER(ans2);
	for (i = 0; i < p_count; i++) {
		*(rans + i) = p_ins[i];
	}
	PROTECT(ans3 = allocVector(INTSXP, s_count));
	rans = INTEGER(ans3);
	for (i = 0; i < s_count; i++) {
		*(rans + i) = s_ats[i];
	}
	PROTECT(ans4 = allocVector(INTSXP, s_count));
	rans = INTEGER(ans4);
	for (i = 0; i < s_count; i++) {
		*(rans + i) = s_ins[i];
	}
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	SET_VECTOR_ELT(ret_list, 2, ans3);
	SET_VECTOR_ELT(ret_list, 3, ans4);
	
	Free(o);
	Free(t);
	Free(p_ins);
	Free(s_ins);
	Free(p_ats);
	Free(s_ats);
	
	UNPROTECT(5);
	
	return ret_list;
}

SEXP alignProfilesAA(SEXP p, SEXP s, SEXP subMatrix, SEXP hecMatrix, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads)
{
	int i, j, k, pos, start, end, *rans, count, z;
	double *pprofile, *sprofile, gp, gs, M, GP, GS, R;
	double max, tot, freq;
	SEXP ans1, ans2, ans3, ans4, dims;
	
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
	int do_HEC, d, size = 29;
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
	int l = lp + ls;
	
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
	double GC[20] = {-0.05,-0.354,0.486,0.016,-0.014,0.611,0.037,0.642,0.107,-1.229,-0.788,0.178,-1.346,-0.718,1.03,1.265,0.441,-0.916,-1.127,-0.728};
	
	// gap modulation factor based on +/- positions of gap
	// half of log-odds because added twice (gap open and close)
	// units in third-bits - may not scale with all substitution matrices
	double GM[320] = { // gap modulate position -8 to +8
		0.1,0.171,0.129,0.094,0.082,0.094,0.085,-0.285,-0.255,-0.102,-0.095,-0.065,-0.013,-0.019,0.014,0.036, // A
		0.078,0.04,0.095,-0.013,0.108,0.064,-0.052,0.073,-0.023,-0.115,-0.165,-0.129,-0.164,-0.143,-0.087,-0.088, // R
		-0.15,-0.152,-0.127,-0.103,0.045,0.165,0.208,0.234,0.384,0.249,0.104,-0.018,-0.176,-0.174,-0.158,-0.101, // N
		-0.211,-0.152,-0.132,-0.123,-0.033,0.224,0.147,0.557,0.663,0.499,0.266,0.143,0.099,0.074,-0.008,-0.023, // D
		0.059,0.025,-0.033,0.038,-0.086,-0.221,-0.374,-0.437,-0.321,-0.432,-0.178,-0.035,-0.067,0.202,0.214,0.271, // C
		0.14,0.033,0.066,0.143,0.164,0.166,0.247,0.232,0.078,0.006,0.006,-0.076,-0.021,-0.018,-0.078,-0.06, // Q
		0.057,0.116,0.157,0.149,0.186,0.345,0.387,0.401,0.386,0.245,0.204,0.214,0.068,0.04,0.027,-0.031, // E
		-0.151,-0.21,-0.203,-0.193,-0.073,0.163,0.361,0.573,0.519,0.292,0.129,-0.01,-0.121,-0.184,-0.195,-0.077, // G
		-0.022,-0.064,-0.092,-0.081,-0.162,-0.109,-0.15,-0.207,-0.396,-0.092,-0.146,-0.143,-0.162,-0.225,-0.19,-0.068, // H
		-0.003,-0.002,-0.16,-0.133,-0.35,-0.589,-0.729,-0.877,-0.635,-0.549,-0.228,-0.166,-0.02,0.037,0.03,0.066, // I
		0.053,0.04,0.042,0.014,-0.056,-0.284,-0.434,-0.549,-0.583,-0.57,-0.355,-0.268,-0.116,-0.063,-0.039,-0.067, // L
		0.02,0.108,0.154,0.207,0.286,0.333,0.428,0.499,0.308,0.27,0.199,0.091,0.06,-0.059,-0.006,-0.017, // K
		-0.228,-0.229,-0.376,-0.343,-0.454,-0.572,-0.735,-0.754,-1.024,-0.929,-0.795,-0.826,-0.522,-0.462,-0.498,-0.484, // M
		-0.014,-0.053,-0.131,-0.057,-0.227,-0.345,-0.52,-0.594,-0.632,-0.603,-0.359,-0.175,-0.159,-0.004,0.009,-0.076, // F
		-0.03,-0.014,0.092,0.149,0.277,0.237,0.48,0.562,0.311,0.644,0.574,0.391,0.343,0.298,0.217,0.146, // P
		0.039,0.019,0.088,0.08,0.216,0.267,0.315,0.24,0.318,0.392,0.23,0.243,0.156,0.14,0.119,0.169, // S
		0.029,-0.029,-0.022,-0.021,-0.007,0.016,0.13,-0.002,-0.018,0.203,0.19,0.215,0.175,0.101,0.094,0.07, // T
		0.106,0.115,-0.055,-0.11,-0.257,-0.406,-0.642,-0.494,-0.583,-0.753,-0.314,-0.146,-0.101,-0.02,0.126,0.074, // W
		-0.01,0.027,-0.028,-0.085,-0.189,-0.368,-0.545,-0.666,-0.548,-0.608,-0.408,-0.228,-0.092,-0.163,-0.103,-0.077, // Y
		0.024,0,0.008,-0.021,-0.268,-0.435,-0.577,-0.616,-0.231,-0.218,-0.085,0.068,0.151,0.211,0.197,0.124, // V
	};
	
	// gap extension based on residues opposing the gap
	double *GCp = Calloc(lp, double); // initialized to zero
	double *GCs = Calloc(ls, double); // initialized to zero
	// gap opening based on local sequence context
	double *GOp = Calloc(lp, double); // initialized to zero
	double *GOs = Calloc(ls, double); // initialized to zero
	// gap opening in the opposing sequence based on runs
	double *GRp = Calloc(lp, double); // initialized to zero
	double *GRs = Calloc(ls, double); // initialized to zero
	
	// start of run of length 2 at position zero
	//          -2          -1            0           +1
	// -0.06983115  0.25875077  -0.86657803  -3.10770430
	
	// start of run of length > 2 at position zero
	//          -2          -1            0           +1
	// 2.038597279 2.859584452  1.913399269 -3.896055295
	
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
				
				for (k = 4; k <= 11; k++) { // -4 to +4
					pos = i + k - 7;
					if (pos < 0 || pos >= lp)
						continue;
					
					GOp[pos] += freq*GM[j*16 + 15 - k];
				}
			}
		}
		
		// do not need to check bounds
		R = pprofile[27 + size*i];
		if (R > 0) { // run of length 2
//			if (i >= 2)
//				GRp[i - 2] += R*-0.0349;
//			if (i >= 1)
//				GRp[i - 1] += R*0.1294;
//			GRp[i] += R*-0.4333;
			if (i < (lp - 1))
				GRp[i + 1] += R*-1.5539;
		}
		R = pprofile[28 + size*i];
		if (R > 0) { // run of length > 2
			if (i >= 2)
				GRp[i - 2] += R*1.0193;
			if (i >= 1)
				GRp[i - 1] += R*1.4298;
			GRp[i] += R*0.9567;
			if (i < (lp - 1))
				GRp[i + 1] += R*-1.9480;
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
				
				for (k = 4; k <= 11; k++) { // -4 to +4
					pos = i + k - 7;
					if (pos < 0 || pos >= ls)
						continue;
					
					GOs[pos] += freq*GM[j*16 + 15 - k];
				}
			}
		}
		
		// do not need to check bounds
		R = sprofile[27 + size*i];
		if (R > 0) { // run of length 2
//			if (i >= 2)
//				GRs[i - 2] += R*-0.0349;
//			if (i >= 1)
//				GRs[i - 1] += R*0.1294;
//			GRs[i] += R*-0.4333;
			if (i < (ls - 1))
				GRs[i + 1] += R*-1.5539;
		}
		R = sprofile[28 + size*i];
		if (R > 0) { // run of length > 2
			if (i >= 2)
				GRs[i - 2] += R*1.0193;
			if (i >= 1)
				GRs[i - 1] += R*1.4298;
			GRs[i] += R*0.9567;
			if (i < (ls - 1))
				GRs[i + 1] += R*-1.9480;
		}
	}
	
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
	
	// create array of normalization factors (coverge^POW)
	double *pnorm = Calloc(lp, double); // initialized to zero
	double *snorm = Calloc(ls, double); // initialized to zero
	for (i = 0; i < lp; i++)
		if (pprofile[26 + size*i] > 0)
			pnorm[i] = pow(pprofile[26 + size*i], POW);
	for (i = 0; i < ls; i++)
		if (sprofile[26 + size*i] > 0)
			snorm[i] = pow(sprofile[26 + size*i], POW);
	
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
	for (k = 1; k < l; k++) {
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
					gp = (GOp[i] + GRs[j])*(2 - 0.5*(sprofile[24 + SIZEJ] + sprofile[25 + SIZEJ]));
				}
				if (i > 0 && *(o + (i - 1)*ls + j) < 0) {
					gs = GE*(1 - 0.3*pprofile[23 + SIZEI])*zip[*(o + (i - 1)*ls + j)*-1] + 0.5*GOs[j]*(pprofile[25 + size*(i - 1)] - pprofile[25 + SIZEI]);
				} else {
					gs = (GOs[j] + GRp[i])*(2 - 0.5*(pprofile[24 + SIZEI] + pprofile[25 + SIZEI]));
				}
				
				// add modulation factor for center gap region
				// amount is based on the opposing sequence
				gp += GCs[j]*(1 - sprofile[23 + SIZEJ]);
				gs += GCp[i]*(1 - pprofile[23 + SIZEI]);
			}
			
			tot = (1 - pprofile[23 + SIZEI])*(1 - sprofile[23 + SIZEJ]); // max of dot-product
			GS = gs*tot;
			GP = gp*tot;
			
			if (pnorm[i] > 0 && snorm[j] > 0) {
				tot = sqrt(pnorm[i]*snorm[j]); // normalization factor
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
						M += pprofile[n + 29 + SIZEI]*sprofile[p + 29 + SIZEJ] * *(hecM + n*d + p);
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
	
	Free(pnorm);
	Free(snorm);
	
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
		maxp = lp;
	} else {
		maxs = ls;
	}
	
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
	Free(Op);
	Free(Os);
	Free(zip);
	Free(GCp);
	Free(GCs);
	Free(GOp);
	Free(GOs);
	Free(GRp);
	Free(GRs);
	Free(pstarts);
	Free(pstops);
	Free(sstarts);
	Free(sstops);
	
	int *t = Calloc(l, int); // initialized to zero
	
	i = maxp - 1;
	j = maxs - 1;
	count = l - 1;
	while ((i > -1) && (j > -1)) {
		*(t + count) = *(o + i*ls + j);
		z = *(o + i*ls + j);
		//*(o + i*ls + j) = 0;
		//*(m + (i + 1)*(ls + 1) + j + 1) = 0;
		if (z <= 0)
			i--;
		if (z >= 0)
			j--;
		count--;
	}
	int minp = i + 2;
	int mins = j + 2;
	
	int *p_ins = Calloc(l, int); // initialized to zero
	int *s_ins = Calloc(l, int); // initialized to zero
	int *p_ats = Calloc(l, int); // initialized to zero
	int *s_ats = Calloc(l, int); // initialized to zero
	
	int p_count, s_count, value;
	
	if (mins < minp) {
		*(s_ins) = minp - 1;
		*(s_ats) = 1;
		p_count = 0;
		s_count = 1;
	} else if (minp < mins) {
		*(p_ins) = mins - 1;
		*(p_ats) = 1;
		p_count = 1;
		s_count = 0;
	} else {
		p_count = 0;
		s_count = 0;
	}
	
	i = minp;
	j = mins;
	k = count + 1;
	while (k < l) {
		value = *(t + k);
		for (count = k + 1; count < l; count++) {
			if (*(t + count)!=value)
				break;
		}
		z = count - k;
		if (value==0) {
			i += z;
			j += z;
		} else if (value > 0) {
			p_ins[p_count] = z;
			p_ats[p_count] = i;
			j += z;
			p_count++;
		} else { // value < 0
			s_ins[s_count] = z;
			s_ats[s_count] = j;
			i += z;
			s_count++;
		}
		k = count;
	}
	
	if (ls > maxs) {
		p_ins[p_count] = ls - maxs;
		p_ats[p_count] = lp + 1;
		p_count++;
	}
	
	if (lp > maxp) {
		s_ins[s_count] = lp - maxp;
		s_ats[s_count] = ls + 1;
		s_count++;
	}
	
	PROTECT(ans1 = allocVector(INTSXP, p_count));
	rans = INTEGER(ans1);
	for (i = 0; i < p_count; i++) {
		*(rans + i) = p_ats[i];
	}
	PROTECT(ans2 = allocVector(INTSXP, p_count));
	rans = INTEGER(ans2);
	for (i = 0; i < p_count; i++) {
		*(rans + i) = p_ins[i];
	}
	PROTECT(ans3 = allocVector(INTSXP, s_count));
	rans = INTEGER(ans3);
	for (i = 0; i < s_count; i++) {
		*(rans + i) = s_ats[i];
	}
	PROTECT(ans4 = allocVector(INTSXP, s_count));
	rans = INTEGER(ans4);
	for (i = 0; i < s_count; i++) {
		*(rans + i) = s_ins[i];
	}
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	SET_VECTOR_ELT(ret_list, 2, ans3);
	SET_VECTOR_ELT(ret_list, 3, ans4);
	
	Free(o);
	Free(t);
	Free(p_ins);
	Free(s_ins);
	Free(p_ats);
	Free(s_ats);
	
	UNPROTECT(5);
	
	return ret_list;
}
