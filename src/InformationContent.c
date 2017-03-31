/****************************************************************************
 *                  Finds Information Content By Position                   *
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

// for math functions
#include <math.h>

// DECIPHER header file
#include "DECIPHER.h"

SEXP informationContent(SEXP p, SEXP nS, SEXP correction)
{
	int i, j;
	double *pprofile;
	double pNseqs = asReal(nS);
	int corr = asInteger(correction);
	
	pprofile = REAL(p);
	int size = 8;
	R_len_t lp = length(p)/size;
	
	// calculate the background character frequency distribution
	double palpha[4] = {0};
	double weight = 0;
	for (i = 0; i < lp; i++) {
		for (j = 0; j < 4; j++) {
			palpha[j] += pprofile[j + size*i]*pprofile[7 + size*i];
		}
		weight += (1 - pprofile[4 + size*i])*pprofile[7 + size*i];
	}
	if (weight > 0) { // non-gap characters present
		for (j = 0; j < 4; j++) {
			palpha[j] /= weight;
		}
	}
	
	double smallN[100] = {
		2.0000,1.4037,1.1179,0.9429,0.8223,0.7332,0.6641,0.6086,0.5630,0.5246,
		0.4919,0.4635,0.4387,0.4167,0.3972,0.3796,0.3637,0.3493,0.3361,0.3240,
		0.3128,0.3025,0.2929,0.2840,0.2757,0.2679,0.2606,0.2538,0.2473,0.2412,
		0.2354,0.2300,0.2248,0.2198,0.2151,0.2107,0.2064,0.2023,0.1984,0.1947,
		0.1911,0.1876,0.1843,0.1812,0.1781,0.1752,0.1723,0.1696,0.1669,0.1644,
		0.1619,0.1595,0.1572,0.1550,0.1528,0.1507,0.1487,0.1467,0.1448,0.1430,
		0.1411,0.1394,0.1377,0.1360,0.1344,0.1328,0.1313,0.1298,0.1283,0.1269,
		0.1255,0.1241,0.1228,0.1215,0.1203,0.1190,0.1178,0.1166,0.1155,0.1144,
		0.1132,0.1122,0.1111,0.1101,0.1091,0.1081,0.1071,0.1061,0.1052,0.1043,
		0.1034,0.1025,0.1016,0.1008,0.0999,0.0991,0.0983,0.0975,0.0967,0.0960
	};
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, lp));
	double *pentropy = REAL(ans);
	
	double freqj;
	int Nseqs;
	for (i = 0; i < lp; i++) {
		pentropy[i] = 0;
		
		// calculate entropy = sum(p*log2(p/q))
		weight = 1 - pprofile[4 + size*i];
		if (weight > 0) {
			for (j = 0; j < 4; j++) {
				freqj = pprofile[j + size*i];
				if (freqj > 0) {
					freqj /= weight;
					pentropy[i] += freqj*log2(freqj/palpha[j]);
				}
			}
		}
		
		if (corr) {
			// apply correction for small N
			Nseqs = (int)(pNseqs*weight*pprofile[7 + size*i] - 1); // index in smallN
			if (Nseqs < 100) {
				if (Nseqs < 0)
					Nseqs = 0;
				pentropy[i] -= smallN[Nseqs];
				if (pentropy[i] < 0)
					pentropy[i] = 0;
			}
		}
	}
	
	UNPROTECT(1);
	return(ans);
}

SEXP informationContentAA(SEXP p, SEXP nS, SEXP correction)
{
	int i, j;
	double *pprofile;
	double pNseqs = asReal(nS);
	int corr = asInteger(correction);
	
	pprofile = REAL(p);
	int size = 29;
	R_len_t lp = length(p)/size;
	
	// calculate the background character frequency distribution
	double palpha[20] = {0};
	double weight = 0;
	for (i = 0; i < lp; i++) {
		for (j = 0; j < 20; j++) {
			palpha[j] += pprofile[j + size*i]*pprofile[26 + size*i];
		}
		weight += (1 - pprofile[23 + size*i])*pprofile[26 + size*i];
	}
	if (weight > 0) { // non-gap characters present
		for (j = 0; j < 20; j++) {
			palpha[j] /= weight;
		}
	}
	
	double smallN[200] = {
		4.3219,3.4227,2.9295,2.5989,2.3552,2.1652,2.0113,1.8833,1.7745,1.6806,
		1.5984,1.5257,1.4608,1.4024,1.3495,1.3013,1.2571,1.2165,1.1789,1.1441,
		1.1117,1.0814,1.0530,1.0264,1.0014,0.9778,0.9554,0.9343,0.9143,0.8952,
		0.8771,0.8599,0.8434,0.8276,0.8126,0.7981,0.7843,0.7710,0.7583,0.7460,
		0.7342,0.7228,0.7118,0.7013,0.6910,0.6812,0.6716,0.6624,0.6534,0.6447,
		0.6363,0.6282,0.6203,0.6126,0.6051,0.5979,0.5908,0.5839,0.5772,0.5707,
		0.5644,0.5582,0.5522,0.5463,0.5406,0.5350,0.5295,0.5242,0.5190,0.5139,
		0.5089,0.5041,0.4993,0.4946,0.4901,0.4856,0.4812,0.4770,0.4728,0.4687,
		0.4646,0.4607,0.4568,0.4530,0.4493,0.4456,0.4420,0.4385,0.4351,0.4317,
		0.4283,0.4250,0.4218,0.4187,0.4155,0.4125,0.4095,0.4065,0.4036,0.4007,
		0.3979,0.3952,0.3924,0.3897,0.3871,0.3845,0.3819,0.3794,0.3769,0.3745,
		0.3721,0.3697,0.3673,0.3650,0.3628,0.3605,0.3583,0.3561,0.3540,0.3518,
		0.3497,0.3477,0.3456,0.3436,0.3416,0.3397,0.3377,0.3358,0.3340,0.3321,
		0.3303,0.3284,0.3266,0.3249,0.3231,0.3214,0.3197,0.3180,0.3163,0.3147,
		0.3131,0.3115,0.3099,0.3083,0.3068,0.3052,0.3037,0.3022,0.3007,0.2992,
		0.2978,0.2964,0.2949,0.2935,0.2921,0.2908,0.2894,0.2881,0.2867,0.2854,
		0.2841,0.2828,0.2815,0.2803,0.2790,0.2778,0.2766,0.2753,0.2741,0.2729,
		0.2718,0.2706,0.2694,0.2683,0.2672,0.2660,0.2649,0.2638,0.2627,0.2616,
		0.2606,0.2595,0.2585,0.2574,0.2564,0.2554,0.2543,0.2533,0.2523,0.2513,
		0.2504,0.2494,0.2484,0.2475,0.2465,0.2456,0.2447,0.2437,0.2428,0.2419
	};
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, lp));
	double *pentropy = REAL(ans);
	
	double freqj;
	int Nseqs;
	for (i = 0; i < lp; i++) {
		pentropy[i] = 0;
		
		// calculate entropy = sum(p*log2(p/q))
		weight = 1 - pprofile[23 + size*i];
		if (weight > 0) {
			for (j = 0; j < 20; j++) { // omit U in 20, O in 21
				freqj = pprofile[j + size*i];
				if (freqj > 0) {
					freqj /= weight;
					pentropy[i] += freqj*log2(freqj/palpha[j]);
				}
			}
		}
		
		if (corr) {
			// apply correction for small N
			Nseqs = (int)(pNseqs*weight*pprofile[26 + size*i] - 1); // index in smallN
			if (Nseqs < 200) {
				if (Nseqs < 0)
					Nseqs = 0;
				pentropy[i] -= smallN[Nseqs];
				if (pentropy[i] < 0)
					pentropy[i] = 0;
			}
		}
	}
	
	UNPROTECT(1);
	return(ans);
}
