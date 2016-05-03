/****************************************************************************
 *                           Calculates FISH dG1                            *
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

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("calculateFISH", probes, targets, PACKAGE="DECIPHER")
SEXP calculateFISH(SEXP probes, SEXP targets)
{
	double dH_DR[4][4] = {
		-11.5, -7.8, -7, -8.3,
		-10.4, -12.8, -16.3, -9.1,
		-8.6, -8, -9.3, -5.9,
		-7.8, -5.5, -9, -7.8
	};
	double dS_DR[4][4] = {
		-36.4, -21.6, -19.7, -23.9,
		-28.4, -31.9, -47.1, -23.5,
		-22.9, -17.1, -23.2, -12.3,
		-23.2, -13.5, -26.1, -21.9
	};
	double dH_DD[4][4] = {
		-7.9, -8.4, -7.8, -7.2,
		-8.5, -8, -10.6, -7.8,
		-8.2, -9.8, -8, -8.4,
		-7.2, -8.2, -8.5, -7.9
	};
	double dS_DD[4][4] = {
		-22.2, -22.4, -21, -20.4,
		-22.7, -19.9, -27.2, -21,
		-22.2, -24.4, -19.9, -22.4,
		-21.3, -22.2, -22.7, -22.2
	};
	double dH_RR[4][4] = {
		-6.6, -10.17, -7.65, -5.76,
		-10.56, -12.21, -7.95, -7.65,
		-13.37, -14.21, -12.21, -10.17,
		-8.11, -13.37, -10.56, -6.6
	};
	double dS_RR[4][4] = {
		-18.38, -26.03, -19.18, -15.67,
		-28.25, -30.02, -19.18, -19.18,
		-35.68, -34.85, -30.02, -26.03,
		-22.59, -35.68, -28.25, -18.38
	};
	
	int n = length(probes);
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, n, 8));
	double *rans = REAL(ans);
	
	for (int i = 0; i < n; i++) {
		// initialize variables
		rans[i] = 1.9; // dH_PM_DR_ini
		rans[n + i] = -3.9; // dS_PM_DR_ini
		rans[2*n + i] = 0; // dH_MM_DR_ini
		rans[3*n + i] = 0; // dS_MM_DR_ini
		rans[4*n + i] = 0; // dH_MM_DD_ini
		rans[5*n + i] = 0; // dS_MM_DD_ini
		rans[6*n + i] = 0; // dH_MM_RR_ini
		rans[7*n + i] = 0; // dS_MM_RR_ini
		
		int l = -1; // last probe position (unset)
		int mm = 0;
		int p;
		
		const char *probe = CHAR(STRING_ELT(probes, i));
		const char *target = CHAR(STRING_ELT(targets, i));
		
		for (int j = 0; j < length(STRING_ELT(probes, i)); j++) {
			switch (probe[j]) {
				case 'A':
				case 'a':
					p = 0;
					break;
				case 'C':
				case 'c':
					p = 1;
					break;
				case 'G':
				case 'g':
					p = 2;
					break;
				case 'T':
				case 't':
					p = 3;
					break;
				case '-':
					p = 4;
					mm = 1;
					continue;
				default:
					error("Non-DNA character found!");
					break;
			}
			
			if (l > -1) {
				rans[i] += dH_DR[l][p];
				rans[n + i] += dS_DR[l][p];
				if (mm || (probe[j] != target[j])) {
					rans[2*n + i] += dH_DR[l][p];
					rans[3*n + i] += dS_DR[l][p];
					rans[4*n + i] += dH_DD[l][p];
					rans[5*n + i] += dS_DD[l][p];
					rans[6*n + i] += dH_RR[l][p];
					rans[7*n + i] += dS_RR[l][p];
				}
			}
			
			if (probe[j] != target[j]) {
				mm = 1;
			} else {
				mm = 0;
			}
			
			l = p;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}
