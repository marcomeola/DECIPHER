/****************************************************************************
 *                Calculate 3' Terminal Mismatch Penalties                  *
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

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("terminalMismatch", probes, targets, cut, mGaps, PACKAGE="DECIPHER")
SEXP terminalMismatch(SEXP p, SEXP t, SEXP cutoff, SEXP mGaps)
{
	int i, j, lp, lt, l, mm, count, end, gaps;
	int n = length(p);
	double *cut;
	const char *probe;
	const char *target;
	cut = REAL(cutoff);
	int maxGaps = asInteger(mGaps);
	
	double *rans;
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, n));
	rans = REAL(ans);
	
	#pragma omp parallel for private(i,j,lp,lt,l,mm,count,end,gaps,probe,target) schedule(guided)
	for (i = 0; i < n; i++) {
		lp = length(STRING_ELT(p, i)) - 1;
		lt = length(STRING_ELT(t, i)) - 1;
		if (lp > lt) {
			l = lt;
		} else {
			l = lp;
		}
		
		if (l <= 0) {
			rans[i] = 0;
			continue;
		}
		
		probe = CHAR(STRING_ELT(p, i));
		target = CHAR(STRING_ELT(t, i));
		
		mm = 0;
		count = 0;
		gaps = 0;
		end = lp;
		for (j = 0; j <= l; j++) {
			switch (probe[lp - j]) {
				case 'A':
				case 'a':
					count++;
					switch (target[lt - lp + j]) {
						case 'T':
						case 't':
							break;
						case '-':
							gaps++;
							mm++;
							break;
						default:
							mm++;
							break;
					}
					break;
				case 'C':
				case 'c':
					count++;
					switch (target[lt - lp + j]) {
						case 'G':
						case 'g':
							break;
						case '-':
							gaps++;
							mm++;
							break;
						default:
							mm++;
							break;
					}
					break;
				case 'G':
				case 'g':
					count++;
					switch (target[lt - lp + j]) {
						case 'C':
						case 'c':
							break;
						case '-':
							gaps++;
							mm++;
							break;
						default:
							mm++;
							break;
					}
					break;
				case 'T':
				case 't':
					count++;
					switch (target[lt - lp + j]) {
						case 'A':
						case 'a':
							break;
						case '-':
							gaps++;
							mm++;
							break;
						default:
							mm++;
							break;
					}
					break;
				case '-':
					if (end == (lp - j)) { // haven't reached primer 3' end yet
						end--;
					} else { // internal to primer
						switch (target[lt - lp + j]) {
							case '-':
								break;
							default:
								gaps++;
								mm++;
								count++;
								break;
						}
					}
					break;
				default:
					break;
			}
			//Rprintf("\nj=%d,count=%d,cutoff=%d,p/t=%c/%c,gaps=%d,mismatches=%d",j,count,(int)ceil(*cut * count),probe[lp - j],target[lt - lp + j],gaps,mm);
			if (gaps > maxGaps)
				break;
			if (mm > ceil(*cut * count)) // rolling distance
				break;
		}
		
		if (gaps > maxGaps) {
			rans[i] = 0;
			continue;
		}
		if (mm > ceil(*cut * count)) { // greater than max distance
			rans[i] = 0;
			continue;
		} else {
			rans[i] = 1;
		}
		
		if (lp > lt) { // 3' ends do not align
			continue;
		}
		
		// penultimate base
		switch (probe[end - 1]) {
			case 'A':
			case 'a':
				// ultimate base
				switch (probe[end]) {
					case 'A':
					case 'a':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .048;
								break;
							case 'C':
							case 'c':
								rans[i] = .099;
								break;
							case 'G':
							case 'g':
								rans[i] = .001;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'C':
					case 'c':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .575;
								break;
							case 'C':
							case 'c':
								rans[i] = .000;
								break;
							case 'T':
							case 't':
								rans[i] = .229;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'G':
					case 'g':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .010;
								break;
							case 'G':
							case 'g':
								rans[i] = .007;
								break;
							case 'T':
							case 't':
								rans[i] = .151;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'T':
					case 't':
						switch (target[lt - end]) {
							case 'C':
							case 'c':
								rans[i] = .087;
								break;
							case 'G':
							case 'g':
								rans[i] = .494;
								break;
							case 'T':
							case 't':
								rans[i] = .029;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						break;
				}
				break;
			case 'C':
			case 'c':
				// ultimate base
				switch (probe[end]) {
					case 'A':
					case 'a':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .004;
								break;
							case 'C':
							case 'c':
								rans[i] = 1;
								break;
							case 'G':
							case 'g':
								rans[i] = .003;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'C':
					case 'c':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .092;
								break;
							case 'C':
							case 'c':
								rans[i] = .033;
								break;
							case 'T':
							case 't':
								rans[i] = .346;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'G':
					case 'g':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .001;
								break;
							case 'G':
							case 'g':
								rans[i] = .088;
								break;
							case 'T':
							case 't':
								rans[i] = .278;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'T':
					case 't':
						switch (target[lt - end]) {
							case 'C':
							case 'c':
								rans[i] = .927;
								break;
							case 'G':
							case 'g':
								rans[i] = .677;
								break;
							case 'T':
							case 't':
								rans[i] = .016;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						break;
				}
				break;
			case 'G':
			case 'g':
				// ultimate base
				switch (probe[end]) {
					case 'A':
					case 'a':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .000;
								break;
							case 'C':
							case 'c':
								rans[i] = .114;
								break;
							case 'G':
							case 'g':
								rans[i] = .002;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'C':
					case 'c':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .042;
								break;
							case 'C':
							case 'c':
								rans[i] = .002;
								break;
							case 'T':
							case 't':
								rans[i] = .006;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'G':
					case 'g':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .000;
								break;
							case 'G':
							case 'g':
								rans[i] = .023;
								break;
							case 'T':
							case 't':
								rans[i] = .022;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'T':
					case 't':
						switch (target[lt - end]) {
							case 'C':
							case 'c':
								rans[i] = .004;
								break;
							case 'G':
							case 'g':
								rans[i] = .332;
								break;
							case 'T':
							case 't':
								rans[i] = .001;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						break;
				}
				break;
			case 'T':
			case 't':
				// ultimate base
				switch (probe[end]) {
					case 'A':
					case 'a':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .001;
								break;
							case 'C':
							case 'c':
								rans[i] = .035;
								break;
							case 'G':
							case 'g':
								rans[i] = .002;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'C':
					case 'c':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .029;
								break;
							case 'C':
							case 'c':
								rans[i] = .001;
								break;
							case 'T':
							case 't':
								rans[i] = .002;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'G':
					case 'g':
						switch (target[lt - end]) {
							case 'A':
							case 'a':
								rans[i] = .003;
								break;
							case 'G':
							case 'g':
								rans[i] = .001;
								break;
							case 'T':
							case 't':
								rans[i] = .019;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					case 'T':
					case 't':
						switch (target[lt - end]) {
							case 'C':
							case 'c':
								rans[i] = .048;
								break;
							case 'G':
							case 'g':
								rans[i] = .037;
								break;
							case 'T':
							case 't':
								rans[i] = .000;
								break;
							case '-':
								rans[i] = 0;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						break;
				}
				break;
			default:
				rans[i] = 0;
				break;
		}
		
		if (l == 1)
			continue;
		
		// penultimate base
		switch (probe[end - 1]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 1]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .210;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .469;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .201;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 1]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .347;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .207;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .317;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 1]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .203;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .224;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .294;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 1]) {
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .413;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .508;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .209;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			default:
				rans[i] = 0;
				break;
		}
		
		if (l == 2)
			continue;
		
		// antepenultimate base
		switch (probe[end - 2]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 2]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .408;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .602;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .401;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 2]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .511;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .405;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .487;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 2]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .402;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .418;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .470;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 2]) {
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .560;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .631;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .407;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			default:
				rans[i] = 0;
				break;
		}
		
		if (l == 3)
			continue;
		
		// preantepenultimate base
		switch (probe[end - 3]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 3]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .605;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .734;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .601;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 3]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .674;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .604;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .658;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 3]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .601;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .612;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .647;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 3]) {
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .706;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .754;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .605;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			default:
				rans[i] = 0;
				break;
		}
		
		if (l == 4)
			continue;
		
		// 5th base
		switch (probe[end - 4]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 4]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .803;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .867;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .800;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 4]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .837;
						}
						break;
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .802;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .829;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 4]) {
					case 'A':
					case 'a':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .801;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .806;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .823;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 4]) {
					case 'C':
					case 'c':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .853;
						}
						break;
					case 'G':
					case 'g':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .877;
						}
						break;
					case 'T':
					case 't':
						if (rans[i] < 1) {
							rans[i] = 0;
						} else {
							rans[i] = .802;
						}
						break;
					case '-':
						rans[i] = 0;
						break;
					default:
						break;
				}
				break;
			default:
				rans[i] = 0;
				break;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}
