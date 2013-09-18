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
SEXP terminalMismatch(SEXP p, SEXP t, SEXP cutoff, SEXP mGaps, SEXP nThreads)
{
	int i, j, lp, lt, l, mm, count, end, gaps;
	int n = length(p);
	double *cut;
	const char *probe;
	const char *target;
	cut = REAL(cutoff);
	int maxGaps = asInteger(mGaps);
	int nthreads = asInteger(nThreads);
	
	double *rans;
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, n));
	rans = REAL(ans);
	
	#pragma omp parallel for private(i,j,lp,lt,l,mm,count,end,gaps,probe,target) schedule(guided) num_threads(nthreads)
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
		
		mm = 0;
		
		// penultimate base of primer (nearest neighbor)
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = .099;
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .001;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = .000;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .229;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .007;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .151;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .494;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .029;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						mm++;
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = 1;
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .003;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = .033;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .346;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .088;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .278;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .677;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .016;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						mm++;
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = .114;
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .002;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = .002;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .006;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .023;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .022;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .332;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .001;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						mm++;
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = .035;
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .002;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'C':
							case 'c':
								rans[i] = .001;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .002;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .001;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .019;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
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
								mm++;
								break;
							case 'G':
							case 'g':
								rans[i] = .037;
								mm++;
								break;
							case 'T':
							case 't':
								rans[i] = .000;
								mm++;
								break;
							case '-':
								rans[i] = 0;
								mm++;
								break;
							default:
								break;
						}
						break;
					default:
						rans[i] = 0;
						mm++;
						break;
				}
				break;
			default:
				rans[i] = 0;
				mm++;
				break;
		}
		
		if (l == 1 || rans[i] == 0)
			continue;
		
		// 2nd base = penultimate base
		switch (probe[end - 1]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 1]) {
					case 'T':
					case 't':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .372;
							mm++;
						}
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 1]) {
					case 'G':
					case 'g':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .372;
							mm++;
						}
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 1]) {
					case 'C':
					case 'c':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .372;
							mm++;
						}
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 1]) {
					case 'A':
					case 'a':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .372;
							mm++;
						}
						break;
				}
				break;
			default:
				rans[i] = 0;
				mm++;
				break;
		}
		
		if (l == 2 || rans[i] == 0)
			continue;
		
		// 3rd base = antepenultimate base
		switch (probe[end - 2]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 2]) {
					case 'T':
					case 't':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .634;
							mm++;
						}
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 2]) {
					case 'G':
					case 'g':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .634;
							mm++;
						}
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 2]) {
					case 'C':
					case 'c':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .634;
							mm++;
						}
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 2]) {
					case 'A':
					case 'a':
						break;
					case '-':
						rans[i] = 0;
						mm++;
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .634;
							mm++;
						}
						break;
				}
				break;
			default:
				if (mm == 1) {
					rans[i] = 0;
				} else {
					rans[i] = .237;
					mm++;
				}
				break;
		}
		
		if (l == 3 || rans[i] == 0)
			continue;
		
		// 4th base = preantepenultimate base
		switch (probe[end - 3]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 3]) {
					case 'T':
					case 't':
						break;
					case '-':
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .018;
							mm++;
						}
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .308;
							mm++;
						}
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 3]) {
					case 'G':
					case 'g':
						break;
					case '-':
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .018;
							mm++;
						}
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .308;
							mm++;
						}
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 3]) {
					case 'C':
					case 'c':
						break;
					case '-':
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .018;
							mm++;
						}
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .308;
							mm++;
						}
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 3]) {
					case 'A':
					case 'a':
						break;
					case '-':
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .018;
							mm++;
						}
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							rans[i] = .308;
							mm++;
						}
						break;
				}
				break;
			default:
				if (mm == 1) {
					rans[i] = 0;
				} else {
					rans[i] = .43;
					mm++;
				}
				break;
		}
		
		if (l == 4 || rans[i] == 0)
			continue;
		
		// 5th base = ultrapreantepenultimate
		switch (probe[end - 4]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 4]) {
					case 'T':
					case 't':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 4]) {
					case 'G':
					case 'g':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 4]) {
					case 'C':
					case 'c':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 4]) {
					case 'A':
					case 'a':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			default:
				if (mm == 1) {
					rans[i] = 0;
				} else {
					mm++;
				}
				break;
		}
		
		if (l == 5 || rans[i] == 0)
			continue;
		
		// 6th base = preultrapreantepenultimate
		switch (probe[end - 5]) {
			case 'A':
			case 'a':
				switch (target[lt - end + 5]) {
					case 'T':
					case 't':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			case 'C':
			case 'c':
				switch (target[lt - end + 5]) {
					case 'G':
					case 'g':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			case 'G':
			case 'g':
				switch (target[lt - end + 5]) {
					case 'C':
					case 'c':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			case 'T':
			case 't':
				switch (target[lt - end + 5]) {
					case 'A':
					case 'a':
						break;
					default:
						if (mm == 1) {
							rans[i] = 0;
						} else {
							mm++;
						}
						break;
				}
				break;
			default:
				if (mm == 1) {
					rans[i] = 0;
				} else {
					mm++;
				}
				break;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}
