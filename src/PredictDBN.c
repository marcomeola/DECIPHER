/****************************************************************************
 *                Predicts 3 State Secondary RNA Structure                  *
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

// for math functions
#include <math.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

static int firstpos(const Chars_holder *P)
{
	int i;
	const char *p;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
		 i < P->length;
		 i++, p++)
	{
		if (!((*p) & 0x10 || (*p) & 0x40)) {
			return i;
		}
	}
	return i;
}

static int lastpos(const Chars_holder *P)
{
	int i;
	const char *p;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
		 i >= 0;
		 i--, p--)
	{
		if (!((*p) & 0x10 || (*p) & 0x40)) {
			return i;
		}
	}
	return i;
}

void Traceback(double *MI, int tot, int *unpaired, int *pos, char *states, char leftSymbol, char rightSymbol, int i, int j) {
	while (j > (i + 1)) { // prevent < 2 base hairpins
//		Rprintf("\ni = %d j = %d MI[%d, %d] = %1.0f", i + 1, j + 1, unpaired[i] + 1, unpaired[j] + 1, MI[unpaired[j]*tot + unpaired[i]]);
		if (MI[unpaired[j]*tot + unpaired[i]] > 1e9) { // bifurcation
			Traceback(MI, tot, unpaired, pos, states, leftSymbol, rightSymbol, MI[unpaired[j]*tot + unpaired[i]] - 1e9 + 1, j);
			Traceback(MI, tot, unpaired, pos, states, leftSymbol, rightSymbol, i, MI[unpaired[j]*tot + unpaired[i]] - 1e9);
			break;
		} else if (MI[unpaired[j]*tot + unpaired[i]] < 0 && MI[unpaired[j]*tot + unpaired[i]] > -1e9) {
			i -= MI[unpaired[j]*tot + unpaired[i]];
		} else if (MI[unpaired[j]*tot + unpaired[i]] < -1e9) {
			j += MI[unpaired[j]*tot + unpaired[i]] + 1e9;
		} else { // base pairing
//			if (states[pos[unpaired[i]]]!='.' || states[pos[unpaired[j]]]!='.')
//				error("crossed-over twice");
			states[pos[unpaired[i]]] = leftSymbol;
			states[pos[unpaired[j]]] = rightSymbol;
			i++;
			j--;
		}
	}
}

SEXP predictDBN(SEXP x, SEXP output, SEXP minOccupancy, SEXP impact, SEXP avgProdCorr, SEXP slope, SEXP shift, SEXP weights, SEXP pseudoknots, SEXP threshold, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	int i, j, k, p, s, d, l;
	
	// type of output
	// 1 = DBN
	// 2 = pairs
	// 3 = scores
	// 4 = structures
	int o = asInteger(output);
	double minO = asReal(minOccupancy);
	double *coef = REAL(impact);
	double APC = asReal(avgProdCorr);
	double sl = asReal(slope);
	double sh = asReal(shift);
	double *w = REAL(weights);
	int pseudo = asInteger(pseudoknots);
	double thresh = asReal(threshold);
	int before, v, last, *rPercentComplete;
	double *rans, soFar;
	int nthreads = asInteger(nThreads);
	SEXP ans, ans_s, percentComplete, utilsPackage;
	v = asLogical(verbose);
	if (v) { // percent complete variables
		soFar = 0;
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	XStringSet_holder x_set;
	x_set = hold_XStringSet(x);
	int x_length = get_length_from_XStringSet_holder(&x_set);
	Chars_holder x_s;
	x_s = get_elt_from_XStringSet_holder(&x_set, 0);
	int width = x_s.length;
	
	// initialize an array of terminal gap lengths
	int *endpoints = Calloc(2*x_length, int); // initialized to zero
	
	for (i = 0; i < x_length; i++) {
		x_s = get_elt_from_XStringSet_holder(&x_set, i);
		endpoints[i] = firstpos(&x_s);
		endpoints[x_length + i] = lastpos(&x_s);
	}
	
	// initialize an array of letter (A, C, G, T/U, other) frequencies
	double *counts = Calloc(5*width, double); // initialized to zero
	
	for (s = 0; s < x_length; s++) {
		x_s = get_elt_from_XStringSet_holder(&x_set, s);
		// assume x_s.length is equal to width
		for (p = endpoints[s]; p <= endpoints[x_length + s]; p++) {
			switch (x_s.ptr[p]) {
				case 1: // A
					counts[5*p] += w[s];
					break;
				case 2: // C
					counts[5*p + 1] += w[s];
					break;
				case 4: // G
					counts[5*p + 2] += w[s];
					break;
				case 8: // T/U
					counts[5*p + 3] += w[s];
					break;
				default: // other
					counts[5*p + 4] += w[s];
					break;
			}
		}
	}
	
	// initialize an array of positions
	int *pos = Calloc(width, int); // initialized to zero
	double sum; // sum of weights
	int tot = 0; // total number of positions >= minO
	for (p = 0; p < width; p++) {
		sum = counts[5*p] + counts[5*p + 1] + counts[5*p + 2] + counts[5*p + 3];
		if ((sum/x_length) >= minO) { // at least minOccupancy
			// normalize to sum of weights
			sum += counts[5*p + 4];
			counts[5*p] /= sum;
			counts[5*p + 1] /= sum;
			counts[5*p + 2] /= sum;
			counts[5*p + 3] /= sum;
			counts[5*p + 4] /= sum;
			
			pos[tot] = p;
			tot++;
			//Rprintf("\np = %d A=%1.2f C=%1.2f G=%1.2f U=%1.2f other=%1.2f", p + 1, counts[5*p], counts[5*p + 1], counts[5*p + 2], counts[5*p + 3], counts[5*p + 4]);
		}
	}
	
	// initialize an array of mutual information
	double *MI = Calloc(tot*tot, double); // initialized to zero
	double *rowMeans = Calloc(tot, double); // initialized to zero
	
	last = tot - 1;
	for (i = 0; i < (tot - 1); i++) {
		#pragma omp parallel for private(j,l,s,x_s) schedule(guided) num_threads(nthreads)
		for (j = i + 1; j < tot; j++) {
			double AU = 0, UA = 0, GC = 0, CG = 0, GU = 0, UG = 0, other = 0, temp = 0, bg;
			l = 0;
			
			for (s = 0; s < x_length; s++) {
				if (i < endpoints[s] || j > endpoints[x_length + s])
					continue;
				
				x_s = get_elt_from_XStringSet_holder(&x_set, s);
				
				int p1, p2;
				switch (x_s.ptr[pos[i]]) {
					case 1: // A
						p1 = 1;
						break;
					case 2: // C
						p1 = 2;
						break;
					case 4: // G
						p1 = 3;
						break;
					case 8: // T/U
						p1 = 4;
						break;
					default: // other
						p1 = 5;
						break;
				}
				switch (x_s.ptr[pos[j]]) {
					case 1: // A
						p2 = 1;
						break;
					case 2: // C
						p2 = 2;
						break;
					case 4: // G
						p2 = 3;
						break;
					case 8: // T/U
						p2 = 4;
						break;
					default: // other
						p2 = 5;
						break;
				}
				
				if (p1==5) { // -
					other += w[s];
				} else if (p2==5) { // -
					other += w[s];
				} else if (p1==1) { // A
					if (p2==4) { // T/U
						AU += w[s];
					} else {
						other += w[s];
					}
				} else if (p1==2) { // C
					if (p2==3) { // G
						CG += w[s];
					} else {
						other += w[s];
					}
				} else if (p1==3) { // G
					if (p2==2) { // C
						GC += w[s];
					} else if (p2==4) { // T/U
						GU += w[s];
					} else {
						other += w[s];
					}
				} else if (p1==4) { // T/U
					if (p2==1) { // A
						UA += w[s];
					} else if (p2==3) { // G
						UG += w[s];
					} else {
						other += w[s];
					}
				}
			}
			
			// normalize to x_length
			AU /= x_length;
			UA /= x_length;
			GC /= x_length;
			CG /= x_length;
			GU /= x_length;
			UG /= x_length;
			other /= x_length; // other does not include terminal gaps
			
			bg = counts[5*pos[i]]*counts[5*pos[j] + 3]; // AU
			if (bg > 0 && AU > 0)
				temp += AU*log2(AU/bg)*coef[0];
			bg = counts[5*pos[i] + 3]*counts[5*pos[j]]; // UA
			if (bg > 0 && UA > 0)
				temp += UA*log2(UA/bg)*coef[0];
			bg = counts[5*pos[i] + 2]*counts[5*pos[j] + 1]; // GC
			if (bg > 0 && GC > 0)
				temp += GC*log2(GC/bg)*coef[1];
			bg = counts[5*pos[i] + 1]*counts[5*pos[j] + 2]; // CG
			if (bg > 0 && CG > 0)
				temp += CG*log2(CG/bg)*coef[1];
			bg = counts[5*pos[i] + 2]*counts[5*pos[j] + 3]; // GU
			if (bg > 0 && GU > 0)
				temp += GU*log2(GU/bg)*coef[2];
			bg = counts[5*pos[i] + 3]*counts[5*pos[j] + 2]; // UG
			if (bg > 0 && UG > 0)
				temp += UG*log2(UG/bg)*coef[2];
			
			temp += other*coef[3];
			MI[i*tot + j] = temp;
			//Rprintf("\ni = %d j = %d MI = %1.2f", pos[i] + 1, pos[j] + 1, MI[i*tot + j]);
			
			rowMeans[j] += temp;
		}
		for (j = i + 1; j < tot; j++)
			rowMeans[i] += MI[i*tot + j];
		
		if (v) { // print the percent completed so far
			soFar = (2*last - i)*(i + 1);
			*rPercentComplete = floor(100*soFar/(last*(last + 1)));
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	Free(endpoints);
	Free(counts);
	
	// apply Average Product Correction (APC)
	double avg;
	if (APC != 0) {
		avg = 0;
		for (i = 0; i < tot; i++) {
			rowMeans[i] /= tot - 1;
			avg += rowMeans[i];
		}
		avg /= tot;
		
		for (i = 0; i < (tot - 1); i++)
			for (j = i + 1; j < tot; j++)
				MI[i*tot + j] -= APC*rowMeans[i]*rowMeans[j]/avg;
	}
	Free(rowMeans);
	
//	for (i = 0; i < tot; i++) {
//		Rprintf("\n");
//		for (j = 0; j < tot; j++)
//			Rprintf("%1.2f ", MI[j*tot + i]);
//	}
	
	// determine the nth highest value
	int n = (int)(tot/2); // max number of MI values that could be paired
	double *vals = Calloc(n, double); // initialized to zero
	for (i = 0; i < n; i++)
		vals[i] = -1e12;
	double minVal = -1e12; // minimum value in vals
	int index = 0; // index of minimum in vals
	// determine the nth largest value in MI
	for (i = 0; i < (tot - 1); i++) {
		for (j = i + 1; j < tot; j++) {
			if (MI[i*tot + j] > minVal) {
				// replace minVal
				vals[index] = MI[i*tot + j];
				minVal = vals[index];
				
				// find the new minVal
				for (p = 0; p < n; p++) {
					if (vals[p] < minVal) {
						minVal = vals[p];
						index = p;
					}
				}
			}
		}
	}
	Free(vals);
	
	// apply sigmoidal transformation
	sh *= minVal; // shift
	sl /= -1*(sh - minVal); // -slope
	for (i = 0; i < (tot - 1); i++) {
		for (j = i + 1; j < tot; j++) {
			MI[i*tot + j] = 1/(1 + exp(sl*log(MI[i*tot + j]/sh)));
			if (ISNAN(MI[i*tot + j])) {
				MI[i*tot + j] = 0;
			} else if (o < 3 && MI[i*tot + j] < thresh) {
				MI[i*tot + j] = 0;
			}
			//Rprintf("\ni = %d j = %d MI = %1.2f", pos[i] + 1, pos[j] + 1, MI[i*tot + j]);
		}
	}
	//Rprintf("\nshift = %1.2f slope = %1.2f minVal = %1.2f", sh, -1*sl, minVal);
	
	char *states;
	int *unpaired, q, *leftMax, *rightMax;
	double *MI2, *rowMax, *colMax;
	double match, left, right, prevL, prevR;
	char leftSymbol, rightSymbol;
	
	if (o < 3) { // perform traceback
		l = tot; // number remaining unpaired
		q = 10; // block size
		states = Calloc(width + 1, char);
		for (i = 0; i < width; i++)
			states[i] = '-';
		for (i = 0; i < tot; i++)
			states[pos[i]] = '.';
		states[width] = '\0'; // end (null terminate) the string
		unpaired = Calloc(tot, int); // initialized to zero
		for (i = 1; i < tot; i++)
			unpaired[i] = i;
		for (p = 0; p <= pseudo; p++) {
			if (p < pseudo) { // copy MI
				MI2 = Calloc(tot*tot, double); // initialized to zero
				for (i = 0; i < (tot - 1); i++)
					for (j = i + 1; j < tot; j++)
						MI2[i*tot + j] = MI[i*tot + j];
			}
			
			n = ceil((double)l/(double)q);
			rowMax = Calloc(l*n, double); // initialized to zero
			colMax = Calloc(l*n, double); // initialized to zero
			for (d = 2; d <= l; d++) {
				#pragma omp parallel for private(i,j,k,match,left,right,prevL,prevR) schedule(guided) num_threads(nthreads)
				for (i = 0; i < (l - d + 1); i++) {
					j = i + d - 1; // i <= j
					
					if (d > 2) { // not along first diagonals
						match = MI[unpaired[i]*tot + unpaired[j]];
						match += MI[(unpaired[i + 1])*tot + unpaired[j - 1]];
					} else {
						match = 0;
					}
					
					left = MI[(unpaired[i + 1])*tot + unpaired[j]];
					if (i < l) {
						prevL = MI[unpaired[j]*tot + unpaired[i + 1]];
						if (prevL < 0 && prevL > -1e9) {
							prevL -= 1;
						} else {
							prevL = -1;
						}
					} else {
						prevL = -1;
					}
					
					right = MI[unpaired[i]*tot + unpaired[j - 1]];
					if (j > 1) {
						prevR = MI[(unpaired[j - 1])*tot + unpaired[i]];
						if (prevR < -1e9) {
							prevR -= 1;
						} else {
							prevR = -1e9 - 1;
						}
					} else {
						prevR = -1e9 - 1;
					}
					
					if (match > left && match > right) {
						MI[unpaired[i]*tot + unpaired[j]] = match;
					} else if (left > right) {
						MI[unpaired[i]*tot + unpaired[j]] = left;
						MI[unpaired[j]*tot + unpaired[i]] = prevL;
					} else {
						MI[unpaired[i]*tot + unpaired[j]] = right;
						MI[unpaired[j]*tot + unpaired[i]] = prevR;
					}
					
					// bifurcation
					k = i + 3;
					while (k <= (j - 4)) {
						if ((k % q)==0) {
							if ((rowMax[i*n + k/q] + colMax[j*n + k/q]) <= MI[unpaired[i]*tot + unpaired[j]]) {
								k += q;
								continue;
							}
						}
						
						if ((MI[unpaired[i]*tot + unpaired[k]] + MI[(unpaired[k + 1])*tot + unpaired[j]]) > MI[unpaired[i]*tot + unpaired[j]]) {
							MI[unpaired[i]*tot + unpaired[j]] = MI[unpaired[i]*tot + unpaired[k]] + MI[(unpaired[k + 1])*tot + unpaired[j]];
							MI[unpaired[j]*tot + unpaired[i]] = k + 1e9;
						}
						k++;
					}
					if (MI[unpaired[i]*tot + unpaired[j]] > rowMax[i*n + j/q])
						rowMax[i*n + j/q] = MI[unpaired[i]*tot + unpaired[j]];
					if (i > 0 && MI[unpaired[i]*tot + unpaired[j]] > colMax[j*n + (i - 1)/q])
						colMax[j*n + (i - 1)/q] = MI[unpaired[i]*tot + unpaired[j]];
				}
			}
			Free(rowMax);
			Free(colMax);
			
			if (p==0) {
				leftSymbol = '(';
				rightSymbol = ')';
			} else if (p==1) {
				leftSymbol = '[';
				rightSymbol = ']';
			} else if (p==2) {
				leftSymbol = '{';
				rightSymbol = '}';
			} else if (p==3) {
				leftSymbol = '<';
				rightSymbol = '>';
			}
			
			Traceback(MI, tot, unpaired, pos, states, leftSymbol, rightSymbol, 0, l - 1);
			
			if (p < pseudo) { // replace MI with the original
				Free(MI);
				MI = MI2;
				
				// excluded paired positions from unpaired
				l = 0;
				for (i = 0; i < tot; i++) {
					if (states[pos[i]]=='.')
						unpaired[l++] = i;
				}
				if (l==0) // all positions paired
					break;
			}
		}
		Free(unpaired);
		
		PROTECT(ans = allocVector(STRSXP, 1));
		SET_STRING_ELT(ans, 0, mkChar(states));
		Free(states);
	} else if (o==3) { // scores
		PROTECT(ans = allocMatrix(REALSXP, 3, width)); // [state][pos]
		rans = REAL(ans);
		for (i = 0; i < 3*width; i++)
			rans[i] = 0;
		
		// determine the largest value in each row/column
		for (i = 0; i < (tot - 1); i++) {
			for (j = i + 1; j < tot; j++) {
				if (MI[i*tot + j] > rans[3*pos[i] + 1])
					rans[3*pos[i] + 1] = MI[i*tot + j];
				if (MI[i*tot + j] > rans[3*pos[j] + 2])
					rans[3*pos[j] + 2] = MI[i*tot + j];
			}
		}
		
		// normalize the scores
		for (i = 0; i < tot; i++) {
			sum = rans[3*pos[i] + 1] + rans[3*pos[i] + 2];
			if (sum > 1) {
				rans[3*pos[i] + 1] /= sum;
				rans[3*pos[i] + 2] /= sum;
			} else {
				rans[3*pos[i]] = 1 - sum;
			}
		}
	} else { // structures
		leftMax = Calloc(tot, int); // initialized to zero
		rightMax = Calloc(tot, int); // initialized to zero
		
		// locate the largest value in each row/column
		for (i = 0; i < (tot - 1); i++) {
			for (j = i + 1; j < tot; j++) {
				if (MI[i*tot + j] > MI[i*tot + leftMax[i]]) // same (; diff )
					leftMax[i] = j; // index = (; value = )
				if (MI[i*tot + j] > MI[rightMax[j]*tot + j]) // same ); diff (
					rightMax[j] = i; // index = ); value = (
			}
		}
		
		PROTECT(ans = allocVector(VECSXP, x_length));
		
		for (s = 0; s < x_length; s++) {
			x_s = get_elt_from_XStringSet_holder(&x_set, s);
			
			n = 0;
			int *nucs = Calloc(width, int); // initialized to zero
			
			for (i = 0; i < width; i++) {
				if (x_s.ptr[i] != 16 && x_s.ptr[i] != 32 && x_s.ptr[i] != 64)
					nucs[i] = n++;
			}
			
			PROTECT(ans_s = allocMatrix(REALSXP, 3, n)); // [state][pos]
			
			rans = REAL(ans_s);
			for (i = 0; i < 3*n; i++)
				rans[i] = 0;
			
			int *anchor = Calloc(tot, int); // initialized to zero
			
			for (i = 0; i < tot; i++) {
				if (((x_s.ptr[pos[i]] & 0x1) && (x_s.ptr[pos[leftMax[i]]] & 0x8)) || // A/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[leftMax[i]]] & 0x1)) || // U/A
					((x_s.ptr[pos[i]] & 0x2) && (x_s.ptr[pos[leftMax[i]]] & 0x4)) || // C/G
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[leftMax[i]]] & 0x2)) || // G/C
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[leftMax[i]]] & 0x8)) || // G/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[leftMax[i]]] & 0x4))) { // U/G
					if (rans[3*nucs[pos[i]] + 1] < MI[i*tot + leftMax[i]])
						rans[3*nucs[pos[i]] + 1] = MI[i*tot + leftMax[i]];
					if (MI[i*tot + leftMax[i]] >= thresh)
						anchor[i] = 1; // left anchor
				} else if (MI[i*tot + leftMax[i]] >= thresh) {
					anchor[i] = -1; // missing left anchor
				}
				if (((x_s.ptr[pos[i]] & 0x1) && (x_s.ptr[pos[rightMax[i]]] & 0x8)) || // A/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[rightMax[i]]] & 0x1)) || // U/A
					((x_s.ptr[pos[i]] & 0x2) && (x_s.ptr[pos[rightMax[i]]] & 0x4)) || // C/G
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[rightMax[i]]] & 0x2)) || // G/C
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[rightMax[i]]] & 0x8)) || // G/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[rightMax[i]]] & 0x4))) { // U/G
					if (rans[3*nucs[pos[i]] + 2] < MI[rightMax[i]*tot + i])
						rans[3*nucs[pos[i]] + 2] = MI[rightMax[i]*tot + i];
					if (MI[rightMax[i]*tot + i] >= thresh)
						anchor[i] = 2; // right anchor
				} else if (MI[rightMax[i]*tot + i] >= thresh) {
					anchor[i] = -2; // missing right anchor
				}
			}
			
			// eliminate single anchors between two unanchored positions
			int last1[2] = {-1};
			int last2[2] = {-1};
			for (i = 0; i < tot; i++) {
				if (anchor[i] > 0) { // anchored
					// shift previous
					last2[0] = last1[0];
					last2[1] = last1[1];
					last1[0] = i;
					last1[1] = anchor[i];
				} else if (anchor[i] < 0) { // unanchored
					if (last1[1] > 0 && last2[1] < 0) { // middle anchor
						last1[1] *= -1; // unanchor
						anchor[last1[0]] = last1[1];
					}
					// shift previous
					last2[0] = last1[0];
					last2[1] = last1[1];
					last1[0] = i;
					last1[1] = anchor[i];
				}
			}
			
			// initialize an array of shifted positions
			int *pos2 = Calloc(width, int); // initialized to zero
			for (i = 0; i < width; i++) {
				pos2[i] = pos[i];
			}
			/*
			int previous = -1;
			int missing = 0;
			int matches, delta;
			for (i = 0; i < tot; i++) {
				if (anchor[i] > 0) { // anchored
					if (missing > 0) { // search for missing anchors
						// shift-in nucleotides up to missing-positions away
						for (delta = 1; delta <= missing; delta++) {
							matches = 0;
							for (j = previous + 1; j < i; j++) {
								// first check not matching in shifted nucs,
								// then shift-in by delta-many nucleotides,
								// now if it matches then +1 & continue to next
							}
							if (matches==missing)
								break;
							matches = 0;
							for (j = i - 1; j > previous; j--) {
								
							}
							if (matches==missing)
								break;
							matches = 0;
							// repeat for -delta
						}
					}
					
					previous = i;
					missing = 0;
				} else if (anchor[i] < 0) { // unanchored
					missing++;
				}
			}
			if (s==92) {
				for (i = 0; i < tot; i++) {
					if (anchor[i] != 0)
						Rprintf("\nnuc = %d anchor = %d", nucs[pos[i]] + 1, anchor[i]);
				}
			}
			*/
			Free(pos2);
			Free(anchor);
			
			// normalize the scores
			for (i = 0; i < tot; i++) {
				sum = rans[3*nucs[pos[i]] + 1] + rans[3*nucs[pos[i]] + 2];
				if (sum > 1) {
					rans[3*nucs[pos[i]] + 1] /= sum;
					rans[3*nucs[pos[i]] + 2] /= sum;
				} else {
					rans[3*nucs[pos[i]]] = 1 - sum;
				}
			}
			
			Free(nucs);
			
			SET_VECTOR_ELT(ans, s, ans_s);
			UNPROTECT(1); // ans_s
		}
		
		Free(leftMax);
		Free(rightMax);
	}
	
	Free(MI);
	Free(pos);
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}
