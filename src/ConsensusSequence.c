/****************************************************************************
 *                         Forms Consensus Sequence                         *
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

static int frontTerminalGaps(const cachedCharSeq *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->seq;
	     i < P->length;
	     i++, p++)
	{
		if ((*p) & 0x10) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGaps(const cachedCharSeq *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->seq + P->length - 1);
	     i >= 0;
	     i--, p--)
	{
		if ((*p) & 0x10) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static void alphabetFrequency(const cachedCharSeq *P, double *bits, int seqLength, int degeneracy, int ignore, int start, int end, double weight)
{
	int j;
	const char *p;
	
	for (j = start, p = (P->seq + start);
		j < (P->length - end);
		j++, p++)
	{
		if (degeneracy==1) { // include degeneracy codes
			// another base counted
			*(bits + 6*seqLength + j) += weight;
			
			// tally the bases into the encoded array
			switch (*p) {
				case 1: // A
					*(bits + 0*seqLength + j) += weight;
					break;
				case 2: // C
					*(bits + 1*seqLength + j) += weight;
					break;
				case 3: // M
					*(bits + 0*seqLength + j) += .5*weight; *(bits + 1*seqLength + j) += .5*weight; // AC
					break;
				case 4: // G
					*(bits + 2*seqLength + j) += weight;
					break;
				case 5: // R
					*(bits + 0*seqLength + j) += .5*weight; *(bits + 2*seqLength + j) += .5*weight; // AG
					break;
				case 6: // S
					*(bits + 1*seqLength + j) += .5*weight; *(bits + 2*seqLength + j) += .5*weight; // CG
					break;
				case 7: // V
					*(bits + 0*seqLength + j) += (double)1/3*weight; *(bits + 1*seqLength + j) += (double)1/3*weight; *(bits + 2*seqLength + j) += (double)1/3*weight; // ACG
					break;
				case 8: // T
					*(bits + 3*seqLength + j) += weight;
					break;
				case 9: // W
					*(bits + 0*seqLength + j) += .5*weight; *(bits + 3*seqLength + j) += .5*weight; // AT
					break;
				case 10: // Y
					*(bits + 1*seqLength + j) += .5*weight; *(bits + 3*seqLength + j) += .5*weight; // CT
					break;
				case 11: // H
					*(bits + 0*seqLength + j) += (double)1/3*weight; *(bits + 1*seqLength + j) += (double)1/3*weight; *(bits + 3*seqLength + j) += (double)1/3*weight; // ACT
					break;
				case 12: // K
					*(bits + 2*seqLength + j) += .5*weight; *(bits + 3*seqLength + j) += .5*weight; // GT
					break;
				case 13: // D
					*(bits + 0*seqLength + j) += (double)1/3*weight; *(bits + 2*seqLength + j) += (double)1/3*weight; *(bits + 3*seqLength + j) += (double)1/3*weight; // AGT
					break;
				case 14: // B
					*(bits + 1*seqLength + j) += (double)1/3*weight; *(bits + 2*seqLength + j) += (double)1/3*weight; *(bits + 3*seqLength + j) += (double)1/3*weight; // CGT
					break;
				case 15: // N
					*(bits + 0*seqLength + j) += .25*weight; *(bits + 1*seqLength + j) += .25*weight; *(bits + 2*seqLength + j) += .25*weight; *(bits + 3*seqLength + j) += .25*weight; // ACGT
					break;
				case 16: // -
					if (ignore==1) { // don't include gaps
						*(bits + 6*seqLength + j) -= weight;
					} else { // include gaps
						*(bits + 4*seqLength + j) += weight;
					}
					break;
				case 32: // +
					if (ignore==1) { // don't include masks
						*(bits + 6*seqLength + j) -= weight;
					} else { // include masks
						*(bits + 5*seqLength + j) += weight;
					}
					break;
				default:
					error("not DNA!");
					break;
			}
		} else { // don't include degeneracy codes
			switch (*p) {
				case 1: // A
					*(bits + 0*seqLength + j) += weight;
					*(bits + 6*seqLength + j) += weight;
					break;
				case 2: // C
					*(bits + 1*seqLength + j) += weight;
					*(bits + 6*seqLength + j) += weight;
					break;
				case 4: // G
					*(bits + 2*seqLength + j) += weight;
					*(bits + 6*seqLength + j) += weight;
					break;
				case 8: // T
					*(bits + 3*seqLength + j) += weight;
					*(bits + 6*seqLength + j) += weight;
					break;
				case 16: // -
					*(bits + 4*seqLength + j) += weight;
					*(bits + 6*seqLength + j) += weight;
					break;
				case 32: // +
					*(bits + 5*seqLength + j) += weight;
					*(bits + 6*seqLength + j) += weight;
					break;
			}
		}
	}
}

static void makeConsensus(double *bits, char *seq, int seqLength, int x_length, double threshold, double minInfo, int tGaps)
{
	int j;
	double percentA, percentC, percentG, percentT, percentGap, percentMask;
	double information;
	
	for (j = 0; j < seqLength; j++) {
		
		if (!tGaps && (*(bits + 6*seqLength + j)==0)) {
			*(seq + j) = '-';
			continue;
		}
		
		// find the percentage of each base in position j
		percentA = (double)(*(bits + 0*seqLength + j))/(*(bits + 6*seqLength + j));	
		percentC = (double)(*(bits + 1*seqLength + j))/(*(bits + 6*seqLength + j));	
		percentG = (double)(*(bits + 2*seqLength + j))/(*(bits + 6*seqLength + j));	
		percentT = (double)(*(bits + 3*seqLength + j))/(*(bits + 6*seqLength + j));	
		percentGap = (double)(*(bits + 4*seqLength + j))/(*(bits + 6*seqLength + j));
		percentMask = (double)(*(bits + 5*seqLength + j))/(*(bits + 6*seqLength + j));
		information = 0;
		
		// determine the most common base over the threshold percentage
		if (percentA >= threshold && percentC >= threshold && percentG >= threshold && percentT >= threshold) {
			if ((percentA + percentC + percentG + percentT) >= percentGap &&
				(percentA + percentC + percentG + percentT) >= percentMask) {
				*(seq + j) = 'N';
				information = percentA + percentC + percentG + percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentC >= threshold && percentG >= threshold && percentT >= threshold) {
			if ((percentC + percentG + percentT) >= percentGap &&
				(percentC + percentG + percentT) >= percentMask) {
				*(seq + j) = 'B';
				information = percentC + percentG + percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentA >= threshold && percentG >= threshold && percentT >= threshold) {
			if ((percentA + percentG + percentT) >= percentGap &&
				(percentA + percentG + percentT) >= percentMask) {
				*(seq + j) = 'D';
				information = percentA + percentG + percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentA >= threshold && percentC >= threshold && percentT >= threshold) {
			if ((percentA + percentC + percentT) >= percentGap &&
				(percentA + percentC + percentT) >= percentMask) {
				*(seq + j) = 'H';
				information = percentA + percentC + percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentA >= threshold && percentC >= threshold && percentG >= threshold) {
			if ((percentA + percentC + percentG) >= percentGap &&
				(percentA + percentC + percentG) >= percentMask) {
				*(seq + j) = 'V';
				information = percentA + percentC + percentG;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentC >= threshold && percentT >= threshold) {
			if ((percentC + percentT) >= percentGap &&
				(percentC + percentT) >= percentMask) {
				*(seq + j) = 'Y';
				information = percentC + percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentG >= threshold && percentT >= threshold) {
			if ((percentG + percentT) >= percentGap &&
				(percentG + percentT) >= percentMask) {
				*(seq + j) = 'K';
				information = percentG + percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentA >= threshold && percentT >= threshold) {
			if ((percentA + percentT) >= percentGap &&
				(percentA + percentT) >= percentMask) {
				*(seq + j) = 'W';
				information = percentA + percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentC >= threshold && percentG >= threshold) {
			if ((percentC + percentG) >= percentGap &&
				(percentC + percentG) >= percentMask) {
				*(seq + j) = 'S';
				information = percentC + percentG;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentA >= threshold && percentG >= threshold) {
			if ((percentA + percentG) >= percentGap &&
				(percentA + percentG) >= percentMask) {
				*(seq + j) = 'R';
				information = percentA + percentG;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentA >= threshold && percentC >= threshold) {
			if ((percentA + percentC) >= percentGap &&
				(percentA + percentC) >= percentMask) {
				*(seq + j) = 'M';
				information = percentA + percentC;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentA >= threshold) {
			if (percentA >= percentGap &&
				percentA >= percentMask) {
				*(seq + j) = 'A';
				information = percentA;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentC >= threshold) {
			if (percentC >= percentGap &&
				percentC >= percentMask) {
				*(seq + j) = 'C';
				information = percentC;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentG >= threshold) {
			if (percentG >= percentGap &&
				percentG >= percentMask) {
				*(seq + j) = 'G';
				information = percentG;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentT >= threshold) {
			if (percentT >= percentGap &&
				percentT >= percentMask) {
				*(seq + j) = 'T';
				information = percentT;
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (percentGap >= threshold ||
				   percentMask >= threshold) {
			if (percentGap >= percentMask) {
				*(seq + j) = '-';
				information = percentGap;
			} else {
				*(seq + j) = '+';
				information = percentMask;
			}
		}
		if (information < minInfo) {
			*(seq + j) = '?'; // no consensus above threshold
		}
	}
}

//ans_start <- .Call("consensusSequence", myDNAStringSet, threshold, ambiguity, minInformation, ignoreNonLetters, terminalGaps, PACKAGE="DECIPHER")
SEXP consensusSequence(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps)
{
	cachedXStringSet x_set;
	cachedCharSeq x_i;
	int x_length, i, j, seqLength, degeneracy, ignore, tGaps;
	double *thresh, *minInfo;
	SEXP consensusSeq;
	
	// initialize the XStringSet
	x_set = cache_XStringSet(x);
	x_length = get_cachedXStringSet_length(&x_set);
	int gapLengths[x_length][2];
	degeneracy = asLogical(ambiguity);
	ignore = asLogical(ignoreNonLetters);
	tGaps = asLogical(terminalGaps);
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_cachedXStringSet_elt(&x_set, i);
		if (x_i.length > seqLength) {
			seqLength = x_i.length;
		}
	}
	// initialize an array of encoded base counts
	double bases[7][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 7; i++) {
		for (j = 0; j < seqLength; j++) {
			bases[i][j] = 0;
		}
	}
	
	// loop through each sequence in the DNAStringSet
	for (i = 0; i < x_length; i++) {
		// extract each ith DNAString from the DNAStringSet
		x_i = get_cachedXStringSet_elt(&x_set, i);
		
		// update the alphabet for this string
		if (!tGaps) { // don't include terminal gaps
			gapLengths[i][0] = frontTerminalGaps(&x_i);
			gapLengths[i][1] = endTerminalGaps(&x_i);
			alphabetFrequency(&x_i, &bases[0][0], seqLength, degeneracy, ignore, gapLengths[i][0], gapLengths[i][1], 1);
		} else { // include terminal gaps
			alphabetFrequency(&x_i, &bases[0][0], seqLength, degeneracy, ignore, 0, 0, 1);
		}
	}
	
	thresh = REAL(threshold);
	minInfo = REAL(minInformation);
	char seq[seqLength + 1]; // last position is for null terminating
	makeConsensus(&bases[0][0], &seq[0], seqLength, x_length, *thresh, *minInfo, tGaps);
	seq[seqLength] = '\0'; // end (null terminate) the string
	
	PROTECT(consensusSeq = allocVector(STRSXP, 1));
	SET_STRING_ELT(consensusSeq, 0, mkChar(seq));
	
	UNPROTECT(1);

	return consensusSeq;
}

//ans_start <- .Call("consensusProfile", myDNAStringSet, weight, PACKAGE="DECIPHER")
SEXP consensusProfile(SEXP x, SEXP weight)
{
	cachedXStringSet x_set;
	cachedCharSeq x_i;
	int x_length, i, j, seqLength;
	SEXP ans;//, subM, ret_list
	double *rans, *w = REAL(weight);//, *m
	
	// initialize the XStringSet
	x_set = cache_XStringSet(x);
	x_length = get_cachedXStringSet_length(&x_set);
	double gapLengths[x_length][2];
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_cachedXStringSet_elt(&x_set, i);
		if (x_i.length > seqLength) {
			seqLength = x_i.length;
		}
	}
	// initialize an array of encoded base counts
	double bases[7][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 7; i++) {
		for (j = 0; j < seqLength; j++) {
			bases[i][j] = 0;
		}
	}
	
	// initialize an array of encoded base counts
	double gaps[2][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 2; i++) {
		for (j = 0; j < seqLength; j++) {
			gaps[i][j] = 0;
		}
	}
	
	// loop through each sequence in the DNAStringSet
	for (i = 0; i < x_length; i++) {
		// extract each ith DNAString from the DNAStringSet
		x_i = get_cachedXStringSet_elt(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i][0] = frontTerminalGaps(&x_i);
		gapLengths[i][1] = endTerminalGaps(&x_i);
		alphabetFrequency(&x_i, &bases[0][0], seqLength, 1, 0, gapLengths[i][0], gapLengths[i][1], w[i]);
		
		for (j = gapLengths[i][0] + 1; j < seqLength - gapLengths[i][1] - 1; j++) {
			if (!(x_i.seq[j - 1] & 0x10) && (x_i.seq[j] & 0x10)) {
				gaps[0][j] += w[i]; // gap opening
			} else if ((x_i.seq[j] & 0x10) && !(x_i.seq[j + 1] & 0x10)) {
				gaps[1][j] += w[i]; // gap closing
			}
		}
	}
	
	/*
	PROTECT(subM = allocMatrix(REALSXP, 4, 4));
	m = REAL(subM);
	
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			*(m + i*4 + j) = 0;
		}
	}
	
	if (x_length > 1) {
		for (i = 0; i < seqLength; i++) {
			*(m) = *(m) + bases[0][i]*(bases[0][i] - 1);
			*(m + 5) = *(m + 5) + bases[1][i]*(bases[1][i] - 1);
			*(m + 10) = *(m + 10) + bases[2][i]*(bases[2][i] - 1);
			*(m + 15) = *(m + 15) + bases[3][i]*(bases[3][i] - 1);
			
			*(m + 1) = *(m + 1) + 2*bases[0][i]*bases[1][i];
			*(m + 2) = *(m + 2) + 2*bases[0][i]*bases[2][i];
			*(m + 3) = *(m + 3) + 2*bases[0][i]*bases[3][i];
			*(m + 6) = *(m + 6) + 2*bases[1][i]*bases[2][i];
			*(m + 7) = *(m + 7) + 2*bases[1][i]*bases[3][i];
			*(m + 11) = *(m + 11) + 2*bases[2][i]*bases[3][i];
		}
	}
	*/
	
	PROTECT(ans = allocMatrix(REALSXP, 7, seqLength));
	rans = REAL(ans);
	
	for (i = 0; i < seqLength; i++) {
		*(rans + i*7 + 0) = bases[0][i]/x_length;
		*(rans + i*7 + 1) = bases[1][i]/x_length;
		*(rans + i*7 + 2) = bases[2][i]/x_length;
		*(rans + i*7 + 3) = bases[3][i]/x_length;
		*(rans + i*7 + 4) = bases[4][i]/x_length;
		*(rans + i*7 + 5) = gaps[0][i]/x_length;
		*(rans + i*7 + 6) = gaps[1][i]/x_length;
	}
	
	//PROTECT(ret_list = allocVector(VECSXP, 2));
	//SET_VECTOR_ELT(ret_list, 0, ans);
	//SET_VECTOR_ELT(ret_list, 1, subM);
	
	UNPROTECT(1);
	
	return ans;
}
