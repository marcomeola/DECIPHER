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

/* for Calloc/Free */
#include <R_ext/RS.h>

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

static int frontTerminalGaps(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
	     i < P->length;
	     i++, p++)
	{
		if ((*p) & 0x10 || (*p) & 0x40) { // gap character ("-" or ".")
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGaps(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
	     i >= 0;
	     i--, p--)
	{
		if ((*p) & 0x10 || (*p) & 0x40) { // gap character ("-" or ".")
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int frontTerminalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
	     i < P->length;
	     i++, p++)
	{
		if (!((*p) ^ 0x2D) || !((*p) ^ 0x2E)) { // gap character ("-" or ".")
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
	     i >= 0;
	     i--, p--)
	{
		if (!((*p) ^ 0x2D) || !((*p) ^ 0x2E)) { // gap character ("-" or ".")
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static void alphabetFrequency(const Chars_holder *P, double *bits, int seqLength, int degeneracy, int ignore, int start, int end, double weight)
{
	int j;
	const char *p;
	
	for (j = start, p = (P->ptr + start);
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
				case 64: // . treated as -
					if (ignore==1) { // don't include gaps
						*(bits + 6*seqLength + j) -= weight;
					} else { // include gaps
						*(bits + 4*seqLength + j) += weight;
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
				case 64: // . treated as -
					*(bits + 4*seqLength + j) += weight;
					*(bits + 6*seqLength + j) += weight;
					break;
				default:
					error("not DNA!");
					break;
			}
		}
	}
}

static void alphabetFrequencyAA(const Chars_holder *P, double *bits, int seqLength, int degeneracy, int ignore, int start, int end, double weight)
{
	int j, i;
	const char *p;
	
	for (j = start, p = (P->ptr + start);
		 j < (P->length - end);
		 j++, p++)
	{
		// another base counted
		*(bits + 25*seqLength + j) += weight;
		
		// tally the bases into the encoded array
		switch (*p) {
			case 65: // A
				*(bits + 0*seqLength + j) += weight;
				break;
			case 82: // R
				*(bits + 1*seqLength + j) += weight;
				break;
			case 78: // N
				*(bits + 2*seqLength + j) += weight;
				break;
			case 68: // D
				*(bits + 3*seqLength + j) += weight;
				break;
			case 67: // C
				*(bits + 4*seqLength + j) += weight;
				break;
			case 81: // Q
				*(bits + 5*seqLength + j) += weight;
				break;
			case 69: // E
				*(bits + 6*seqLength + j) += weight;
				break;
			case 71: // G
				*(bits + 7*seqLength + j) += weight;
				break;
			case 72: // H
				*(bits + 8*seqLength + j) += weight;
				break;
			case 73: // I
				*(bits + 9*seqLength + j) += weight;
				break;
			case 76: // L
				*(bits + 10*seqLength + j) += weight;
				break;
			case 75: // K
				*(bits + 11*seqLength + j) += weight;
				break;
			case 77: // M
				*(bits + 12*seqLength + j) += weight;
				break;
			case 70: // F
				*(bits + 13*seqLength + j) += weight;
				break;
			case 80: // P
				*(bits + 14*seqLength + j) += weight;
				break;
			case 83: // S
				*(bits + 15*seqLength + j) += weight;
				break;
			case 84: // T
				*(bits + 16*seqLength + j) += weight;
				break;
			case 87: // W
				*(bits + 17*seqLength + j) += weight;
				break;
			case 89: // Y
				*(bits + 18*seqLength + j) += weight;
				break;
			case 86: // V
				*(bits + 19*seqLength + j) += weight;
				break;
			case 85: // U
				*(bits + 20*seqLength + j) += weight;
				break;
			case 79: // O
				*(bits + 21*seqLength + j) += weight;
				break;
			case 66: // B = N or D
				if (degeneracy==1) { // include degeneracy codes
					*(bits + 2*seqLength + j) += 0.5*weight;
					*(bits + 3*seqLength + j) += 0.5*weight;
				} else { // don't include degeneracy codes
					*(bits + 25*seqLength + j) -= weight;
				}
				break;
			case 90: // Z = Q or E
				if (degeneracy==1) { // include degeneracy codes
					*(bits + 5*seqLength + j) += 0.5*weight;
					*(bits + 6*seqLength + j) += 0.5*weight;
				} else { // don't include degeneracy codes
					*(bits + 25*seqLength + j) -= weight;
				}
				break;
			case 74: // J = I or L
				if (degeneracy==1) { // include degeneracy codes
					*(bits + 9*seqLength + j) += 0.5*weight;
					*(bits + 10*seqLength + j) += 0.5*weight;
				} else { // don't include degeneracy codes
					*(bits + 25*seqLength + j) -= weight;
				}
				break;
			case 88: // X = any letter
				if (degeneracy==1) { // include degeneracy codes
					for (i = 0; i < 20; i++) {
						// split between the 20 canonical amino acids
						*(bits + i*seqLength + j) += weight/20;
					}
				} else { // don't include degeneracy codes
					*(bits + 25*seqLength + j) -= weight;
				}
				break;
			case 42: // * (stop)
				*(bits + 22*seqLength + j) += weight;
				break;
			case 45: // -
				if (ignore==1) { // don't include gaps
					*(bits + 25*seqLength + j) -= weight;
				} else { // include gaps
					*(bits + 23*seqLength + j) += weight;
				}
				break;
			case 43: // +
				if (ignore==1) { // don't include masks
					*(bits + 25*seqLength + j) -= weight;
				} else { // include masks
					*(bits + 24*seqLength + j) += weight;
				}
				break;
			case 46: // . treated as -
				if (ignore==1) { // don't include gaps
					*(bits + 25*seqLength + j) -= weight;
				} else { // include gaps
					*(bits + 23*seqLength + j) += weight;
				}
				break;
			default:
				error("not AA!");
				break;
		}
	}
}

static void adjustFrequency(const char p, double *bits, int seqLength, int degeneracy, int ignore, int j, double weight)
{
	if (degeneracy==1) { // include degeneracy codes
		// another base counted
		*(bits + 6*seqLength + j) += weight;
		
		// tally the bases into the encoded array
		switch (p) {
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
			case 64: // . treated as -
				if (ignore==1) { // don't include gaps
					*(bits + 6*seqLength + j) -= weight;
				} else { // include gaps
					*(bits + 4*seqLength + j) += weight;
				}
				break;
			default:
				error("not DNA!");
				break;
		}
	} else { // don't include degeneracy codes
		switch (p) {
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
			case 64: // . treated as -
				*(bits + 4*seqLength + j) += weight;
				*(bits + 6*seqLength + j) += weight;
				break;
			default:
				error("not DNA!");
				break;
		}
	}
}

static void adjustFrequencyAA(const char p, double *bits, int seqLength, int degeneracy, int ignore, int j, double weight)
{
	int i;
	
	// another base counted
	*(bits + 25*seqLength + j) += weight;
	
	// tally the bases into the encoded array
	switch (p) {
		case 65: // A
			*(bits + 0*seqLength + j) += weight;
			break;
		case 82: // R
			*(bits + 1*seqLength + j) += weight;
			break;
		case 78: // N
			*(bits + 2*seqLength + j) += weight;
			break;
		case 68: // D
			*(bits + 3*seqLength + j) += weight;
			break;
		case 67: // C
			*(bits + 4*seqLength + j) += weight;
			break;
		case 81: // Q
			*(bits + 5*seqLength + j) += weight;
			break;
		case 69: // E
			*(bits + 6*seqLength + j) += weight;
			break;
		case 71: // G
			*(bits + 7*seqLength + j) += weight;
			break;
		case 72: // H
			*(bits + 8*seqLength + j) += weight;
			break;
		case 73: // I
			*(bits + 9*seqLength + j) += weight;
			break;
		case 76: // L
			*(bits + 10*seqLength + j) += weight;
			break;
		case 75: // K
			*(bits + 11*seqLength + j) += weight;
			break;
		case 77: // M
			*(bits + 12*seqLength + j) += weight;
			break;
		case 70: // F
			*(bits + 13*seqLength + j) += weight;
			break;
		case 80: // P
			*(bits + 14*seqLength + j) += weight;
			break;
		case 83: // S
			*(bits + 15*seqLength + j) += weight;
			break;
		case 84: // T
			*(bits + 16*seqLength + j) += weight;
			break;
		case 87: // W
			*(bits + 17*seqLength + j) += weight;
			break;
		case 89: // Y
			*(bits + 18*seqLength + j) += weight;
			break;
		case 86: // V
			*(bits + 19*seqLength + j) += weight;
			break;
		case 85: // U
			*(bits + 20*seqLength + j) += weight;
			break;
		case 79: // O
			*(bits + 21*seqLength + j) += weight;
			break;
		case 66: // B = N or D
			if (degeneracy==1) { // include degeneracy codes
				*(bits + 2*seqLength + j) += 0.5*weight;
				*(bits + 3*seqLength + j) += 0.5*weight;
			} else { // don't include degeneracy codes
				*(bits + 25*seqLength + j) -= weight;
			}
			break;
		case 90: // Z = Q or E
			if (degeneracy==1) { // include degeneracy codes
				*(bits + 5*seqLength + j) += 0.5*weight;
				*(bits + 6*seqLength + j) += 0.5*weight;
			} else { // don't include degeneracy codes
				*(bits + 25*seqLength + j) -= weight;
			}
			break;
		case 74: // J = I or L
			if (degeneracy==1) { // include degeneracy codes
				*(bits + 9*seqLength + j) += 0.5*weight;
				*(bits + 10*seqLength + j) += 0.5*weight;
			} else { // don't include degeneracy codes
				*(bits + 25*seqLength + j) -= weight;
			}
			break;
		case 88: // X = any letter
			if (degeneracy==1) { // include degeneracy codes
				for (i = 0; i < 20; i++) {
					// split between the 20 canonical amino acids
					*(bits + i*seqLength + j) += weight/20;
				}
			} else { // don't include degeneracy codes
				*(bits + 25*seqLength + j) -= weight;
			}
			break;
		case 42: // * (stop)
			*(bits + 22*seqLength + j) += weight;
			break;
		case 45: // -
			if (ignore==1) { // don't include gaps
				*(bits + 25*seqLength + j) -= weight;
			} else { // include gaps
				*(bits + 23*seqLength + j) += weight;
			}
			break;
		case 43: // +
			if (ignore==1) { // don't include masks
				*(bits + 25*seqLength + j) -= weight;
			} else { // include masks
				*(bits + 24*seqLength + j) += weight;
			}
			break;
		case 46: // . treated as -
			if (ignore==1) { // don't include gaps
				*(bits + 25*seqLength + j) -= weight;
			} else { // include gaps
				*(bits + 23*seqLength + j) += weight;
			}
			break;
		default:
			error("not AA!");
			break;
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
		if (percentA > threshold && percentC > threshold && percentG > threshold && percentT > threshold) {
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
		} else if (percentC > threshold && percentG > threshold && percentT > threshold) {
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
		} else if (percentA > threshold && percentG > threshold && percentT > threshold) {
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
		} else if (percentA > threshold && percentC > threshold && percentT > threshold) {
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
		} else if (percentA > threshold && percentC > threshold && percentG > threshold) {
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
		} else if (percentC > threshold && percentT > threshold) {
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
		} else if (percentG > threshold && percentT > threshold) {
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
		} else if (percentA > threshold && percentT > threshold) {
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
		} else if (percentC > threshold && percentG > threshold) {
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
		} else if (percentA > threshold && percentG > threshold) {
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
		} else if (percentA > threshold && percentC > threshold) {
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
		} else if (percentA > threshold) {
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
		} else if (percentC > threshold) {
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
		} else if (percentG > threshold) {
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
		} else if (percentT > threshold) {
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
		} else if (percentGap > threshold ||
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

static void makeConsensusAA(double *bits, char *seq, int seqLength, int x_length, double threshold, double minInfo, int tGaps)
{
	int j, i;
	double information, AAs[23], S, percentGap, percentMask;
	
	for (j = 0; j < seqLength; j++) {
		
		if (!tGaps && (*(bits + 25*seqLength + j)==0)) {
			*(seq + j) = '-';
			continue;
		}
		
		// find the percentage of each base in position j
		S = 0;
		for (i = 0; i < 23; i++) {
			AAs[i] = (double)(*(bits + i*seqLength + j))/(*(bits + 25*seqLength + j));
			S += AAs[i];
		}
		percentGap = (double)(*(bits + 23*seqLength + j))/(*(bits + 25*seqLength + j));
		percentMask = (double)(*(bits + 24*seqLength + j))/(*(bits + 25*seqLength + j));
		information = 0;
		
		// determine the most common base over the threshold percentage
		if (AAs[0] >= threshold) {
			if (AAs[0] >= percentGap &&
				AAs[0] >= percentMask) {
				*(seq + j) = 'A';
				information = AAs[0];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[1] >= threshold) {
			if (AAs[1] >= percentGap &&
				AAs[1] >= percentMask) {
				*(seq + j) = 'R';
				information = AAs[1];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[2] >= threshold) {
			if (AAs[2] >= percentGap &&
				AAs[2] >= percentMask) {
				*(seq + j) = 'N';
				information = AAs[2];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[3] >= threshold) {
			if (AAs[3] >= percentGap &&
				AAs[3] >= percentMask) {
				*(seq + j) = 'D';
				information = AAs[3];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[4] >= threshold) {
			if (AAs[4] >= percentGap &&
				AAs[4] >= percentMask) {
				*(seq + j) = 'C';
				information = AAs[4];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[5] >= threshold) {
			if (AAs[5] >= percentGap &&
				AAs[5] >= percentMask) {
				*(seq + j) = 'Q';
				information = AAs[5];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[6] >= threshold) {
			if (AAs[6] >= percentGap &&
				AAs[6] >= percentMask) {
				*(seq + j) = 'E';
				information = AAs[6];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[7] >= threshold) {
			if (AAs[7] >= percentGap &&
				AAs[7] >= percentMask) {
				*(seq + j) = 'G';
				information = AAs[7];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[8] >= threshold) {
			if (AAs[8] >= percentGap &&
				AAs[8] >= percentMask) {
				*(seq + j) = 'H';
				information = AAs[8];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[9] >= threshold) {
			if (AAs[9] >= percentGap &&
				AAs[9] >= percentMask) {
				*(seq + j) = 'I';
				information = AAs[9];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[10] >= threshold) {
			if (AAs[10] >= percentGap &&
				AAs[10] >= percentMask) {
				*(seq + j) = 'L';
				information = AAs[10];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[11] >= threshold) {
			if (AAs[11] >= percentGap &&
				AAs[11] >= percentMask) {
				*(seq + j) = 'K';
				information = AAs[11];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[12] >= threshold) {
			if (AAs[12] >= percentGap &&
				AAs[12] >= percentMask) {
				*(seq + j) = 'M';
				information = AAs[12];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[13] >= threshold) {
			if (AAs[13] >= percentGap &&
				AAs[13] >= percentMask) {
				*(seq + j) = 'F';
				information = AAs[13];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[14] >= threshold) {
			if (AAs[14] >= percentGap &&
				AAs[14] >= percentMask) {
				*(seq + j) = 'P';
				information = AAs[14];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[15] >= threshold) {
			if (AAs[15] >= percentGap &&
				AAs[15] >= percentMask) {
				*(seq + j) = 'S';
				information = AAs[15];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[16] >= threshold) {
			if (AAs[16] >= percentGap &&
				AAs[16] >= percentMask) {
				*(seq + j) = 'T';
				information = AAs[16];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[17] >= threshold) {
			if (AAs[17] >= percentGap &&
				AAs[17] >= percentMask) {
				*(seq + j) = 'W';
				information = AAs[17];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[18] >= threshold) {
			if (AAs[18] >= percentGap &&
				AAs[18] >= percentMask) {
				*(seq + j) = 'Y';
				information = AAs[18];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[19] >= threshold) {
			if (AAs[19] >= percentGap &&
				AAs[19] >= percentMask) {
				*(seq + j) = 'V';
				information = AAs[19];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[20] >= threshold) {
			if (AAs[20] >= percentGap &&
				AAs[20] >= percentMask) {
				*(seq + j) = 'U';
				information = AAs[20];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[21] >= threshold) {
			if (AAs[21] >= percentGap &&
				AAs[21] >= percentMask) {
				*(seq + j) = 'O';
				information = AAs[21];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (AAs[22] >= threshold) {
			if (AAs[22] >= percentGap &&
				AAs[22] >= percentMask) {
				*(seq + j) = '*';
				information = AAs[22];
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if ((AAs[2] + AAs[3]) >= threshold) {
			if ((AAs[2] + AAs[3]) >= percentGap &&
				(AAs[2] + AAs[3]) >= percentMask) {
				*(seq + j) = 'B'; // B = N (2) or D (3)
				information = (AAs[2] + AAs[3]);
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if ((AAs[5] + AAs[6]) >= threshold) {
			if ((AAs[5] + AAs[6]) >= percentGap &&
				(AAs[5] + AAs[6]) >= percentMask) {
				*(seq + j) = 'Z'; // Z = Q (5) or E (6)
				information = (AAs[5] + AAs[6]);
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if ((AAs[9] + AAs[10]) >= threshold) {
			if ((AAs[9] + AAs[10]) >= percentGap &&
				(AAs[9] + AAs[10]) >= percentMask) {
				*(seq + j) = 'J'; // J = I (9) or L (10)
				information = (AAs[9] + AAs[10]);
			} else {
				if (percentGap >= percentMask) {
					*(seq + j) = '-';
					information = percentGap;
				} else {
					*(seq + j) = '+';
					information = percentMask;
				}
			}
		} else if (S >= threshold) {
			if (S >= percentGap &&
				S >= percentMask) {
				*(seq + j) = 'X';
				information = S;
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
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, seqLength, degeneracy, ignore, tGaps;
	double *thresh, *minInfo;
	SEXP consensusSeq;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	//int gapLengths[x_length][2];
	int gapL, gapR;
	degeneracy = asLogical(ambiguity);
	ignore = asLogical(ignoreNonLetters);
	tGaps = asLogical(terminalGaps);
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (x_i.length > seqLength) {
			seqLength = x_i.length;
		}
	}
	/*
	// initialize an array of encoded base counts
	double bases[7][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 7; i++) {
		for (int j = 0; j < seqLength; j++) {
			bases[i][j] = 0;
		}
	}
	*/
	// initialize an array of encoded base counts
	double *bases = Calloc(7*seqLength, double); // initialized to zero
	
	// loop through each sequence in the DNAStringSet
	for (i = 0; i < x_length; i++) {
		// extract each ith DNAString from the DNAStringSet
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		if (!tGaps) { // don't include terminal gaps
			//gapLengths[i][0] = frontTerminalGaps(&x_i);
			//gapLengths[i][1] = endTerminalGaps(&x_i);
			gapL = frontTerminalGaps(&x_i);
			gapR = endTerminalGaps(&x_i);
			//alphabetFrequency(&x_i, &bases[0][0], seqLength, degeneracy, ignore, gapLengths[i][0], gapLengths[i][1], 1);
			alphabetFrequency(&x_i, bases, seqLength, degeneracy, ignore, gapL, gapR, 1);
		} else { // include terminal gaps
			//alphabetFrequency(&x_i, &bases[0][0], seqLength, degeneracy, ignore, 0, 0, 1);
			alphabetFrequency(&x_i, bases, seqLength, degeneracy, ignore, 0, 0, 1);
		}
	}
	
	thresh = REAL(threshold);
	minInfo = REAL(minInformation);
	char seq[seqLength + 1]; // last position is for null terminating
	//makeConsensus(&bases[0][0], &seq[0], seqLength, x_length, *thresh, *minInfo, tGaps);
	makeConsensus(bases, &seq[0], seqLength, x_length, *thresh, *minInfo, tGaps);
	seq[seqLength] = '\0'; // end (null terminate) the string
	
	PROTECT(consensusSeq = allocVector(STRSXP, 1));
	SET_STRING_ELT(consensusSeq, 0, mkChar(seq));
	
	Free(bases);
	
	UNPROTECT(1);

	return consensusSeq;
}

//ans_start <- .Call("consensusSequenceAA", myAAStringSet, threshold, ambiguity, minInformation, ignoreNonLetters, terminalGaps, PACKAGE="DECIPHER")
SEXP consensusSequenceAA(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, seqLength, degeneracy, ignore, tGaps;
	double thresh, minInfo;
	SEXP consensusSeq;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	//int gapLengths[x_length][2];
	int gapL, gapR;
	degeneracy = asLogical(ambiguity);
	ignore = asLogical(ignoreNonLetters);
	tGaps = asLogical(terminalGaps);
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (x_i.length > seqLength) {
			seqLength = x_i.length;
		}
	}
	/*
	// initialize an array of encoded base counts
	double bases[26][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 26; i++) {
		for (int j = 0; j < seqLength; j++) {
			bases[i][j] = 0;
		}
	}
	 */
	// initialize an array of encoded base counts
	double *bases = Calloc(26*seqLength, double); // initialized to zero
	
	// loop through each sequence in the AAStringSet
	for (i = 0; i < x_length; i++) {
		// extract each ith AAString from the AAStringSet
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		if (!tGaps) { // don't include terminal gaps
			//gapLengths[i][0] = frontTerminalGapsAA(&x_i);
			//gapLengths[i][1] = endTerminalGapsAA(&x_i);
			gapL = frontTerminalGapsAA(&x_i);
			gapR = endTerminalGapsAA(&x_i);
			//alphabetFrequencyAA(&x_i, &bases[0][0], seqLength, degeneracy, ignore, gapLengths[i][0], gapLengths[i][1], 1);
			alphabetFrequencyAA(&x_i, bases, seqLength, degeneracy, ignore, gapL, gapR, 1);
		} else { // include terminal gaps
			//alphabetFrequencyAA(&x_i, &bases[0][0], seqLength, degeneracy, ignore, 0, 0, 1);
			alphabetFrequencyAA(&x_i, bases, seqLength, degeneracy, ignore, 0, 0, 1);
		}
	}
	
	thresh = asReal(threshold);
	thresh = 1 - thresh;
	minInfo = asReal(minInformation);
	char seq[seqLength + 1]; // last position is for null terminating
	//makeConsensusAA(&bases[0][0], &seq[0], seqLength, x_length, thresh, minInfo, tGaps);
	makeConsensusAA(bases, &seq[0], seqLength, x_length, thresh, minInfo, tGaps);
	seq[seqLength] = '\0'; // end (null terminate) the string
	
	PROTECT(consensusSeq = allocVector(STRSXP, 1));
	SET_STRING_ELT(consensusSeq, 0, mkChar(seq));
	
	Free(bases);
	
	UNPROTECT(1);
	
	return consensusSeq;
}

//ans_start <- .Call("consensusProfile", myDNAStringSet, weight, PACKAGE="DECIPHER")
SEXP consensusProfile(SEXP x, SEXP weight)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, seqLength;
	SEXP ans;//, subM, ret_list
	double *rans, *w = REAL(weight), sum, tot = 0;//, *m
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	//int gapLengths[x_length][2];
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (x_i.length > seqLength) {
			seqLength = x_i.length;
		}
	}
	
	/*
	// initialize an array of encoded base counts
	double bases[7][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 7; i++) {
		for (j = 0; j < seqLength; j++) {
			bases[i][j] = 0;
		}
	}
	
	// initialize an array for gap opening and closing
	double gaps[2][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 2; i++) {
		for (j = 0; j < seqLength; j++) {
			gaps[i][j] = 0;
		}
	}
	
	// initialize an array of positional weights
	double totW[seqLength + 1];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (j = 0; j <= seqLength; j++)
		totW[j] = 0;
	*/
	// initialize an array of encoded base counts
	double *bases = Calloc(7*seqLength, double); // initialized to zero
	// initialize an array for gap opening and closing
	double *gaps = Calloc(2*seqLength, double); // initialized to zero
	// initialize an array of positional weights
	double *totW = Calloc(seqLength + 1, double); // initialized to zero
	// initialize an array of terminal gap lengths
	int *gapLengths = Calloc(x_length*2, int); // initialized to zero
	
	// loop through each sequence in the DNAStringSet
	for (i = 0; i < x_length; i++) {
		if (w[i]==0)
			continue;
		
		// extract each ith DNAString from the DNAStringSet
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i*2] = frontTerminalGaps(&x_i);
		if (gapLengths[i*2]==x_i.length) // all gaps
			continue;
		gapLengths[i*2 + 1] = endTerminalGaps(&x_i);
		alphabetFrequency(&x_i, bases, seqLength, 1, 0, gapLengths[i*2], gapLengths[i*2 + 1], w[i]);
		
		totW[gapLengths[i*2]] += w[i];
		totW[seqLength - gapLengths[i*2 + 1]] -= w[i];
		
		for (j = gapLengths[i*2] + 1; j < seqLength - gapLengths[i*2 + 1] - 1; j++) {
			if (!(x_i.ptr[j - 1] & 0x10 || x_i.ptr[j - 1] & 0x40) && (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40)) {
				gaps[2*j] += w[i]; // gap opening
			}
			if ((x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40) && !(x_i.ptr[j + 1] & 0x10 || x_i.ptr[j + 1] & 0x40)) {
				gaps[2*j + 1] += w[i]; // gap closing
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
	
	PROTECT(ans = allocMatrix(REALSXP, 8, seqLength));
	rans = REAL(ans);
	
	for (i = 0; i < seqLength; i++) {
		tot += totW[i];
		if (tot==0) {
			for (j = 0; j < 8; j++)
				*(rans + i*8 + j) = 0;
		} else {
			*(rans + i*8 + 0) = bases[i]/tot;
			*(rans + i*8 + 1) = bases[1*seqLength + i]/tot;
			*(rans + i*8 + 2) = bases[2*seqLength + i]/tot;
			*(rans + i*8 + 3) = bases[3*seqLength + i]/tot;
			*(rans + i*8 + 4) = bases[4*seqLength + i]/tot;
			*(rans + i*8 + 5) = gaps[i*2]/tot;
			*(rans + i*8 + 6) = gaps[i*2 + 1]/tot;
			sum = *(rans + i*8) + *(rans + i*8 + 1) + *(rans + i*8 + 2) + *(rans + i*8 + 3) + *(rans + i*8 + 4);
			if (sum > 0) { // normalize the profile
				*(rans + i*8) /= sum;
				*(rans + i*8 + 1) /= sum;
				*(rans + i*8 + 2) /= sum;
				*(rans + i*8 + 3) /= sum;
				*(rans + i*8 + 4) /= sum;
			}
			*(rans + i*8 + 7) = tot/x_length;
		}
	}
	
	Free(bases);
	Free(gaps);
	Free(totW);
	Free(gapLengths);
	
	//PROTECT(ret_list = allocVector(VECSXP, 2));
	//SET_VECTOR_ELT(ret_list, 0, ans);
	//SET_VECTOR_ELT(ret_list, 1, subM);
	
	UNPROTECT(1);
	
	return ans;
}

//ans_start <- .Call("consensusProfileAA", myAAStringSet, weight, NULL, PACKAGE="DECIPHER")
SEXP consensusProfileAA(SEXP x, SEXP weight, SEXP structs)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, seqLength;
	SEXP ans, elmt, dims;
	double *rans, *w = REAL(weight), sum, tot = 0;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	//int gapLengths[x_length][2];
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (x_i.length > seqLength) {
			seqLength = x_i.length;
		}
	}
	
	// initialize an array of HEC counts
	double *HEC, *s;
	int do_HEC, n, l, d, size = 27;
	if (length(structs)==x_length) {
		do_HEC = 1;
		elmt = VECTOR_ELT(structs, 0);
		PROTECT(dims = GET_DIM(elmt));
		d = INTEGER(dims)[0];
		UNPROTECT(1);
		HEC = Calloc(d*seqLength, double);
		size += d;
	} else {
		do_HEC = 0;
		d = 0;
	}
	
	/*
	// initialize an array of encoded base counts
	double bases[26][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 26; i++) {
		for (j = 0; j < seqLength; j++) {
			bases[i][j] = 0;
		}
	}
	
	// initialize an array for gap opening and closing
	double gaps[2][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 2; i++) {
		for (j = 0; j < seqLength; j++) {
			gaps[i][j] = 0;
		}
	}
	
	// initialize an array of positional weights
	double totW[seqLength + 1];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (j = 0; j <= seqLength; j++)
		totW[j] = 0;
	*/
	// initialize an array of encoded base counts
	double *bases = Calloc(26*seqLength, double); // initialized to zero
	// initialize an array for gap opening and closing
	double *gaps = Calloc(2*seqLength, double); // initialized to zero
	// initialize an array of positional weights
	double *totW = Calloc(seqLength + 1, double); // initialized to zero
	// initialize an array of terminal gap lengths
	int *gapLengths = Calloc(x_length*2, int); // initialized to zero
	
	// loop through each sequence in the AAStringSet
	for (i = 0; i < x_length; i++) {
		if (w[i]==0)
			continue;
		
		// extract each ith AAString from the AAStringSet
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i*2] = frontTerminalGapsAA(&x_i);
		if (gapLengths[i*2]==x_i.length) // all gaps
			continue;
		gapLengths[i*2 + 1] = endTerminalGapsAA(&x_i);
		alphabetFrequencyAA(&x_i, bases, seqLength, 1, 0, gapLengths[i*2], gapLengths[i*2 + 1], w[i]);
		
		totW[gapLengths[i*2]] += w[i];
		totW[seqLength - gapLengths[i*2 + 1]] -= w[i];
		
		for (j = gapLengths[i*2] + 1; j < seqLength - gapLengths[i*2 + 1] - 1; j++) {
			if ((x_i.ptr[j - 1] ^ 0x2D && x_i.ptr[j - 1] ^ 0x2E) && (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E))) {
				gaps[2*j] += w[i]; // gap opening
			}
			if ((x_i.ptr[j + 1] ^ 0x2D && x_i.ptr[j + 1] ^ 0x2E) && (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E))) {
				gaps[2*j + 1] += w[i]; // gap closing
			}
		}
		
		if (do_HEC) {
			elmt = VECTOR_ELT(structs, i);
			s = REAL(elmt);
			l = length(elmt);
			n = 0;
			for (j = gapLengths[i*2]; j < seqLength - gapLengths[i*2 + 1]; j++) {
				if (x_i.ptr[j] != 45 && x_i.ptr[j] != 46 && x_i.ptr[j] != 43) {
					if (n + d > l)
						error("Structure does not match the sequence.");
					for (k = 0; k < d; k++)
						HEC[j + k*seqLength] += s[n++]*w[i];
				}
			}
		}
	}
	
	PROTECT(ans = allocMatrix(REALSXP, size, seqLength));
	rans = REAL(ans);
	
	for (i = 0; i < seqLength; i++) {
		tot += totW[i];
		
		if (tot==0) {
			for (j = 0; j < size; j++)
				*(rans + i*size + j) = 0;
		} else {
			for (j = 0; j < 24; j++) {
				*(rans + i*size + j) = bases[j*seqLength + i]/tot;
			}
			sum = 0;
			for (j = 0; j < 24; j++) {
				sum += *(rans + i*size + j);
			}
			if (sum > 0) {
				for (j = 0; j < 24; j++) {
					*(rans + i*size + j) /= sum;
				}
			}
			*(rans + i*size + 24) = gaps[i*2]/tot;
			*(rans + i*size + 25) = gaps[i*2 + 1]/tot;
			*(rans + i*size + 26) = tot/x_length;
			
			if (do_HEC) {
				for (j = 0; j < d; j++) {
					*(rans + i*size + j + 27) = HEC[i + j*seqLength]/tot;
				}
			}
		}
	}
	
	Free(bases);
	Free(gaps);
	Free(totW);
	Free(gapLengths);
	if (do_HEC)
		Free(HEC);
	
	UNPROTECT(1);
	
	return ans;
}

// returns the sum of substitution scores
// for each alignment column [start-end]
SEXP colScores(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP weights)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, k, i, j, seqLength;
	double *subM = REAL(subMatrix);
	double *w = REAL(weights);
	double GO = asReal(go); // gap opening
	double GE = asReal(ge); // gap extension
	double weight, total, prev, curr, prev_total;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	seqLength = x_i.length;
	
	// initialize an array of encoded base counts
	double *bases = Calloc(7*seqLength, double); // initialized to zero
	// initialize an array of terminal gap lengths
	int *gapLengths = Calloc(x_length*2, int); // initialized to zero
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i*2] = frontTerminalGaps(&x_i);
		if (gapLengths[i*2]==x_i.length) // all gaps
			continue;
		gapLengths[i*2 + 1] = endTerminalGaps(&x_i);
		alphabetFrequency(&x_i, bases, seqLength, 1, 0, gapLengths[i*2], gapLengths[i*2 + 1], w[i]);
	}
	
	SEXP ans;
	double *rans;
	PROTECT(ans = allocVector(REALSXP, seqLength));
	rans = REAL(ans);
	
	prev = 0; // gaps in previous position
	prev_total = 0; // number of letters in previous position
	for (k = 0; k < seqLength; k++) {
		*(rans + k) = 0;
		total = 0;
		for (i = 0; i < 4; i++) {
			total += bases[i*seqLength + k];
			for (j = i; j < 4; j++) {
				weight = (i==j) ? 0.5 : 1;
				weight *= (bases[j*seqLength + k] - ((i==j) ? 1 : 0));
				weight *= bases[i*seqLength + k];
				if (weight > 0)
					*(rans + k) += *(subM + i*4 + j)*weight;
			}
		}
		
		curr = bases[4*seqLength + k]; // number gapped
		if (curr > prev) {
			*(rans + k) += GO*((curr - prev)*total); // gap opening
			*(rans + k) += GE*(prev*total); // gap extension
		} else if (curr < prev) {
			*(rans + k) += GO*((prev - curr)*prev_total); // gap closing
			*(rans + k) += GE*(curr*total); // gap extension
		} else {
			*(rans + k) += GE*(curr*total); // gap extension
		}
		
		prev = curr;
		prev_total = total;
	}
	
	Free(bases);
	Free(gapLengths);
	
	UNPROTECT(1);
	
	return ans;
}

// returns the sum of substitution scores
// for each alignment column [start-end]
SEXP colScoresAA(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP weights, SEXP structs, SEXP hecMatrix)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, k, i, j, seqLength;
	double *subM = REAL(subMatrix);
	double *w = REAL(weights);
	double GO = asReal(go); // gap opening
	double GE = asReal(ge); // gap extension
	double weight, total, prev, curr, prev_total;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	seqLength = x_i.length;
	
	// initialize an array of HEC counts
	SEXP elmt, dims;
	double *HEC, *s;
	int do_HEC, n, l, d;
	double *hecM = REAL(hecMatrix);
	if (length(structs)==x_length) {
		do_HEC = 1;
		elmt = VECTOR_ELT(structs, 0);
		PROTECT(dims = GET_DIM(elmt));
		d = INTEGER(dims)[0];
		UNPROTECT(1);
	} else {
		do_HEC = 0;
	}
	
	// initialize an array of encoded base counts
	double *bases = Calloc(26*seqLength, double); // initialized to zero
	// initialize an array of terminal gap lengths
	int *gapLengths = Calloc(x_length*2, int); // initialized to zero
	if (do_HEC) // initialize an array of structure counts
		HEC = Calloc(d*seqLength, double); // initialized to zero
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i*2] = frontTerminalGapsAA(&x_i);
		if (gapLengths[i*2]==x_i.length) // all gaps
			continue;
		gapLengths[i*2 + 1] = endTerminalGapsAA(&x_i);
		alphabetFrequencyAA(&x_i, bases, seqLength, 1, 0, gapLengths[i*2], gapLengths[i*2 + 1], w[i]);
		
		if (do_HEC) {
			elmt = VECTOR_ELT(structs, i);
			s = REAL(elmt);
			l = length(elmt);
			n = 0;
			for (j = gapLengths[i*2]; j < seqLength - gapLengths[i*2 + 1]; j++) {
				if (x_i.ptr[j] != 45 && x_i.ptr[j] != 46 && x_i.ptr[j] != 43) {
					if (n + d > l)
						error("Structure does not match the sequence.");
					for (k = 0; k < d; k++)
						HEC[j + k*seqLength] += s[n++]*w[i];
				}
			}
		}
	}
	
	SEXP ans;
	double *rans;
	PROTECT(ans = allocVector(REALSXP, seqLength));
	rans = REAL(ans);
	
	prev = 0; // gaps in previous position
	prev_total = 0; // number of letters in previous position
	for (k = 0; k < seqLength; k++) {
		*(rans + k) = 0;
		total = 0;
		for (i = 0; i < 20; i++) {
			total += bases[i*seqLength + k];
			for (j = i; j < 20; j++) {
				weight = (i==j) ? 0.5 : 1;
				weight *= (bases[j*seqLength + k] - ((i==j) ? 1 : 0));
				weight *= bases[i*seqLength + k];
				if (weight > 0)
					*(rans + k) += *(subM + i*21 + j)*weight;
			}
		}
		
		if (do_HEC) {
			for (i = 0; i < d; i++) {
				for (j = i; j < d; j++) {
					weight = (i==j) ? 0.5 : 1;
					weight *= (HEC[j*seqLength + k] - ((i==j) ? 1 : 0));
					weight *= HEC[i*seqLength + k];
					if (weight > 0)
						*(rans + k) += *(hecM + i*d + j)*weight;
				}
			}
		}
		
		curr = bases[23*seqLength + k]; // number gapped
		if (curr > prev) {
			*(rans + k) += GO*((curr - prev)*total); // gap opening
			*(rans + k) += GE*(prev*total); // gap extension
		} else if (curr < prev) {
			*(rans + k) += GO*((prev - curr)*prev_total); // gap closing
			*(rans + k) += GE*(curr*total); // gap extension
		} else {
			*(rans + k) += GE*(curr*total); // gap extension
		}
		
		prev = curr;
		prev_total = total;
	}
	
	Free(bases);
	Free(gapLengths);
	if (do_HEC)
		Free(HEC);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP shiftGaps(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP gl, SEXP sc, SEXP thresh, SEXP weights)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, seqLength, weight, tot, pos, bound, end;
	char p;
	unsigned long long int temp;
	int maxSize = sizeof(unsigned long long int)*8; // number of bits in an int
	double *subM = REAL(subMatrix);
	double *w = REAL(weights);
	double GO = asReal(go); // gap opening
	double GE = asReal(ge); // gap extension
	double GL = asReal(gl); // gap-letter mismatch
	double SC = asReal(sc); // shift cost
	double threshold = asReal(thresh); // increase in score to commit changes
	double fractionGaps, total, sum;
	int maxShift = maxSize - 1; // max distance to shift gap events to the left or right
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	seqLength = x_i.length; // sequences are assumed to be aligned (equal length)
	
	// return the number of gap-event shifts made to the sequences
	SEXP ans;
	double *rans;
	PROTECT(ans = allocVector(REALSXP, 1));
	rans = REAL(ans);
	*(rans) = 0; // initialize to zero changes
	
	// initialize an array of encoded base counts
	double *bases = Calloc(7*seqLength, double); // initialized to zero
	// initialize an array of encoded gap events
	unsigned long long int *gaps = Calloc(seqLength, unsigned long long int); // initialized to zero
	// initialize an array of terminal gap lengths
	int *gapLengths = Calloc(x_length*2, int); // initialized to zero
	// initialize an array of gap counts of each length
	int *gapCount = Calloc(maxSize + 1, int); // initialized to zero
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i*2] = frontTerminalGaps(&x_i);
		if (gapLengths[i*2]==x_i.length) // all gaps
			continue;
		gapLengths[i*2 + 1] = endTerminalGaps(&x_i);
		alphabetFrequency(&x_i, bases, seqLength, 1, 0, gapLengths[i*2], gapLengths[i*2 + 1], w[i]);
		
		// bit-encode the length of gaps in an integer vector:
		// 1's for each length of gap closing at that position
		pos = 0;
		for (j = gapLengths[i*2] + 1; j < seqLength - gapLengths[i*2 + 1] - 1; j++) {
			if (pos==0 && // not inside a gap - check for opening
				!(x_i.ptr[j - 1] & 0x10 || x_i.ptr[j - 1] & 0x40) && (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40)) {
				pos = j; // gap opening
			}
			if (pos > 0 && // inside a gap - check for closing
				!(x_i.ptr[j + 1] & 0x10 || x_i.ptr[j + 1] & 0x40) && (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40)) {
				pos = j - pos + 1; // change in position
				if (pos > maxSize) // bigger gap than can be stored
					pos = maxSize;
				gaps[j] |= (unsigned long long int)1 << pos; // encode length of gap
				pos = 0;
			}
		}
	}
	
	// initialize an array of column scores
	double *scores = Calloc(seqLength, double); // initialized to zero
	
	// populate column scores and count gap events
	for (k = 0; k < seqLength; k++) {
		total = 0;
		for (i = 0; i < 4; i++) {
			total += bases[i*seqLength + k];
			for (j = i; j < 4; j++) {
				sum = (i==j) ? 1 : 2;
				sum *= (bases[j*seqLength + k] - ((i==j) ? 1 : 0));
				sum *= bases[i*seqLength + k];
				if (sum > 0)
					scores[k] += *(subM + i*4 + j)*sum;
			}
		}
		scores[k] /= x_length;
		if ((int)total > 1)
			scores[k] /= total - 1; // normalize score
		
		// assume all non-proteogenic letters are gaps
		// penalize maximally for being half gaps
		fractionGaps = (double)(x_length - total)/x_length;
		fractionGaps -= 0.5;
		fractionGaps = 1 - 4*fractionGaps*fractionGaps;
		scores[k] += GL*fractionGaps; // gap-letter mismatch
		
		if (gaps[k] > 0) { // add cost for any gap events
			weight = 0; // length of gap
			temp = gaps[k]; // gets zeroed-out
			while (temp > 0 && weight <= maxSize) {
				if ((temp & 1)==1) {
					scores[k] += weight*GE + GO;
					gapCount[weight]++;
				}
				temp >>= 1;
				weight++;
			}
		}
	}
	
	// ALGORITHM:
	// for each gap size 1 -> maxSize
	//  continue if gapCount[size]==0 (no gaps)
	//  for each sequence
	//   record the number of sequences with a gap of that size at that position:
	//    for each position with a gap of that size
	//     if there is a gap in that position of that size
	//      increment the count for that position
	//  for each position with gaps of that size
	//   order by number of gaps increasing
	//   find all sequences to shift
	//   find the boundaries of the shift
	//    shift gap left and right up to some max distance
	//     stop when the end or another gap is reached
	//     record scores for each shift
	//   if a better score is available
	//    shift gap to that position in the sequence(s)
	//    if gaps were merged then increment gaps and gapCount
	//    if gaps remained the same length then do nothing
	//   gapCount[size]--
	//   break when gapCount[size]==0
	// return the number of changes made
	
	int count, correct, size, min, left, right;
	for (size = 1; size <= maxSize; size++) { // each gap size
		if (gapCount[size]==0) // no gaps of this size
			continue;
		
		// initialize a vector of gap counts of each length
		int *gapNumber = Calloc(gapCount[size], int); // initialized to zero
		// initialize a vector of positions where gap events end (close)
		int *position = Calloc(gapCount[size], int); // initialized to zero
		count = 0;
		for (j = size; j < (seqLength - 1); j++) { // each position
			if ((gaps[j] & ((unsigned long long int)1 << size)) != 0) { // gap of correct size ending at this position
				for (i = 0; i < x_length; i++) { // each sequence
					x_i = get_elt_from_XStringSet_holder(&x_set, i);
					
					correct = 0;
					if (!(x_i.ptr[j + 1] & 0x10 || x_i.ptr[j + 1] & 0x40) && // next position is a letter
						(x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40) && // this position is a gap
						!(x_i.ptr[j - size] & 0x10 || x_i.ptr[j - size] & 0x40)) { // past size-away is a letter
						correct = 1;
						for (k = size; k > 1; k--) {
							if (!(x_i.ptr[j - k + 1] & 0x10 || x_i.ptr[j - k + 1] & 0x40)) { // letter in this position
								correct = 0;
								break;
							}
						}
					}
					if (correct)
						gapNumber[count]++;
				}
				
				position[count] = j;
				count++;
				
				if (count==gapCount[size]) // reached the expected number of gaps
					break;
			}
		}
		
		while (gapCount[size] > 0) {
			// find the gap event with the fewest sequences
			min = 0;
			for (j = 1; j < gapCount[size]; j++)
				if (gapNumber[j] < gapNumber[min])
					min = j;
			if (gapNumber[min]==0) // no gap events found
				break;
			
			// find all sequences to with this gap event
			int *seqNumbers = Calloc(x_length, int); // initialized to zero
			j = position[min];
			count = 0;
			for (i = 0; i < x_length; i++) { // each sequence
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				
				correct = 0;
				if (!(x_i.ptr[j + 1] & 0x10 || x_i.ptr[j + 1] & 0x40) && // next position is a letter
					(x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40) && // this position is a gap
					!(x_i.ptr[j - size] & 0x10 || x_i.ptr[j - size] & 0x40)) { // past size-away is a letter
					correct = 1;
					for (k = size; k > 1; k--) {
						if (!(x_i.ptr[j - k + 1] & 0x10 || x_i.ptr[j - k + 1] & 0x40)) { // letter in this position
							correct = 0;
							break;
						}
					}
				}
				
				if (correct) {
					seqNumbers[count] = i;
					count++;
				}
			}
			gapNumber[min] = count; // reset the number of gaps to account for previous mergers
			
			// first and last position cannot be changed
			left = position[min] - size; // max possible shift left
			left = (left > maxShift) ? maxShift : left;
			right = seqLength - position[min] - 2; // max possible shift right
			right = (right > maxShift) ? maxShift : right;
			for (i = 0; i < gapNumber[min]; i++) { // each sequence
				x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
				
				// prevent a gap from shifting on top of another gap
				
				for (j = 1; j <= left; j++) { // shift offset
					pos = position[min] - size - j;
					if ((x_i.ptr[pos] & 0x10 || x_i.ptr[pos] & 0x40)) { // position is a gap
						left = j; // shift this far left
						break;
					}
				}
				
				for (j = 1; j <= right; j++) { // shift offset
					pos = position[min] + j + 1;
					if ((x_i.ptr[pos] & 0x10 || x_i.ptr[pos] & 0x40)) { // position is a gap
						right = j; // shift this far right
						break;
					}
				}
			}
			
			int bestShift = 0; // negative for left shift, positive for right shift
			double bestScore = 0; // change in score
			
			double totScore = 0; // score for this subset
			int subset = left + size + 1; // additional position is for the far left
			
			// initialize an array of encoded base counts
			double *basesSubset = Calloc(7*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			unsigned long long int *gapsSubset = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			double *scoresSubset = Calloc(subset, double); // initialized to zero
			// initialize an array of encoded base counts
			double *basesSubsetLeftSaved = Calloc(7*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			unsigned long long int *gapsSubsetLeftSaved = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			double *scoresLeftSaved = Calloc(subset, double); // initialized to zero
			
			// copy subset of counts from bases to basesSubset
			for (i = 0; i <= left + size; i++) { // shift offset
				end = position[min] - i;
				for (j = 0; j < 4; j++)
					basesSubset[j*subset + subset - i - 1] = bases[j*seqLength + end];
				gapsSubset[subset - i - 1] = gaps[end];
				scoresSubset[subset - i - 1] = scores[end];
				totScore += scores[end];
			}
			
			// shift the gap event to the left one position at a time
			for (j = 1; j <= left; j++) { // shift offset
				// simulate swapping the left and right bases and gaps
				// right position is (position[min] - j + 1) in sequence and (subset - j) in basesSubset
				// left position is (position[min] - j - size + 1) in sequence and (subset - j - size) in basesSubset
				
				for (i = 0; i < gapNumber[min]; i++) { // each sequence
					x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
										
					// move left position to the right
					// subtract from left position
					adjustFrequency(x_i.ptr[position[min] - j - size + 1], basesSubset, subset, 1, 0, subset - j - size, -1*w[seqNumbers[i]]);
					// add to right position
					adjustFrequency(x_i.ptr[position[min] - j - size + 1], basesSubset, subset, 1, 0, subset - j, w[seqNumbers[i]]);
				}
				
				// re-encode gaps
				// past gap closing is (subset - j) in gapsSubset
				// new gap closing is (subset - j - 1) in gapsSubset
				
				if (j==1) { // clear the original gap of size
					gapsSubset[subset - j] &= ~((unsigned long long int)1 << size);
				} else { // restore the original gaps
					gapsSubset[subset - j] = gaps[position[min] - j + 1];
				}
				// add the new gap
				if (j == left && // j is maximally left
					gapsSubset[subset - j - 1 - size]) { // there is an adjacent gap closing
					
					// initialize a vector of gap sizes to confirm in the adjacent position
					int *confirm = Calloc(maxSize, int); // initialized to zero
					
					// possibly merge gaps
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						bound = position[min] - j - size - maxSize;
						bound = (bound < 0) ? 0 : bound;
						for (k = position[min] - j - size; k >= bound; k--) {
							if (!(x_i.ptr[k] & 0x10 || x_i.ptr[k] & 0x40)) { // letter in this position
								pos = position[min] - j - k; // start - finish - 1 + size
								if (pos - size < maxSize) // need to confirm gap of this size
									confirm[pos - size] = 1;
								if (pos > maxSize)
									pos = maxSize;
								gapsSubset[subset - j - 1] |= ((unsigned long long int)1 << pos);
								break;
							}
						}
					}
					
					int countConfirm = 0;
					for (i = 1; i < maxSize - 1; i++) {
						if (confirm[i] > 0) {
							confirm[countConfirm] = i; // refill with gap lengths to confirm
							countConfirm++;
						}
					}
					
					if (countConfirm > 0) {
						count = 0; // tracks with sequence in seqNumbers
						for (i = 0; i < x_length; i++) {
							if (count < gapNumber[min] &&
								i == seqNumbers[count]) {
								// do not check shifted sequences
								count++;
								continue;
							}
							
							x_i = get_elt_from_XStringSet_holder(&x_set, i);
							
							pos = position[min] - j - size + 1;
							if ((x_i.ptr[pos] & 0x10 || x_i.ptr[pos] & 0x40)) // gap in this position
								continue; // not a gap closing event
							
							int confirmed = 0; // the number of left gaps in this sequence
							bound = position[min] - j - size - maxSize;
							bound = (bound < 0) ? 0 : bound;
							for (pos = position[min] - j - size; pos >= bound; pos--) {
								if (!(x_i.ptr[pos] & 0x10 || x_i.ptr[pos] & 0x40)) { // letter in this position
									break;
								} else {
									confirmed++;
								}
							}
							
							if (confirmed) { // at least one gap ending at the adjacent gap position
								for (k = 0; k < countConfirm; k++) {
									if (confirm[k]==confirmed) { // this gap size was confirmed
										for (pos = k; pos < countConfirm; pos++)
											confirm[pos] = confirm[pos + 1]; // shift confirm left
										countConfirm--;
										break;
									}
								}
							}
							
							if (countConfirm==0) // no more gap events to confirm
								break;
						}
						
						// clear gaps that were not confirmed
						for (k = 0; k < countConfirm; k++)
							gapsSubset[subset - j - 1 - size] &= ~((unsigned long long int)1 << confirm[k]);
					}
					
					Free(confirm);
				} else {
					gapsSubset[subset - j - 1] |= ((unsigned long long int)1 << size);
				}
				
				// re-compute score in the changed positions
				for (pos = 0; pos < subset; pos++) {
					if (!(pos==(subset - j) ||
						  pos==(subset - j - 1) ||
						  pos==(subset - j - size) ||
						  pos==(subset - j - 1 - size)))
						continue; // nothing changed in this position
					
					scoresSubset[pos] = 0; // reset the score in this column
					total = 0;
					for (i = 0; i < 4; i++) {
						total += basesSubset[i*subset + pos];
						for (k = i; k < 4; k++) {
							sum = (i==k) ? 1 : 2;
							sum *= (basesSubset[k*subset + pos] - ((i==k) ? 1 : 0));
							sum *= basesSubset[i*subset + pos];
							if (sum > 0)
								scoresSubset[pos] += *(subM + i*4 + k)*sum;
						}
					}
					scoresSubset[pos] /= x_length;
					if ((int)total > 1)
						scoresSubset[pos] /= total - 1; // normalize score
					
					// assume all non-proteogenic letters are gaps
					// penalize maximally for being half gaps
					fractionGaps = (double)(x_length - total)/x_length;
					fractionGaps -= 0.5;
					fractionGaps = 1 - 4*fractionGaps*fractionGaps;
					scoresSubset[pos] += GL*fractionGaps; // gap-letter mismatch
					
					if (gapsSubset[pos] > 0) { // add cost for any gap events
						weight = 0; // length of gap
						temp = gapsSubset[pos]; // gets zeroed-out
						while (temp > 0 && weight <= maxSize) {
							if ((temp & 1)==1) {
								scoresSubset[pos] += weight*GE + GO;
							}
							temp >>= 1;
							weight++;
						}
					}
				}
				
				double currentScore = 0;
				for (k = 0; k < subset; k++) {
					currentScore += scoresSubset[k];
				}
				
				currentScore += SC*j; // add shift penalty (offset cost)
				
				// if score improved then replace
				if ((currentScore - totScore) > bestScore) {
					bestScore = currentScore - totScore;
					bestShift = -1*j;
					memcpy(basesSubsetLeftSaved, basesSubset, (7*subset) * sizeof(double));
					memcpy(gapsSubsetLeftSaved, gapsSubset, subset * sizeof(unsigned long long int));
					memcpy(scoresLeftSaved, scoresSubset, subset * sizeof(double));
				}
			}
			
			Free(basesSubset);
			Free(gapsSubset);
			Free(scoresSubset);
			
			// repeat process for right shifts
			
			totScore = 0; // score for this subset
			// subset covers from just left of the gap (position[min] - size) to (position[min] + right + maxSize)
			subset = right + size + maxSize; // additional maxSize positions are for the right side
			if ((position[min] - size + 1 + subset) > seqLength)
				subset = seqLength - (position[min] - size) - 1;
			
			// initialize an array of encoded base counts
			basesSubset = Calloc(7*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			gapsSubset = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			scoresSubset = Calloc(subset, double); // initialized to zero
			// initialize an array of encoded base counts
			double *basesSubsetRightSaved = Calloc(7*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			unsigned long long int *gapsSubsetRightSaved = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			double *scoresRightSaved = Calloc(subset, double); // initialized to zero
			
			// copy subset of counts from bases to basesSubset
			for (i = 0; i < subset; i++) {
				end = position[min] - size + i;
				for (j = 0; j < 4; j++)
					basesSubset[j*subset + i] = bases[j*seqLength + end];
				gapsSubset[i] = gaps[end];
				scoresSubset[i] = scores[end];
				totScore += scores[end];
			}
			
			// shift the gap event to the right one position at a time
			for (j = 1; j <= right; j++) { // shift offset
				// simulate swapping the left and right bases and gaps
				
				for (i = 0; i < gapNumber[min]; i++) { // each sequence
					x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
					
					// move left position to the right
					// subtract from right position
					adjustFrequency(x_i.ptr[position[min] + j], basesSubset, subset, 1, 0, j + size, -1*w[seqNumbers[i]]);
					// add to left position
					adjustFrequency(x_i.ptr[position[min] + j], basesSubset, subset, 1, 0, j, w[seqNumbers[i]]);
				}
				
				// re-encode gaps
				// past gap closing is (j - 1 + size) in gapsSubset
				// new gap closing is (j + size) in gapsSubset
				
				if (j==1) { // clear the original gap of size
					gapsSubset[j - 1 + size] &= ~((unsigned long long int)1 << size);
				} else { // restore the original gaps
					gapsSubset[j - 1 + size] = gaps[position[min] + j - 1];
				}
				// add the new gap
				if (j == right) { // j is maximally right
					// check for an adjacent gap opening
					
					// initialize a vector of gap sizes to confirm in the adjacent position
					int *confirm = Calloc(maxSize, int); // initialized to zero
					
					// possibly merge gaps
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						bound = position[min] + j + maxSize;
						bound = (bound > seqLength) ? seqLength : bound;
						for (k = position[min] + j + 1; k < bound; k++) {
							if (!(x_i.ptr[k] & 0x10 || x_i.ptr[k] & 0x40)) { // letter in this position
								pos = k - (position[min] + j + 1) + size; // finish - start + size
								if (pos - size < maxSize) // need to confirm gap of this size
									confirm[pos - size] = 1;
								if (pos > maxSize)
									pos = maxSize;
								gapsSubset[j + pos] |= ((unsigned long long int)1 << pos);
								break;
							}
						}
					}
					
					int countConfirm = 0;
					for (i = 1; i < maxSize - 1; i++) {
						if (confirm[i] > 0) {
							confirm[countConfirm] = i; // refill with gap lengths to confirm
							countConfirm++;
						}
					}
					
					if (countConfirm > 0) {
						count = 0; // tracks with sequence in seqNumbers
						for (i = 0; i < x_length; i++) {
							if (count < gapNumber[min] &&
								i == seqNumbers[count]) {
								// do not check shifted sequences
								count++;
								continue;
							}
							
							x_i = get_elt_from_XStringSet_holder(&x_set, i);
							
							pos = position[min] + j;
							if (!(!(x_i.ptr[pos] & 0x10 || x_i.ptr[pos] & 0x40) && // letter in this position
								(x_i.ptr[pos + 1] & 0x10 || x_i.ptr[pos + 1] & 0x40))) // gap in this position
								continue; // not a gap closing event
							
							int confirmed = 1; // the number of right gaps in this sequence
							bound = position[min] + j + maxSize;
							bound = (bound > seqLength) ? seqLength : bound;
							for (pos = position[min] + j + 2; pos < bound; pos++) {
								if (!(x_i.ptr[pos] & 0x10 || x_i.ptr[pos] & 0x40)) { // letter in this position
									break;
								} else {
									confirmed++;
								}
							}
							
							if (confirmed) { // at least one gap ending at the adjacent gap position
								for (k = 0; k < countConfirm; k++) {
									if (confirm[k]==confirmed) { // this gap size was confirmed
										for (pos = k; pos < countConfirm; pos++)
											confirm[pos] = confirm[pos + 1]; // shift confirm left
										countConfirm--;
										break;
									}
								}
							}
							
							if (countConfirm==0) // no more gap events to confirm
								break;
						}
						
						// clear gaps that were not confirmed
						for (k = 0; k < countConfirm; k++)
							gapsSubset[j + confirm[k] + size] &= ~((unsigned long long int)1 << confirm[k]);
					}
					
					Free(confirm);
				} else {
					gapsSubset[j + size] |= ((unsigned long long int)1 << size);
				}
				
				// re-compute score in the changed positions
				for (pos = 0; pos < subset; pos++) {
					if (!((j == right && pos > j) ||
						  pos==j ||
						  pos==(j + size) ||
						  pos==(j - 1 + size)))
						continue; // nothing changed in this position
					
					scoresSubset[pos] = 0; // reset the score in this column
					total = 0;
					for (i = 0; i < 4; i++) {
						total += basesSubset[i*subset + pos];
						for (k = i; k < 4; k++) {
							sum = (i==k) ? 1 : 2;
							sum *= (basesSubset[k*subset + pos] - ((i==k) ? 1 : 0));
							sum *= basesSubset[i*subset + pos];
							if (sum > 0)
								scoresSubset[pos] += *(subM + i*4 + k)*sum;
						}
					}
					scoresSubset[pos] /= x_length;
					if ((int)total > 1)
						scoresSubset[pos] /= total - 1; // normalize score
					
					// assume all non-proteogenic letters are gaps
					// penalize maximally for being half gaps
					fractionGaps = (double)(x_length - total)/x_length;
					fractionGaps -= 0.5;
					fractionGaps = 1 - 4*fractionGaps*fractionGaps;
					scoresSubset[pos] += GL*fractionGaps; // gap-letter mismatch
					
					if (gapsSubset[pos] > 0) { // add cost for any gap events
						weight = 0; // length of gap
						temp = gapsSubset[pos]; // gets zeroed-out
						while (temp > 0 && weight <= maxSize) {
							if ((temp & 1)==1) {
								scoresSubset[pos] += weight*GE + GO;
							}
							temp >>= 1;
							weight++;
						}
					}
				}
				
				double currentScore = 0;
				for (k = 0; k < subset; k++) {
					currentScore += scoresSubset[k];
				}
				
				currentScore += SC*j; // add shift penalty (offset cost)
				
				// if score improved then replace
				if ((currentScore - totScore) > bestScore) {
					bestScore = currentScore - totScore;
					bestShift = j;
					memcpy(basesSubsetRightSaved, basesSubset, (7*subset) * sizeof(double));
					memcpy(gapsSubsetRightSaved, gapsSubset, subset * sizeof(unsigned long long int));
					memcpy(scoresRightSaved, scoresSubset, subset * sizeof(double));
				}
			}
			
			Free(basesSubset);
			Free(gapsSubset);
			Free(scoresSubset);
			
			// commit the best shift if above threshold
			if (bestScore > threshold) {
				*(rans) = *(rans) + 1; // another change made
				if (bestShift < 0) { // left shift
					// apply shift to sequences
					bestShift *= -1;
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						// move left position to the right
						for (j = 1; j <= bestShift; j++) { // shift offset
							p = x_i.ptr[position[min] - j + 1]; // right position
							*((char *)x_i.ptr + position[min] - j + 1) = *((char *)x_i.ptr + position[min] - j - size + 1);
							*((char *)x_i.ptr + position[min] - j - size + 1) = p; // left position
						}
					}
					
					// overwrite original information with saved
					subset = left + size + 1;
					for (i = 0; i <= left + size; i++) { // shift offset
						end = position[min] - i;
						for (j = 0; j < 4; j++)
							bases[j*seqLength + end] = basesSubsetLeftSaved[j*subset + subset - i - 1];
						
						// adjust gapCount +/- 1
						if (gaps[end] != gapsSubsetLeftSaved[subset - i - 1]) {
							temp = gaps[end] ^ gapsSubsetLeftSaved[subset - i - 1];
							
							// ignore gap events <= size
							weight = size + 1; // size of gap event
							temp >>= weight;
							while (temp > 0 && weight <= maxSize) {
								if ((temp & 1)==1) {
									if (((gaps[end] >> weight) & 1)==1) { // originally there was a gap event
										gapCount[weight]--;
									} else { // there is a new gap event of weight size
										gapCount[weight]++;
									}
								}
								temp >>= 1;
								weight++;
							}
							
							gaps[end] = gapsSubsetLeftSaved[subset - i - 1];
						}
						scores[end] = scoresLeftSaved[subset - i - 1];
					}
				} else { // right shift
					// apply shift to sequences
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						// move right position to the left
						for (j = 1; j <= bestShift; j++) { // shift offset
							p = x_i.ptr[position[min] + j - size]; // left position
							*((char *)x_i.ptr + position[min] + j - size) = *((char *)x_i.ptr + position[min] + j);
							*((char *)x_i.ptr + position[min] + j) = p; // right position
						}
					}
					
					// overwrite original information with saved
					subset = right + size + maxSize;
					if ((position[min] - size + 1 + subset) > seqLength)
						subset = seqLength - (position[min] - size) - 1;
					for (i = 0; i < subset; i++) { // shift offset
						end = position[min] - size + i;
						for (j = 0; j < 4; j++)
							bases[j*seqLength + end] = basesSubsetRightSaved[j*subset + i];
						
						// adjust gapCount +/- 1
						if (gaps[end] != gapsSubsetRightSaved[i]) {
							temp = gaps[end] ^ gapsSubsetRightSaved[i];
							
							// ignore gap events <= size
							weight = size + 1; // size of gap event
							temp >>= weight;
							while (temp > 0 && weight <= maxSize) {
								if ((temp & 1)==1) {
									if (((gaps[end] >> weight) & 1)==1) { // originally there was a gap event
										gapCount[weight]--;
									} else { // there is a new gap event of weight size
										gapCount[weight]++;
									}
								}
								temp >>= 1;
								weight++;
							}
							
							gaps[end] = gapsSubsetRightSaved[i];
						}
						scores[end] = scoresRightSaved[i];
					}
				}
			}
			
			Free(basesSubsetLeftSaved);
			Free(gapsSubsetLeftSaved);
			Free(scoresLeftSaved);
			Free(basesSubsetRightSaved);
			Free(gapsSubsetRightSaved);
			Free(scoresRightSaved);
			
			// eliminate this gap event from further consideration
			for (j = min; j < (gapCount[size] - 1); j++) {
				gapNumber[j] = gapNumber[j + 1];
				position[j] = position[j + 1];
			}
			gapCount[size]--;
			Free(seqNumbers);
		}
		
		Free(gapNumber);
		Free(position);
	}
	
	Free(bases);
	Free(gaps);
	Free(gapLengths);
	Free(scores);
	Free(gapCount);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP shiftGapsAA(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP gl, SEXP sc, SEXP thresh, SEXP weights)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, k, seqLength, weight, tot, pos, bound, end;
	char p;
	unsigned long long int temp;
	int maxSize = sizeof(unsigned long long int)*8; // number of bits in an int
	double *subM = REAL(subMatrix);
	double *w = REAL(weights);
	double GO = asReal(go); // gap opening
	double GE = asReal(ge); // gap extension
	double GL = asReal(gl); // gap-letter mismatch
	double SC = asReal(sc); // shift cost
	double threshold = asReal(thresh); // increase in score to commit changes
	double fractionGaps, total, sum;
	int maxShift = maxSize - 1; // max distance to shift gap events to the left or right
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	seqLength = x_i.length; // sequences are assumed to be aligned (equal length)
	
	// return the number of gap-event shifts made to the sequences
	SEXP ans;
	double *rans;
	PROTECT(ans = allocVector(REALSXP, 1));
	rans = REAL(ans);
	*(rans) = 0; // initialize to zero changes
	
	// initialize an array of encoded base counts
	double *bases = Calloc(26*seqLength, double); // initialized to zero
	// initialize an array of encoded gap events
	unsigned long long int *gaps = Calloc(seqLength, unsigned long long int); // initialized to zero
	// initialize an array of terminal gap lengths
	int *gapLengths = Calloc(x_length*2, int); // initialized to zero
	// initialize an array of gap counts of each length
	int *gapCount = Calloc(maxSize + 1, int); // initialized to zero
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i*2] = frontTerminalGapsAA(&x_i);
		if (gapLengths[i*2]==x_i.length) // all gaps
			continue;
		gapLengths[i*2 + 1] = endTerminalGapsAA(&x_i);
		alphabetFrequencyAA(&x_i, bases, seqLength, 1, 0, gapLengths[i*2], gapLengths[i*2 + 1], w[i]);
		
		// bit-encode the length of gaps in an integer vector:
		// 1's for each length of gap closing at that position
		pos = 0;
		for (j = gapLengths[i*2] + 1; j < seqLength - gapLengths[i*2 + 1] - 1; j++) {
			if (pos==0 && // not inside a gap - check for opening
				(x_i.ptr[j - 1] ^ 0x2D && x_i.ptr[j - 1] ^ 0x2E) && (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E))) {
				pos = j; // gap opening
			}
			if (pos > 0 && // inside a gap - check for closing
				(x_i.ptr[j + 1] ^ 0x2D && x_i.ptr[j + 1] ^ 0x2E) && (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E))) {
				pos = j - pos + 1; // change in position
				if (pos > maxSize) // bigger gap than can be stored
					pos = maxSize;
				gaps[j] |= (unsigned long long int)1 << pos; // encode length of gap
				pos = 0;
			}
		}
	}
	
	// initialize an array of column scores
	double *scores = Calloc(seqLength, double); // initialized to zero
	
	// populate column scores and count gap events
	for (k = 0; k < seqLength; k++) {
		total = 0;
		for (i = 0; i < 20; i++) {
			total += bases[i*seqLength + k];
			for (j = i; j < 20; j++) {
				sum = (i==j) ? 1 : 2;
				sum *= (bases[j*seqLength + k] - ((i==j) ? 1 : 0));
				sum *= bases[i*seqLength + k];
				if (sum > 0)
					scores[k] += *(subM + i*21 + j)*sum;
			}
		}
		scores[k] /= x_length;
		if ((int)total > 1)
			scores[k] /= total - 1; // normalize score
		
		// assume all non-proteogenic letters are gaps
		// penalize maximally for being half gaps
		fractionGaps = (double)(x_length - total)/x_length;
		fractionGaps -= 0.5;
		fractionGaps = 1 - 4*fractionGaps*fractionGaps;
		scores[k] += GL*fractionGaps; // gap-letter mismatch
//		Rprintf("\nk %d total %d percent %1.2f adjusted %1.2f penalty %1.2f", k, total, (double)(x_length - total)/x_length, fractionGaps, GL*fractionGaps);
		
		if (gaps[k] > 0) { // add cost for any gap events
			weight = 0; // length of gap
			temp = gaps[k]; // gets zeroed-out
			while (temp > 0 && weight <= maxSize) {
				if ((temp & 1)==1) {
					scores[k] += weight*GE + GO;
					gapCount[weight]++;
				}
				temp >>= 1;
				weight++;
			}
		}
		
//		scores[k] *= total; // weight by number of bases
//		Rprintf("\nk %d score %f", k, scores[k]);
	}
	
//	for (k = 0; k <= maxSize; k++) {
//		Rprintf("\nk = %d gapCount = %d", k, gapCount[k]);
//	}
	
	// ALGORITHM:
	// for each gap size 1 -> maxSize
	//  continue if gapCount[size]==0 (no gaps)
	//  for each sequence
	//   record the number of sequences with a gap of that size at that position:
	//    for each position with a gap of that size
	//     if there is a gap in that position of that size
	//      increment the count for that position
	//  for each position with gaps of that size
	//   order by number of gaps increasing
	//   find all sequences to shift
	//   find the boundaries of the shift
	//    shift gap left and right up to some max distance
	//     stop when the end or another gap is reached
	//     record scores for each shift
	//   if a better score is available
	//    shift gap to that position in the sequence(s)
	//    if gaps were merged then increment gaps and gapCount
	//    if gaps remained the same length then do nothing
	//   gapCount[size]--
	//   break when gapCount[size]==0
	// return the number of changes made
	
	int count, correct, size, min, left, right;
	for (size = 1; size <= maxSize; size++) { // each gap size
		if (gapCount[size]==0) // no gaps of this size
			continue;
//		Rprintf("\n\nSize = %d", size);
		
		// initialize a vector of gap counts of each length
		int *gapNumber = Calloc(gapCount[size], int); // initialized to zero
		// initialize a vector of positions where gap events end (close)
		int *position = Calloc(gapCount[size], int); // initialized to zero
		count = 0;
		for (j = size; j < (seqLength - 1); j++) { // each position
			if ((gaps[j] & ((unsigned long long int)1 << size)) != 0) { // gap of correct size ending at this position
//				Rprintf("\nPosition = %d\n", j);
				for (i = 0; i < x_length; i++) { // each sequence
					x_i = get_elt_from_XStringSet_holder(&x_set, i);
					
					correct = 0;
					if ((x_i.ptr[j + 1] ^ 0x2D && x_i.ptr[j + 1] ^ 0x2E) && // next position is a letter
						(!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E)) && // this position is a gap
						(x_i.ptr[j - size] ^ 0x2D && x_i.ptr[j - size] ^ 0x2E)) { // past size-away is a letter
						correct = 1;
						for (k = size; k > 1; k--) {
							if (x_i.ptr[j - k + 1] ^ 0x2D && x_i.ptr[j - k + 1] ^ 0x2E) { // letter in this position
								correct = 0;
								break;
							}
						}
					}
					
//					if (correct)
//						Rprintf("Sequence = %d\n", i);
					if (correct)
						gapNumber[count]++;
				}
				
				position[count] = j;
				count++;
				
				if (count==gapCount[size]) // reached the expected number of gaps
					break;
			}
		}
//		if (count < gapCount[size])
//			error("ERROR:  MISSING GAPS!!!!!!!!!!!!");
		
		while (gapCount[size] > 0) {
			// find the gap event with the fewest sequences
			min = 0;
			for (j = 1; j < gapCount[size]; j++)
				if (gapNumber[j] < gapNumber[min])
					min = j;
			if (gapNumber[min]==0) // no gap events found
				break;
			
//			Rprintf("\n\nSize = %d Position = %d\n", size, position[min]);
			// find all sequences to with this gap event
//			int *seqNumbers = Calloc(gapNumber[min], int); // initialized to zero
			int *seqNumbers = Calloc(x_length, int); // initialized to zero
			j = position[min];
			count = 0;
			for (i = 0; i < x_length; i++) { // each sequence
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				
				correct = 0;
				if ((x_i.ptr[j + 1] ^ 0x2D && x_i.ptr[j + 1] ^ 0x2E) && // next position is a letter
					(!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E)) && // this position is a gap
					(x_i.ptr[j - size] ^ 0x2D && x_i.ptr[j - size] ^ 0x2E)) { // past size-away is a letter
					correct = 1;
					for (k = size; k > 1; k--) {
						if (x_i.ptr[j - k + 1] ^ 0x2D && x_i.ptr[j - k + 1] ^ 0x2E) { // letter in this position
							correct = 0;
							break;
						}
					}
				}
				
				if (correct) {
					seqNumbers[count] = i;
					count++;
//					Rprintf("Sequence = %d\n", i);
				}
			}
			gapNumber[min] = count; // reset the number of gaps to account for previous mergers
			
			// first and last position cannot be changed
			left = position[min] - size; // max possible shift left
			left = (left > maxShift) ? maxShift : left;
			right = seqLength - position[min] - 2; // max possible shift right
			right = (right > maxShift) ? maxShift : right;
			for (i = 0; i < gapNumber[min]; i++) { // each sequence
				x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
				
				// prevent a gap from shifting on top of another gap
				
				for (j = 1; j <= left; j++) { // shift offset
					pos = position[min] - size - j;
					if (!(x_i.ptr[pos] ^ 0x2D) || !(x_i.ptr[pos] ^ 0x2E)) { // position is a gap
						left = j; // shift this far left
						break;
					}
				}
//				Rprintf("Sequence = %d shift left = %d\n", i, left);
				
				for (j = 1; j <= right; j++) { // shift offset
					pos = position[min] + j + 1;
					if (!(x_i.ptr[pos] ^ 0x2D) || !(x_i.ptr[pos] ^ 0x2E)) { // position is a gap
						right = j; // shift this far right
						break;
					}
				}
//				Rprintf("Sequence = %d shift right = %d\n", i, right);
			}
			
			int bestShift = 0; // negative for left shift, positive for right shift
			double bestScore = 0; // change in score
			
			double totScore = 0; // score for this subset
			int subset = left + size + 1; // additional position is for the far left
			
			// initialize an array of encoded base counts
			double *basesSubset = Calloc(26*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			unsigned long long int *gapsSubset = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			double *scoresSubset = Calloc(subset, double); // initialized to zero
			// initialize an array of encoded base counts
			double *basesSubsetLeftSaved = Calloc(26*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			unsigned long long int *gapsSubsetLeftSaved = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			double *scoresLeftSaved = Calloc(subset, double); // initialized to zero
			
			// copy subset of counts from bases to basesSubset
			for (i = 0; i <= left + size; i++) { // shift offset
				end = position[min] - i;
				for (j = 0; j < 20; j++)
					basesSubset[j*subset + subset - i - 1] = bases[j*seqLength + end];
//				Rprintf("%d ", (int)scores[end]);
				gapsSubset[subset - i - 1] = gaps[end];
				scoresSubset[subset - i - 1] = scores[end];
				totScore += scores[end];
			}
//			Rprintf("\ntotScore = %f", totScore);
			
			// shift the gap event to the left one position at a time
			for (j = 1; j <= left; j++) { // shift offset
				// simulate swapping the left and right bases and gaps
				// right position is (position[min] - j + 1) in sequence and (subset - j) in basesSubset
				// left position is (position[min] - j - size + 1) in sequence and (subset - j - size) in basesSubset
				
				for (i = 0; i < gapNumber[min]; i++) { // each sequence
//					Rprintf("\nSeq %d", seqNumbers[i]);
					x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
					
//					// move right position to the left
//					// subtract from right position
//					adjustFrequencyAA(x_i.ptr[position[min] - j + 1], basesSubset, subset, 1, 0, subset - j, -1*w[seqNumbers[i]]);
//					// add to left position
//					adjustFrequencyAA(x_i.ptr[position[min] - j + 1], basesSubset, subset, 1, 0, subset - j - size, w[seqNumbers[i]]);
					
					// move left position to the right
					// subtract from left position
					adjustFrequencyAA(x_i.ptr[position[min] - j - size + 1], basesSubset, subset, 1, 0, subset - j - size, -1*w[seqNumbers[i]]);
					// add to right position
					adjustFrequencyAA(x_i.ptr[position[min] - j - size + 1], basesSubset, subset, 1, 0, subset - j, w[seqNumbers[i]]);
				}
//				//Rprintf("\nsubtracted right %d from subset %d", position[min] - j + 1, subset - j);
//				//Rprintf("\nadded right %d to subset %d", position[min] - j + 1, subset - j - size);
//				Rprintf("\nsubtracted left %d from subset %d", position[min] - j - size + 1, subset - j - size);
//				Rprintf("\nadded left %d to subset %d", position[min] - j - size + 1, subset - j);
				
				// re-encode gaps
				// past gap closing is (subset - j) in gapsSubset
				// new gap closing is (subset - j - 1) in gapsSubset
				
				if (j==1) { // clear the original gap of size
					gapsSubset[subset - j] &= ~((unsigned long long int)1 << size);
				} else { // restore the original gaps
					gapsSubset[subset - j] = gaps[position[min] - j + 1];
				}
				// add the new gap
//				Rprintf("\nchecked position = %d", subset - j - 2);
				if (j == left && // j is maximally left
					gapsSubset[subset - j - 1 - size]) { // there is an adjacent gap closing
					
					// initialize a vector of gap sizes to confirm in the adjacent position
					int *confirm = Calloc(maxSize, int); // initialized to zero
					
					// possibly merge gaps
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
//						Rprintf("here1");
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						bound = position[min] - j - size - maxSize;
						bound = (bound < 0) ? 0 : bound;
						for (k = position[min] - j - size; k >= bound; k--) {
//							Rprintf("here2");
							if (x_i.ptr[k] ^ 0x2D && x_i.ptr[k] ^ 0x2E) { // letter in this position
//								Rprintf("here3");
								pos = position[min] - j - k; // start - finish - 1 + size
//								Rprintf("\nposition = %d", pos); // length of gap
								if (pos - size < maxSize) // need to confirm gap of this size
									confirm[pos - size] = 1;
								if (pos > maxSize)
									pos = maxSize;
								gapsSubset[subset - j - 1] |= ((unsigned long long int)1 << pos);
								break;
							}
						}
					}
					
					int countConfirm = 0;
					for (i = 1; i < maxSize - 1; i++) {
						if (confirm[i] > 0) {
							confirm[countConfirm] = i; // refill with gap lengths to confirm
							countConfirm++;
						}
					}
					
					if (countConfirm > 0) {
						count = 0; // tracks with sequence in seqNumbers
						for (i = 0; i < x_length; i++) {
							if (count < gapNumber[min] &&
								i == seqNumbers[count]) {
								// do not check shifted sequences
								count++;
								continue;
							}
							
							x_i = get_elt_from_XStringSet_holder(&x_set, i);
							
							pos = position[min] - j - size + 1;
//							Rprintf("\nSequence %d checking gap closing at %d", i, pos);
							if (!(x_i.ptr[pos] ^ 0x2D) || !(x_i.ptr[pos] ^ 0x2E)) // gap in this position
								continue; // not a gap closing event
							
							int confirmed = 0; // the number of left gaps in this sequence
							bound = position[min] - j - size - maxSize;
							bound = (bound < 0) ? 0 : bound;
							for (pos = position[min] - j - size; pos >= bound; pos--) {
								if (x_i.ptr[pos] ^ 0x2D && x_i.ptr[pos] ^ 0x2E) { // letter in this position
//									Rprintf("\nSequence %d confirmed gap of length %d", i, confirmed);
									break;
								} else {
									confirmed++;
								}
							}
							
							if (confirmed) { // at least one gap ending at the adjacent gap position
								for (k = 0; k < countConfirm; k++) {
									if (confirm[k]==confirmed) { // this gap size was confirmed
										for (pos = k; pos < countConfirm; pos++)
											confirm[pos] = confirm[pos + 1]; // shift confirm left
										countConfirm--;
										break;
									}
								}
							}
							
							if (countConfirm==0) // no more gap events to confirm
								break;
						}
						
						// clear gaps that were not confirmed
						for (k = 0; k < countConfirm; k++)
							gapsSubset[subset - j - 1 - size] &= ~((unsigned long long int)1 << confirm[k]);
					}
					
					Free(confirm);
				} else {
					gapsSubset[subset - j - 1] |= ((unsigned long long int)1 << size);
				}
				
//				Rprintf("\nshift %d\n", j);
//				for (i = 0; i < subset; i++)
//					Rprintf("%d ", (int)gapsSubset[i]);
				
				// re-compute score in the changed positions
				for (pos = 0; pos < subset; pos++) {
					if (!(pos==(subset - j) ||
						pos==(subset - j - 1) ||
						pos==(subset - j - size) ||
						pos==(subset - j - 1 - size)))
						continue; // nothing changed in this position
					
					scoresSubset[pos] = 0; // reset the score in this column
					total = 0;
//					Rprintf("\npos %d: ", pos);
					for (i = 0; i < 20; i++) {
//						Rprintf("%1.2f ", basesSubset[i*subset + pos]);
						total += basesSubset[i*subset + pos];
						for (k = i; k < 20; k++) {
							sum = (i==k) ? 1 : 2;
							sum *= (basesSubset[k*subset + pos] - ((i==k) ? 1 : 0));
							sum *= basesSubset[i*subset + pos];
							if (sum > 0)
								scoresSubset[pos] += *(subM + i*21 + k)*sum;
						}
					}
					scoresSubset[pos] /= x_length;
					if ((int)total > 1)
						scoresSubset[pos] /= total - 1; // normalize score
					
					// assume all non-proteogenic letters are gaps
					// penalize maximally for being half gaps
					fractionGaps = (double)(x_length - total)/x_length;
					fractionGaps -= 0.5;
					fractionGaps = 1 - 4*fractionGaps*fractionGaps;
					scoresSubset[pos] += GL*fractionGaps; // gap-letter mismatch
					
					if (gapsSubset[pos] > 0) { // add cost for any gap events
						weight = 0; // length of gap
						temp = gapsSubset[pos]; // gets zeroed-out
						while (temp > 0 && weight <= maxSize) {
							if ((temp & 1)==1) {
								scoresSubset[pos] += weight*GE + GO;
							}
							temp >>= 1;
							weight++;
						}
					}
					
//					scoresSubset[pos] *= total; // weight by number of bases
				}
				
				double currentScore = 0;
				for (k = 0; k < subset; k++) {
					currentScore += scoresSubset[k];
//					Rprintf("%d ", (int)scoresSubset[k]);
				}
//				Rprintf("\nshift %d currentScore %f", -1*j, currentScore);
				
				currentScore += SC*j; // add shift penalty (offset cost)
				
				// if score improved then replace
				if ((currentScore - totScore) > bestScore) {
//					Rprintf("\nshift %d score %f", -1*j, currentScore - totScore);
					bestScore = currentScore - totScore;
					bestShift = -1*j;
					memcpy(basesSubsetLeftSaved, basesSubset, (26*subset) * sizeof(double));
					memcpy(gapsSubsetLeftSaved, gapsSubset, subset * sizeof(unsigned long long int));
					memcpy(scoresLeftSaved, scoresSubset, subset * sizeof(double));
				}
			}
			
			Free(basesSubset);
			Free(gapsSubset);
			Free(scoresSubset);
			
			// repeat process for right shifts
			
			totScore = 0; // score for this subset
			// subset covers from just left of the gap (position[min] - size) to (position[min] + right + maxSize)
			subset = right + size + maxSize; // additional maxSize positions are for the right side
			if ((position[min] - size + 1 + subset) > seqLength)
				subset = seqLength - (position[min] - size) - 1;
			
			// initialize an array of encoded base counts
			basesSubset = Calloc(26*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			gapsSubset = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			scoresSubset = Calloc(subset, double); // initialized to zero
			// initialize an array of encoded base counts
			double *basesSubsetRightSaved = Calloc(26*subset, double); // initialized to zero
			// initialize an array of encoded gap events
			unsigned long long int *gapsSubsetRightSaved = Calloc(subset, unsigned long long int); // initialized to zero
			// initialize a vector of column scores
			double *scoresRightSaved = Calloc(subset, double); // initialized to zero
			
			// copy subset of counts from bases to basesSubset
			for (i = 0; i < subset; i++) {
				end = position[min] - size + i;
				for (j = 0; j < 20; j++)
					basesSubset[j*subset + i] = bases[j*seqLength + end];
//				Rprintf("%d ", (int)scores[end]);
				gapsSubset[i] = gaps[end];
				scoresSubset[i] = scores[end];
				totScore += scores[end];
			}
//			Rprintf("\ntotScore = %f", totScore);
			
			// shift the gap event to the right one position at a time
			for (j = 1; j <= right; j++) { // shift offset
				// simulate swapping the left and right bases and gaps
				
				for (i = 0; i < gapNumber[min]; i++) { // each sequence
//					Rprintf("\nSeq %d", seqNumbers[i]);
					x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
					
					// move left position to the right
					// subtract from right position
					adjustFrequencyAA(x_i.ptr[position[min] + j], basesSubset, subset, 1, 0, j + size, -1*w[seqNumbers[i]]);
					// add to left position
					adjustFrequencyAA(x_i.ptr[position[min] + j], basesSubset, subset, 1, 0, j, w[seqNumbers[i]]);
				}
//				Rprintf("\nsubtracted right %d from subset %d", position[min] + j, j + size);
//				Rprintf("\nadded right %d to subset %d", position[min] + j, j);
				
				// re-encode gaps
				// past gap closing is (j - 1 + size) in gapsSubset
				// new gap closing is (j + size) in gapsSubset
				
				if (j==1) { // clear the original gap of size
					gapsSubset[j - 1 + size] &= ~((unsigned long long int)1 << size);
				} else { // restore the original gaps
					gapsSubset[j - 1 + size] = gaps[position[min] + j - 1];
				}
				// add the new gap
				if (j == right) { // j is maximally right
					// check for an adjacent gap opening
					
					// initialize a vector of gap sizes to confirm in the adjacent position
					int *confirm = Calloc(maxSize, int); // initialized to zero
					
					// possibly merge gaps
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
//						Rprintf("here1");
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						bound = position[min] + j + maxSize;
						bound = (bound > seqLength) ? seqLength : bound;
						for (k = position[min] + j + 1; k < bound; k++) {
//							Rprintf("here2");
							if (x_i.ptr[k] ^ 0x2D && x_i.ptr[k] ^ 0x2E) { // letter in this position
//								Rprintf("here3");
								pos = k - (position[min] + j + 1) + size; // finish - start + size
//								Rprintf("\nposition = %d", pos); // length of gap
								if (pos - size < maxSize) // need to confirm gap of this size
									confirm[pos - size] = 1;
								if (pos > maxSize)
									pos = maxSize;
								gapsSubset[j + pos] |= ((unsigned long long int)1 << pos);
								break;
							}
						}
					}
					
					int countConfirm = 0;
					for (i = 1; i < maxSize - 1; i++) {
						if (confirm[i] > 0) {
							confirm[countConfirm] = i; // refill with gap lengths to confirm
							countConfirm++;
						}
					}
					
					if (countConfirm > 0) {
						count = 0; // tracks with sequence in seqNumbers
						for (i = 0; i < x_length; i++) {
							if (count < gapNumber[min] &&
								i == seqNumbers[count]) {
								// do not check shifted sequences
								count++;
								continue;
							}
							
							x_i = get_elt_from_XStringSet_holder(&x_set, i);
							
							pos = position[min] + j;
//							Rprintf("\nSequence %d checking gap closing at %d", i, pos);
							if (!((x_i.ptr[pos] ^ 0x2D && x_i.ptr[pos] ^ 0x2E) && // letter in this position
								(!(x_i.ptr[pos + 1] ^ 0x2D) || !(x_i.ptr[pos + 1] ^ 0x2E)))) // gap in this position
								continue; // not a gap closing event
							
							int confirmed = 1; // the number of right gaps in this sequence
							bound = position[min] + j + maxSize;
							bound = (bound > seqLength) ? seqLength : bound;
							for (pos = position[min] + j + 2; pos < bound; pos++) {
								if (x_i.ptr[pos] ^ 0x2D && x_i.ptr[pos] ^ 0x2E) { // letter in this position
//									Rprintf("\nSequence %d confirmed gap of length %d", i, confirmed);
									break;
								} else {
									confirmed++;
								}
							}
							
							if (confirmed) { // at least one gap ending at the adjacent gap position
								for (k = 0; k < countConfirm; k++) {
									if (confirm[k]==confirmed) { // this gap size was confirmed
										for (pos = k; pos < countConfirm; pos++)
											confirm[pos] = confirm[pos + 1]; // shift confirm left
										countConfirm--;
										break;
									}
								}
							}
							
							if (countConfirm==0) // no more gap events to confirm
								break;
						}
						
						// clear gaps that were not confirmed
						for (k = 0; k < countConfirm; k++)
							gapsSubset[j + confirm[k] + size] &= ~((unsigned long long int)1 << confirm[k]);
					}
					
					Free(confirm);
				} else {
					gapsSubset[j + size] |= ((unsigned long long int)1 << size);
				}
				
//				Rprintf("\nshift %d\n", j);
//				for (i = 0; i < subset; i++)
//					Rprintf("%d ", (int)gapsSubset[i]);
				
				// re-compute score in the changed positions
				for (pos = 0; pos < subset; pos++) {
					if (!((j == right && pos > j) ||
						pos==j ||
						pos==(j + size) ||
						pos==(j - 1 + size)))
						continue; // nothing changed in this position
					
					scoresSubset[pos] = 0; // reset the score in this column
					total = 0;
//					Rprintf("\npos %d: ", pos);
					for (i = 0; i < 20; i++) {
//						Rprintf("%1.2f ", basesSubset[i*subset + pos]);
						total += basesSubset[i*subset + pos];
						for (k = i; k < 20; k++) {
							sum = (i==k) ? 1 : 2;
							sum *= (basesSubset[k*subset + pos] - ((i==k) ? 1 : 0));
							sum *= basesSubset[i*subset + pos];
							if (sum > 0)
								scoresSubset[pos] += *(subM + i*21 + k)*sum;
						}
					}
					scoresSubset[pos] /= x_length;
					if ((int)total > 1)
						scoresSubset[pos] /= total - 1; // normalize score
					
					// assume all non-proteogenic letters are gaps
					// penalize maximally for being half gaps
					fractionGaps = (double)(x_length - total)/x_length;
					fractionGaps -= 0.5;
					fractionGaps = 1 - 4*fractionGaps*fractionGaps;
					scoresSubset[pos] += GL*fractionGaps; // gap-letter mismatch
					
					if (gapsSubset[pos] > 0) { // add cost for any gap events
						weight = 0; // length of gap
						temp = gapsSubset[pos]; // gets zeroed-out
						while (temp > 0 && weight <= maxSize) {
							if ((temp & 1)==1) {
								scoresSubset[pos] += weight*GE + GO;
							}
							temp >>= 1;
							weight++;
						}
					}
					
//					scoresSubset[pos] *= total; // weight by number of bases
				}
				
				double currentScore = 0;
				for (k = 0; k < subset; k++) {
					currentScore += scoresSubset[k];
//					Rprintf("%d ", (int)scoresSubset[k]);
				}
//				Rprintf("\nshift %d currentScore %f", j, currentScore);
				
				currentScore += SC*j; // add shift penalty (offset cost)
				
				// if score improved then replace
				if ((currentScore - totScore) > bestScore) {
//					Rprintf("\nshift %d score %f", j, currentScore - totScore);
					bestScore = currentScore - totScore;
					bestShift = j;
					memcpy(basesSubsetRightSaved, basesSubset, (26*subset) * sizeof(double));
					memcpy(gapsSubsetRightSaved, gapsSubset, subset * sizeof(unsigned long long int));
					memcpy(scoresRightSaved, scoresSubset, subset * sizeof(double));
				}
			}
			
			Free(basesSubset);
			Free(gapsSubset);
			Free(scoresSubset);
			
			// commit the best shift if above threshold
			if (bestScore > threshold) {
				*(rans) = *(rans) + 1; // another change made
				if (bestShift < 0) { // left shift
					// apply shift to sequences
					bestShift *= -1;
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						// move left position to the right
						for (j = 1; j <= bestShift; j++) { // shift offset
							p = x_i.ptr[position[min] - j + 1]; // right position
							*((char *)x_i.ptr + position[min] - j + 1) = *((char *)x_i.ptr + position[min] - j - size + 1);
							*((char *)x_i.ptr + position[min] - j - size + 1) = p; // left position
//							Rprintf("\nSequence %d from pos %d to pos %d", seqNumbers[i], position[min] - j - size + 1, position[min] - j + 1);
						}
					}
					
					// overwrite original information with saved
					subset = left + size + 1;
					for (i = 0; i <= left + size; i++) { // shift offset
						end = position[min] - i;
						for (j = 0; j < 20; j++)
							bases[j*seqLength + end] = basesSubsetLeftSaved[j*subset + subset - i - 1];
						
						// adjust gapCount +/- 1
						if (gaps[end] != gapsSubsetLeftSaved[subset - i - 1]) {
							temp = gaps[end] ^ gapsSubsetLeftSaved[subset - i - 1];
//							Rprintf("\n%d %d %d", gaps[end], gapsSubsetLeftSaved[subset - i - 1], temp);
							
							// ignore gap events <= size
							weight = size + 1; // size of gap event
							temp >>= weight;
							while (temp > 0 && weight <= maxSize) {
								if ((temp & 1)==1) {
									if (((gaps[end] >> weight) & 1)==1) { // originally there was a gap event
										gapCount[weight]--;
//										Rprintf("\nGap of size %d lost in position %d", weight, end);
									} else { // there is a new gap event of weight size
										gapCount[weight]++;
//										Rprintf("\nGap of size %d gained in position %d", weight, end);
									}
								}
								temp >>= 1;
								weight++;
							}
							
							gaps[end] = gapsSubsetLeftSaved[subset - i - 1];
						}
						scores[end] = scoresLeftSaved[subset - i - 1];
					}
				} else { // right shift
					// apply shift to sequences
					for (i = 0; i < gapNumber[min]; i++) { // each sequence
						x_i = get_elt_from_XStringSet_holder(&x_set, seqNumbers[i]);
						
						// move right position to the left
						for (j = 1; j <= bestShift; j++) { // shift offset
							p = x_i.ptr[position[min] + j - size]; // left position
							*((char *)x_i.ptr + position[min] + j - size) = *((char *)x_i.ptr + position[min] + j);
							*((char *)x_i.ptr + position[min] + j) = p; // right position
//							Rprintf("\nSequence %d from pos %d to pos %d", seqNumbers[i], position[min] + j, position[min] + j - size);
						}
					}
					
					// overwrite original information with saved
					subset = right + size + maxSize;
					if ((position[min] - size + 1 + subset) > seqLength)
						subset = seqLength - (position[min] - size) - 1;
					for (i = 0; i < subset; i++) { // shift offset
						end = position[min] - size + i;
						for (j = 0; j < 20; j++)
							bases[j*seqLength + end] = basesSubsetRightSaved[j*subset + i];
						
						// adjust gapCount +/- 1
						if (gaps[end] != gapsSubsetRightSaved[i]) {
							temp = gaps[end] ^ gapsSubsetRightSaved[i];
//							Rprintf("\n%d %d %d", gaps[end], gapsSubsetRightSaved[i], temp);
							
							// ignore gap events <= size
							weight = size + 1; // size of gap event
							temp >>= weight;
							while (temp > 0 && weight <= maxSize) {
								if ((temp & 1)==1) {
									if (((gaps[end] >> weight) & 1)==1) { // originally there was a gap event
										gapCount[weight]--;
//										Rprintf("\nGap of size %d lost in position %d", weight, end);
									} else { // there is a new gap event of weight size
										gapCount[weight]++;
//										Rprintf("\nGap of size %d gained in position %d", weight, end);
									}
								}
								temp >>= 1;
								weight++;
							}
							
							gaps[end] = gapsSubsetRightSaved[i];
						}
						scores[end] = scoresRightSaved[i];
					}
				}
			}
			
			Free(basesSubsetLeftSaved);
			Free(gapsSubsetLeftSaved);
			Free(scoresLeftSaved);
			Free(basesSubsetRightSaved);
			Free(gapsSubsetRightSaved);
			Free(scoresRightSaved);
			
			// eliminate this gap event from further consideration
			for (j = min; j < (gapCount[size] - 1); j++) {
				gapNumber[j] = gapNumber[j + 1];
				position[j] = position[j + 1];
			}
			gapCount[size]--;
			Free(seqNumbers);
		}
		
		Free(gapNumber);
		Free(position);
	}
	
	Free(bases);
	Free(gaps);
	Free(gapLengths);
	Free(scores);
	Free(gapCount);
	
	UNPROTECT(1);
	
	return ans;
}
