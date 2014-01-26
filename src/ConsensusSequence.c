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

static int frontTerminalGaps(const Chars_holder *P)
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

static int endTerminalGaps(const Chars_holder *P)
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

static int frontTerminalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->seq;
	     i < P->length;
	     i++, p++)
	{
		if (!((*p) ^ 0x2D)) { // gap character
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
	for (i = (P->length - 1), p = (P->seq + P->length - 1);
	     i >= 0;
	     i--, p--)
	{
		if (!((*p) ^ 0x2D)) { // gap character
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

static void alphabetFrequencyAA(const Chars_holder *P, double *bits, int seqLength, int degeneracy, int ignore, int start, int end, double weight)
{
	int j, i;
	const char *p;
	
	for (j = start, p = (P->seq + start);
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
			default:
				error("not AA!");
				break;
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

static void makeConsensusAA(double *bits, char *seq, int seqLength, int x_length, double threshold, double minInfo, int tGaps)
{
	int j, i, X;
	double information, AAs[22], S, percentGap, percentMask;
	
	for (j = 0; j < seqLength; j++) {
		
		if (!tGaps && (*(bits + 25*seqLength + j)==0)) {
			*(seq + j) = '-';
			continue;
		}
		
		// find the percentage of each base in position j
		S = 0;
		for (i = 0; i < 22; i++) {
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
	int x_length, i, j, seqLength, degeneracy, ignore, tGaps;
	double *thresh, *minInfo;
	SEXP consensusSeq;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	int gapLengths[x_length][2];
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
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
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

//ans_start <- .Call("consensusSequenceAA", myDNAStringSet, threshold, ambiguity, minInformation, ignoreNonLetters, terminalGaps, PACKAGE="DECIPHER")
SEXP consensusSequenceAA(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, seqLength, degeneracy, ignore, tGaps;
	double *thresh, *minInfo;
	SEXP consensusSeq;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	int gapLengths[x_length][2];
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
	// initialize an array of encoded base counts
	double bases[26][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 26; i++) {
		for (j = 0; j < seqLength; j++) {
			bases[i][j] = 0;
		}
	}
	
	// loop through each sequence in the DNAStringSet
	for (i = 0; i < x_length; i++) {
		// extract each ith DNAString from the DNAStringSet
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		if (!tGaps) { // don't include terminal gaps
			gapLengths[i][0] = frontTerminalGapsAA(&x_i);
			gapLengths[i][1] = endTerminalGapsAA(&x_i);
			alphabetFrequencyAA(&x_i, &bases[0][0], seqLength, degeneracy, ignore, gapLengths[i][0], gapLengths[i][1], 1);
		} else { // include terminal gaps
			alphabetFrequencyAA(&x_i, &bases[0][0], seqLength, degeneracy, ignore, 0, 0, 1);
		}
	}
	
	thresh = REAL(threshold);
	*thresh = 1 - *thresh;
	minInfo = REAL(minInformation);
	char seq[seqLength + 1]; // last position is for null terminating
	makeConsensusAA(&bases[0][0], &seq[0], seqLength, x_length, *thresh, *minInfo, tGaps);
	seq[seqLength] = '\0'; // end (null terminate) the string
	
	PROTECT(consensusSeq = allocVector(STRSXP, 1));
	SET_STRING_ELT(consensusSeq, 0, mkChar(seq));
	
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
	double *rans, *w = REAL(weight), sum;//, *m
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	double gapLengths[x_length][2];
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
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
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i][0] = frontTerminalGaps(&x_i);
		gapLengths[i][1] = endTerminalGaps(&x_i);
		alphabetFrequency(&x_i, &bases[0][0], seqLength, 1, 0, gapLengths[i][0], gapLengths[i][1], w[i]);
		
		for (j = gapLengths[i][0] + 1; j < seqLength - gapLengths[i][1] - 1; j++) {
			if (!(x_i.seq[j - 1] & 0x10) && (x_i.seq[j] & 0x10)) {
				gaps[0][j] += w[i]; // gap opening
			}
			if ((x_i.seq[j] & 0x10) && !(x_i.seq[j + 1] & 0x10)) {
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
		sum = *(rans + i*7 + 0) + *(rans + i*7 + 1) + *(rans + i*7 + 2) + *(rans + i*7 + 3) + *(rans + i*7 + 4);
		if (sum > 0) { // normalize the profile
			*(rans + i*7 + 0) /= sum;
			*(rans + i*7 + 1) /= sum;
			*(rans + i*7 + 2) /= sum;
			*(rans + i*7 + 3) /= sum;
			*(rans + i*7 + 4) /= sum;
		}
	}
	
	//PROTECT(ret_list = allocVector(VECSXP, 2));
	//SET_VECTOR_ELT(ret_list, 0, ans);
	//SET_VECTOR_ELT(ret_list, 1, subM);
	
	UNPROTECT(1);
	
	return ans;
}

//ans_start <- .Call("consensusProfileAA", myAAStringSet, weight, PACKAGE="DECIPHER")
SEXP consensusProfileAA(SEXP x, SEXP weight)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i, j, seqLength;
	SEXP ans;
	double *rans, *w = REAL(weight), sum;
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	double gapLengths[x_length][2];
	
	// find the longest length XString
	seqLength = 0;
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (x_i.length > seqLength) {
			seqLength = x_i.length;
		}
	}
	// initialize an array of encoded base counts
	double bases[26][seqLength];
	// 2D arrays cannot be set to zero at initialization so loops are needed
	for (i = 0; i < 26; i++) {
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
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// update the alphabet for this string
		gapLengths[i][0] = frontTerminalGapsAA(&x_i);
		gapLengths[i][1] = endTerminalGapsAA(&x_i);
		// for AA degeneracy = 0, otherwise results in pooly aligned X's
		alphabetFrequencyAA(&x_i, &bases[0][0], seqLength, 0, 0, gapLengths[i][0], gapLengths[i][1], w[i]);
		
		for (j = gapLengths[i][0] + 1; j < seqLength - gapLengths[i][1] - 1; j++) {
			if ((x_i.seq[j - 1] ^ 0x2D) && !(x_i.seq[j] ^ 0x2D)) {
				gaps[0][j] += w[i]; // gap opening
			}
			if (!(x_i.seq[j] ^ 0x2D) && (x_i.seq[j + 1] ^ 0x2D)) {
				gaps[1][j] += w[i]; // gap closing
			}
		}
	}
	
	PROTECT(ans = allocMatrix(REALSXP, 26, seqLength));
	rans = REAL(ans);
	
	for (i = 0; i < seqLength; i++) {
		for (j = 0; j < 24; j++) {
			*(rans + i*26 + j) = bases[j][i]/x_length;
		}
		sum = 0;
		for (j = 0; j < 24; j++) {
			sum += *(rans + i*26 + j);
		}
		if (sum > 0) {
			for (j = 0; j < 24; j++) {
				*(rans + i*26 + j) /= sum;
			}
		}
		*(rans + i*26 + 24) = gaps[0][i]/x_length;
		*(rans + i*26 + 25) = gaps[1][i]/x_length;
	}
	
	UNPROTECT(1);
	
	return ans;
}
