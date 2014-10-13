/****************************************************************************
 *                            Designs Probe Set                             *
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

SEXP designProbes(SEXP x, SEXP max_pl, SEXP min_pl, SEXP max_c, SEXP numMMs, SEXP numPs, SEXP st, SEXP en, SEXP max_ov, SEXP h_percent, SEXP min_f, SEXP max_f, SEXP minS, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	XStringSet_holder x_set;
	Chars_holder x_i, x_l;
	int x_length; //seqLength;
	int before, v, *rPercentComplete;
	int max_probeLength = asInteger(max_pl);
	int min_probeLength = asInteger(min_pl);
	int maxCombos = asInteger(max_c);
	double hyb_percent = asReal(h_percent); // % formamide at hybridization
	double min_fam = asReal(min_f); // % minimum melting formamide
	double max_fam = asReal(max_f); // % maximum melting formamide
	double dGini = 1.96 + hyb_percent*0.17314; // dG initial + m1 * FAm
	int nMMs = asInteger(numMMs); // number of top MM probes to record (> 0)
	int i, j, k, l, p, n, o;
	int polyT = 20; // maximum length of poly-T base
	double *rans, *rMMs;
	SEXP ans, MMs, probes, perms, percentComplete, utilsPackage;
	int nthreads = asInteger(nThreads);
	
	v = asLogical(verbose);
	if (v) { // percent complete variables
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	// NN lookup table
	double NN[4][4] = {
		-0.816507461,-2.5401714,-1.647430026,-1.184658548
		,-1.854740485,-2.479102613,-2.826248182,-1.647430026
		,-2.48761723,-4.694133177,-2.479102613,-2.5401714
		,-0.495794417,-2.48761723,-1.854740485,-0.816507461
	};
	
	double PM[4][4] = {
		-0.141370102,-0.439805276,-0.285236035,-0.205111781
		,-0.321129768,-0.429231826,-0.48933661,-0.285236035
		,-0.430706047,-0.812742218,-0.429231826,-0.439805276
		,-0.085841845,-0.430706047,-0.321129768,-0.141370102
	};
	
	double sMM[4][5][5][4] = {
		0,0,0,0
		,1.545032445,1.254355018,1.491691514,1.329138183
		,1.150635633,0.582415494,1.075877275,1.187937642
		,1.203555051,1.001540513,0.864287715,0.717125848
		,0.75,0.65,0.69,0.78
		,0.630005348,0.18553379,0.730763505,0.709272397
		,0,0,0,0
		,0.856582783,-0.143236405,0.716721488,0.603652831
		,0.851622883,0.653168672,0.676545316,1.187937642
		,0.75,0.65,0.69,0.78
		,1.231861002,0.746214538,1.087821916,0.989140748
		,1.822113278,1.270687029,1.336192565,1.364584949
		,0,0,0,0
		,1.443665704,1.385046493,1.256013166,1.329138183
		,0.75,0.65,0.69,0.78
		,1.478009492,0.882097231,1.20450984,1.061002478
		,1.496720812,0.846496194,0.967868114,0.989140748
		,0.766581547,-0.024857805,0.50754303,0.709272397
		,0,0,0,0
		,0.75,0.65,0.69,0.78
		,0.75,0.65,0.69,0.78
		,0.75,0.65,0.69,0.78
		,0.75,0.65,0.69,0.78
		,0.76,0.65,0.69,0.78
		,0,0,0,0
		,0,0,0,0
		,1.295827995,0.84547091,0.91019099,1.256013166
		,0.755889609,0.241428373,0.396379912,0.676545316
		,0.99945386,0.740323132,0.435659206,0.864287715
		,0.65,0.55,0.48,0.69
		,0.843147406,0.101248351,0.49063599,0.50754303
		,0,0,0,0
		,1.0651638,0.249934344,0.699352949,0.716721488
		,0.871921533,0.59458138,0.396379912,1.075877275
		,0.65,0.56,0.49,0.69
		,1.07531714,0.318907854,0.653287717,0.967868114
		,1.099899195,0.730184613,0.661798984,1.336192565
		,0,0,0,0
		,1.45897431,1.318532145,0.91019099,1.491691514
		,0.65,0.56,0.49,0.69
		,1.242135174,0.894838095,1.108555445,1.20450984
		,0.911428974,0.524430101,0.653287717,1.087821916
		,0.503209827,0.274849491,0.49063599,0.730763505
		,0,0,0,0
		,0.65,0.55,0.48,0.69
		,0.65,0.55,0.48,0.69
		,0.65,0.56,0.49,0.69
		,0.65,0.56,0.49,0.69
		,0.65,0.55,0.48,0.69
		,0,0,0,0
		,0,0,0,0
		,1.100661785,0.969784756,1.318532145,1.385046493
		,0.565895968,-0.060347902,0.59458138,0.653168672
		,0.782168488,0.788161238,0.740323132,1.001540513
		,0.68,0.46,0.55,0.65
		,0.468913405,-0.469855984,0.274849491,-0.024857805
		,0,0,0,0
		,0.258195131,-0.70438632,0.249934344,-0.143236405
		,0.502914193,-0.060347902,0.241428373,0.582415494
		,0.68,0.47,0.56,0.65
		,0.584083861,0.258975454,0.524430101,0.846496194
		,0.968040559,0.797499702,0.730184613,1.270687029
		,0,0,0,0
		,1.081040749,0.969784756,0.84547091,1.254355018
		,0.68,0.47,0.56,0.65
		,1.048553951,0.728354541,0.894838095,0.882097231
		,0.88611252,0.258975454,0.318907854,0.746214538
		,0.239520858,-0.469855984,0.101248351,0.18553379
		,0,0,0,0
		,0.68,0.46,0.55,0.65
		,0.68,0.46,0.55,0.65
		,0.68,0.47,0.56,0.65
		,0.68,0.47,0.56,0.65
		,0.68,0.46,0.55,0.65
		,0,0,0,0
		,0,0,0,0
		,1.566899704,1.081040749,1.45897431,1.443665704
		,0.976725675,0.502914193,0.871921533,0.851622883
		,1.482046826,0.782168488,0.99945386,1.203555051
		,0.85,0.68,0.65,0.76
		,0.798628781,0.239520858,0.503209827,0.766581547
		,0,0,0,0
		,1.141098246,0.258195131,1.0651638,0.856582783
		,0.976725675,0.565895968,0.755889609,1.150635633
		,0.85,0.68,0.65,0.75
		,1.125403302,0.88611252,0.911428974,1.496720812
		,1.68169282,0.968040559,1.099899195,1.822113278
		,0,0,0,0
		,1.566899704,1.100661785,1.295827995,1.545032445
		,0.85,0.68,0.65,0.75
		,1.35948517,1.048553951,1.242135174,1.478009492
		,1.125403302,0.584083861,1.07531714,1.231861002
		,0.798628781,0.468913405,0.843147406,0.630005348
		,0,0,0,0
		,0.85,0.68,0.65,0.75
		,0.85,0.68,0.65,0.75
		,0.85,0.68,0.65,0.75
		,0.85,0.68,0.65,0.75
		,0.85,0.68,0.65,0.75
		,0,0,0,0	
	};
	
	// initialize the XStringSet
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	int numProbes = asInteger(numPs);
	PROTECT(probes = allocVector(STRSXP, x_length*numProbes)); // each probe
	PROTECT(perms = allocVector(STRSXP, x_length*numProbes)); // all permuations
	int cols = 9;
	PROTECT(ans = allocMatrix(REALSXP, x_length*numProbes, cols)); // probe's info
	rans = REAL(ans); // OTU, Start, Length, FAm, Score
	PROTECT(MMs = allocMatrix(REALSXP, x_length*numProbes, nMMs*2)); // worst mismatches
	rMMs = REAL(MMs); // first half of columns are hyb_effs, others are MM_seqs
	double mScore = asReal(minS);
	double minScore;
	// initialize to minimum score
	for (n = 0; n < ((x_length - 1)*numProbes + numProbes - 1 + (cols - 1)*x_length*numProbes + 1); n++)
		rans[n] = mScore;
	// initialize mismatches matrix
	for (n = 0; n < ((x_length - 1)*numProbes + numProbes - 1 + (2*nMMs - 1)*x_length*numProbes + 1); n++)
		rMMs[n] = -1;
	
	//seqLength = get_elt_from_XStringSet_holder(&x_set, 0).length;
	int start = asInteger(st) - 1; // first position AFTER forward primer
	int end = asInteger(en) - 1;//seqLength; // last position BEFORE reverse primer
	int maxOverlap = asInteger(max_ov);
	
	// loop through each sequence in the DNAStringSet
	for (i = 0; i < x_length; i++) {
		minScore = mScore;
		
		// extract each ith DNAString from the DNAStringSet
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		// make an array of position, ambiguity
		int pos[(end - start)][2];
		int count = 0;
		// determine the non-gap pos and ambiguities
		for (j = start; j < end; j++) {
			if (!(x_i.seq[j] & 0x10 || x_i.seq[j] & 0x40)) {
				pos[count][0] = j;
				switch (x_i.seq[j]) {
					case 3:
					case 5:
					case 6:
					case 9:
					case 10:
					case 12:
						pos[count][1] = 2;
						break;
					case 7:
					case 11:
					case 13:
					case 14:
						pos[count][1] = 3;
						break;
					case 15:
						pos[count][1] = 4;
						break;
					default:
						pos[count][1] = 1;
						break;
				}
				//pos[count][1] = __builtin_popcount(x_i.seq[j]);
				
				//Rprintf("%d (%d)",pos[count][0],pos[count][1]);
				count++;
			}
		}
		
		#pragma omp parallel for \
			private(x_l, j, k, l, p, n, o) \
			default(shared) \
			schedule(guided) \
			num_threads(nthreads)
		for (j = 0; j <= (count - max_probeLength); j++) {
			// score each possible oligo
			int combos = 1; // number of ambiguities
			int alt; // whether to alternate
			double score, specScore, lowScore, thisScore;
			double dG, dGmin;
			int MM, num, thisStart, thisEnd, lastCycle, thisCycle, cycles;
			int pers[maxCombos][max_probeLength]; // all permutations of the oligo
			for (k = 0; k < max_probeLength; k++)
				for (p = 0; p < maxCombos; p++)
					pers[p][k] = 0;
			int lengths[maxCombos]; // the length of each oligo
			double FAms[maxCombos]; // the FAm of each oligo
			int c; // count of oligos that meet criteria
			double min_hyb_eff, lowest_hyb_eff, hyb_eff; // minimum hybridization efficiency
			double hyb_effs[nMMs]; // hybridization efficiency of worst mismatches
			int MM_seqs[nMMs]; // sequence number of worst mismatches
			
			if (j > 0) // not first possible position
				if (x_i.seq[pos[j - 1][0]] & 0x1) // previous position is A
					continue; // already checked case with poly-T tail
			
			for (p = 0; p < maxCombos; p++)
				lengths[p] = 0;
			
			for (k = 0; k < max_probeLength; k++) {
				if (combos*pos[j + k][1] > maxCombos) {
					combos = maxCombos + 1;
					break;
				}
				if ((pos[j + k][1] > 1) && // ambiguity AND
					(!((combos & 1) ^ // XOR
					  (pos[j + k][1] & 1)) || // both odd or even OR
					((combos % pos[j + k][1]) == 0))) { // remainder is zero
					// alternate order every combo
					alt = 1;
				} else {
					alt = 0;
				}
				switch (x_i.seq[pos[j + k][0]]) {
					case 1: // A
						for (p = 0; p < maxCombos; p++) {
							pers[p][k] = 0;
						}
						break;
					case 2: // C
						for (p = 0; p < maxCombos; p++) {
							pers[p][k] = 1;
						}
						break;
					case 3: // M = A or C
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = (p + (int)(p/combos)) & 1;
							} else {
								pers[p][k] = p & 1;
							}
						}
						break;
					case 4: // G
						for (p = 0; p < maxCombos; p++) {
							pers[p][k] = 2;
						}
						break;
					case 5: // R = A or G
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = ((p + (int)(p/combos)) & 1)*2;
							} else {
								pers[p][k] = (p & 1)*2;
							}
						}
						break;
					case 6: // S = C or G
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = ((p + (int)(p/combos)) & 1) + 1;
							} else {
								pers[p][k] = (p & 1) + 1;
							}
						}
						break;
					case 7: // V = A or C or G
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = (p + (int)(p/combos)) % 3;
							} else {
								pers[p][k] = p % 3;
							}
						}
						break;
					case 8: // T
						for (p = 0; p < maxCombos; p++) {
							pers[p][k] = 3;
						}
						break;
					case 9: // W = A or T
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = ((p + (int)(p/combos)) & 1)*3;
							} else {
								pers[p][k] = (p & 1)*3;
							}
						}
						break;
					case 10: // Y = C or T
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = ((p + (int)(p/combos)) & 1)*2 + 1;
							} else {
								pers[p][k] = (p & 1)*2 + 1;
							}
						}
						break;
					case 11: // H = A or C or T
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = (p + (int)(p/combos)) % 3;
							} else {
								pers[p][k] = p % 3;
							}
							if (pers[p][k]==2)
								pers[p][k]++;
						}
						break;
					case 12: // K = G or T
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = ((p + (int)(p/combos)) & 1) + 2;
							} else {
								pers[p][k] = (p & 1) + 2;
							}
						}
						break;
					case 13: // D = A or G or T
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = (p + (int)(p/combos)) % 3;
							} else {
								pers[p][k] = p % 3;
							}
							if (pers[p][k] > 0)
								pers[p][k]++;
						}
						break;
					case 14: // B = C or G or T
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = ((p + (int)(p/combos)) % 3) + 1;
							} else {
								pers[p][k] = (p % 3) + 1;
							}
						}
						break;
					case 15: // N = A or C or G or T
						for (p = 0; p < maxCombos; p++) {
							if (alt) {
								pers[p][k] = (p + (int)(p/combos)) & 3;
							} else {
								pers[p][k] = p & 3;
							}
						}
						break;
				}
				combos *= pos[j + k][1];
				
				c = 0;
				for (p = 0; p < maxCombos; p++) {
					if (k==0) {
						FAms[p] = -27.9943679147442;
						lengths[p] = 0;
					} else {
						if (lengths[p]==0) {
							FAms[p] -= NN[pers[p][k - 1]][pers[p][k]];
							if (FAms[p] > min_fam) {
								lengths[p] = k;
								c++;
							}
						} else {
							c++;
						}
					}
				}
				
				if (c==maxCombos) // all FAms > min_fam
					break;
			}
			
			// enforce criteria for probe
			if (combos > maxCombos)
				continue;
			
			c = 0;
			for (p = 0; p < combos; p++) {
				if ((lengths[p] < (min_probeLength - 1)) || (lengths[p] >= max_probeLength))
					break;
				if ((FAms[p] < min_fam) || (FAms[p] > max_fam))
					break;
				c++;
				//Rprintf("\n%d\t%d\t%d",i,lengths[p],(int)FAms[p]);
			}
			
			if (c < combos) // not all oligos meet criteria
				continue;
			
			// calculate score
			specScore = 100; // overall score for target site
			for (o = 0; o < nMMs; o++) {
				hyb_effs[o] = 0;
				MM_seqs[o] = -1;
			}
			
			lowest_hyb_eff = .0005; // > 0% when rounded
			for (n = 0; n < combos; n++) { // for each permutation
				score = 100; // score for this permutation
				for (l = 0; l < x_length; l++) { // for each sequence
					if (i==l)
						continue;
					x_l = get_elt_from_XStringSet_holder(&x_set, l);
					MM = 0; // number of mismatches
					dG = dGini; // dG of probe
					dGmin = dGini; // lowest (best) dG along probe length
					c = lengths[n]; // tracks position in probe
					int I[2] = {0};
					int II[2] = {0};
					int III[2] = {0};
					int IV[2] = {0};
					for (p = pos[j + lengths[n]][0]; p >= start; p--) {
						if (!(x_i.seq[p] & 0x10 || x_i.seq[p] & 0x40) || !(x_l.seq[p] & 0x10 || x_l.seq[p] & 0x40)) {
							if (c < 0) { // tail of probe
								if (!(x_l.seq[p] & 0x10 || x_l.seq[p] & 0x40)) {
									c--;
									if (c < -6) // more than 5 consecutive matches
										break;
									I[0] = 0; // poly-T base (matches A)
									if (x_l.seq[p] & 0x1) { // A
										I[1] = -1;
									} else if (x_l.seq[p] & 0x4) { // G
										I[1] = 2;
									} else if (x_l.seq[p] & 0x8) { // T
										I[1] = 3;
									} else { // C
										I[1] = 1;
									} // no gap possible
								} else { // target site is gap
									continue;
								}
							} else {
								if (x_i.seq[p] & x_l.seq[p]) {
									I[0] = pers[n][c];
									c--;
									
									I[1] = -1; // perfect match
								} else if (!(x_i.seq[p] & 0x10 || x_i.seq[p] & 0x40)) { // mismatch
									I[0] = pers[n][c];
									c--;
									
									if (x_l.seq[p] & 0x4) { // G
										I[1] = 2;
									} else if (x_l.seq[p] & 0x1) { // A
										I[1] = 0;
									} else if (x_l.seq[p] & 0x8) { // T
										I[1] = 3;
									} else if (x_l.seq[p] & 0x2) { // C
										I[1] = 1;
									} else { // gap
										I[1] = 4;
									}
								}  else { // gap mismatch
									I[0] = 4;
									
									if (x_l.seq[p] & 0x4) { // G
										I[1] = 2;
									} else if (x_l.seq[p] & 0x1) { // A
										I[1] = 0;
									} else if (x_l.seq[p] & 0x8) { // T
										I[1] = 3;
									} else { // C
										I[1] = 1;
									} // no gap possible
								}
							}
							
							if (I[1] == -1) { // perfect match first (this) position
								if (II[1] == -1) { // perfect match second (next) position
									if ((dG != dGini) && (III[1] != -1) && (IV[1] == -1)) { // previous gap mismatch
										dG += sMM[II[0]][III[0]][III[1]][IV[0]];
										if (dG > dGini) // will not anneal
											dG = dGini;
									}
									dG += PM[I[0]][II[0]];
								} else if (dG != dGini) {
									if (III[1] == -1) { // perfect match in third position
										if (IV[1] == -1) { // perfect match in fourth position
											// add dG on next nearest neighbor pair
											MM++;
										} else { // mismatch in fourth position
											break; // mismatch every other position (becomes 3 MM)
										}
									} else { // mismatch in third position
										if ((II[0] == 4) || (III[0] == 4) || // mismatch pair contains a gap OR
											(IV[1] != -1)) { // three mismatches
											break; // too destabalizing
										} else { // tandem mismatch without any gaps
											// add dG of first MM pair
											if ((I[0] == 0) || (I[0] == 3)) { // AT or TA
												if (((II[0] == 1) && (II[1] == 3)) ||
													((II[0] == 0) && (II[1] == 2))) { // (G = C)T or (T = A)G
													dG += 0.652612195035552;
												} else if (((II[0] == 1) && (II[1] == 0)) ||
														   ((II[0] == 3) && (II[1] == 2))) { // (G = C)A or (A = T)G
													dG += 0.740231877389849;
												} else if ((II[0] == 1) && (II[1] == 2)) { // (G = C)G
													dG += 0.598920459808572;
												} else {
													dG += 1.09676430308156;
												}
											} else { // GC or CG
												if (((II[0] == 1) && (II[1] == 3)) ||
													((II[0] == 0) && (II[1] == 2))) { // (G = C)T or (T = A)G
													dG += 0.552804210078507;
												} else if (((II[0] == 1) && (II[1] == 0)) ||
														   ((II[0] == 3) && (II[1] == 2))) { // (G = C)A or (A = T)G
													dG += 0.375710199571309;
												} else if ((II[0] == 1) && (II[1] == 2)) { // (G = C)G
													dG += 0.797930010128296;
												} else {
													dG += 0.835053560606288;
												}
											}
											// add dG of second MM pair
											if ((IV[0] == 0) || (IV[0] == 3)) { // AT or TA
												if (((III[0] == 1) && (III[1] == 3)) ||
													((III[0] == 0) && (III[1] == 2))) { // (G = C)T or (T = A)G
													dG += 0.652612195035552;
												} else if (((III[0] == 1) && (III[1] == 0)) ||
														   ((III[0] == 3) && (III[1] == 2))) { // (G = C)A or (A = T)G
													dG += 0.740231877389849;
												} else if ((III[0] == 1) && (III[1] == 2)) { // (G = C)G
													dG += 0.598920459808572;
												} else {
													dG += 1.09676430308156;
												}
											} else { // GC or CG
												if (((III[0] == 1) && (III[1] == 3)) ||
													((III[0] == 0) && (III[1] == 2))) { // (G = C)T or (T = A)G
													dG += 0.552804210078507;
												} else if (((III[0] == 1) && (III[1] == 0)) ||
														   ((III[0] == 3) && (III[1] == 2))) { // (G = C)A or (A = T)G
													dG += 0.375710199571309;
												} else if ((III[0] == 1) && (III[1] == 2)) { // (G = C)G
													dG += 0.797930010128296;
												} else {
													dG += 0.835053560606288;
												}
											}
										}
										MM += 2;
									}
									
									if (dG > dGini) // will not anneal
										dG = dGini;
									
									if (MM > 2)
										break;
								}
							}
							
							if (dGmin > dG)
								dGmin = dG;
							
							//if (i==0)
							//	if (l==1)
							//		if (pos[j][0]==2059)
							//			Rprintf("\n%d\t%d\t%d",(int)(100*dG),I[0],I[1]);
							
							// rotate positions
							IV[0] = III[0];
							IV[1] = III[1];
							III[0] = II[0];
							III[1] = II[1];
							II[0] = I[0];
							II[1] = I[1];
						}
					}
					
					hyb_eff = 0.00995176243169644*exp(dGmin/-.626234565)/(1 + 0.00995176243169644*exp(dGmin/-.626234565));
					score -= 100*hyb_eff;
					
					if (dGmin > dG)
						dGmin = dG;
					
					if (hyb_eff > lowest_hyb_eff) {
						// record sequence number if top hybridization efficiency
						min_hyb_eff = 1.01; // initialize to 101%
						num = -1;
						if (n > 0) {
							for (o = 0; o < nMMs; o++) {
								if (MM_seqs[o] == l) { // another permutation
									if (hyb_eff > hyb_effs[o]) {
										num = o;
										break; // record this permutation
									} else { // already recorded
										num = -1;
										break; // do not rerecord
									}
								}
								if (hyb_effs[o] < min_hyb_eff) {
									if (hyb_eff > hyb_effs[o])
										num = o;
									lowest_hyb_eff = min_hyb_eff;
									min_hyb_eff = hyb_effs[o];
								}
								if (MM_seqs[o] == -1)
									break;
							}
						} else {
							for (o = 0; o < nMMs; o++) {
								if (hyb_effs[o] < min_hyb_eff) {
									if (hyb_eff > hyb_effs[o])
										num = o;
									lowest_hyb_eff = min_hyb_eff;
									min_hyb_eff = hyb_effs[o];
								}
								if (MM_seqs[o] == -1)
									break;
							}
						}
						if (num != -1) {
							hyb_effs[num] = hyb_eff;//-1*(dGmin - hyb_percent*0.17314);
							MM_seqs[num] = l;
							if (hyb_eff < lowest_hyb_eff)
								lowest_hyb_eff = hyb_eff;
						}
					}
				}
				
				if (score < specScore) // worst scoring permutation
					specScore = score;
				if (specScore < minScore)
					break;
			}
			
			if (specScore < minScore)
				continue;
			
			#pragma omp critical
			{
				// determine if probe is top scorer
				lowScore = 101;
				num = -1; // 0 to (numProbes - 1)
				for (p = 0; p < numProbes; p++) {
					thisScore = rans[i*numProbes + p + 4*x_length*numProbes];
					thisStart = rans[i*numProbes + p + 1*x_length*numProbes];
					thisEnd = thisStart + rans[i*numProbes + p + 2*x_length*numProbes] - 1;
					if (((j + 1) > (thisEnd  - maxOverlap)) || // start > end
						((j + k + 1) < (thisStart + maxOverlap))) { // end < start
						// less than maxOverlap between probes
						if (specScore >= thisScore) {
							if (thisScore <= lowScore) {
								// only replace lowest scored probe
								lowScore = thisScore;
								num = p;
							}
						}
					} else { // probes are overlapping
						if (specScore == thisScore) { // same score
							// probes have same score
							if (k + 1 < rans[i*numProbes + p + 2*x_length*numProbes]) {
								// pick shorter probe
								num = p;
								lowScore = 101;
								break;
							} else { // equal or longer probe
								num = -1; // do not replace
								lowScore = 101;
								break;
							}
						} else if (specScore > thisScore) { // higher score
							num = p;
							lowScore = 101;
							break;
						} else { // lower score
							num = -1; // do not replace
							lowScore = 101;
							break;
						}
					}
				}
				
				if (lowScore != 101)
					minScore = lowScore;
				
				if (num != -1) {
					// record probe if top scorer
					rans[i*numProbes + num] = i;
					rans[i*numProbes + num + x_length*numProbes] = j + 1;
					rans[i*numProbes + num + 2*x_length*numProbes] = k + 1;
					rans[i*numProbes + num + 3*x_length*numProbes] = combos;
					rans[i*numProbes + num + 4*x_length*numProbes] = specScore;
					rans[i*numProbes + num + 5*x_length*numProbes] = pos[j][0] + 1;
					rans[i*numProbes + num + 6*x_length*numProbes] = pos[j + k][0] + 1;
					double lowest_fam = max_fam;
					for (o = 0; o < combos; o++)
						if (FAms[o] < lowest_fam)
							lowest_fam = FAms[o];
					rans[i*numProbes + num + 7*x_length*numProbes] = lowest_fam;
					rans[i*numProbes + num + 8*x_length*numProbes] = 0.995176243169644*exp((lowest_fam*-0.17314 + hyb_percent*0.17314 + -2.8869448607588)/-0.626234565)/(1 + 0.00995176243169644*exp((lowest_fam*-0.17314 + hyb_percent*0.17314 + -2.8869448607588)/-0.626234565));
					
					char probe[k + 2]; // last position is for null terminating
					for (p = 0; p <= k; p++)
						probe[p] = DNAdecode(x_i.seq[pos[j + p][0]]);
					probe[k + 1] = '\0'; // end (null terminate) the string
					SET_STRING_ELT(probes, (i*numProbes + num), mkChar(probe));
					
					char perm[(k + 2 + polyT)*combos]; // last position is for null terminating
					c = 0;
					for (p = 0; p < combos; p++) {
						if (p != 0) {
							perm[c] = ',';
							c++;
						}
						
						// design cycles are counted 3' to 5'
						k = lengths[p];
						switch (pers[p][k]) {
							case 0:
								lastCycle = 3;
								perm[c] = 'T';
								break;
							case 1:
								lastCycle = 2;
								perm[c] = 'G';
								break;
							case 2:
								lastCycle = 1;
								perm[c] = 'C';
								break;
							case 3:
								lastCycle = 0;
								perm[c] = 'A';
								break;
						}
						c++;
						cycles = 4;
						
						for (k = lengths[p] - 1; k >= -1*polyT; k--) { // max poly-T length is 20 bases
							if (k < 0) { // poly-T base
								if (cycles > 148) // max of 148 design cycles
									break;
								thisCycle = 3;
								perm[c] = 'T';
							} else { // probe
								switch (pers[p][k]) {
									case 0:
										thisCycle = 3;
										perm[c] = 'T';
										break;
									case 1:
										thisCycle = 2;
										perm[c] = 'G';
										break;
									case 2:
										thisCycle = 1;
										perm[c] = 'C';
										break;
									case 3:
										thisCycle = 0;
										perm[c] = 'A';
										break;
								}
							}
							if (thisCycle < lastCycle) {
								cycles += lastCycle - thisCycle;
							} else {
								cycles += 4 + lastCycle - thisCycle;
							}
							lastCycle = thisCycle;
							c++;
						}
					}
					perm[c] = '\0'; // end (null terminate) the string
					SET_STRING_ELT(perms, (i*numProbes + num), mkChar(perm));
					
					for (o = 0; o < nMMs; o++) {
						rMMs[i*numProbes + num + o*x_length*numProbes] = hyb_effs[o]*100;
						rMMs[i*numProbes + num + (o + nMMs)*x_length*numProbes] = MM_seqs[o];
					}
				}
			}
		}
		
		if (v) {
			*rPercentComplete = floor(100*(double)(i + 1)/x_length);
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(ret_list, 0, ans);
	SET_VECTOR_ELT(ret_list, 1, probes);
	SET_VECTOR_ELT(ret_list, 2, perms);
	SET_VECTOR_ELT(ret_list, 3, MMs);
	
	if (v) {
		UNPROTECT(7);
	} else {
		UNPROTECT(5);
	}
	
	return ret_list;
}
