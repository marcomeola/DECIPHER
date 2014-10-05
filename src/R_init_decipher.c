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

/*
 * -- REGISTRATION OF THE .Call ENTRY POINTS ---
 */
static const R_CallMethodDef callMethods[] = { // method call, pointer, num args
	{"consensusSequence", (DL_FUNC) &consensusSequence, 6},
	{"consensusSequenceAA", (DL_FUNC) &consensusSequenceAA, 6},
	{"clusterNJ", (DL_FUNC) &clusterNJ, 5},
	{"reclusterNJ", (DL_FUNC) &reclusterNJ, 2},
	{"clusterUPGMA", (DL_FUNC) &clusterUPGMA, 6},
	{"reclusterUPGMA", (DL_FUNC) &reclusterUPGMA, 2},
	{"clusterML", (DL_FUNC) &clusterML, 6},
	{"distMatrix", (DL_FUNC) &distMatrix, 9},
	{"gaps", (DL_FUNC) &gaps, 2},
	{"designProbes", (DL_FUNC) &designProbes, 16},
	{"commonGaps", (DL_FUNC) &commonGaps, 1},
	{"multiMatch", (DL_FUNC) &multiMatch, 3},
	{"multiMatchUpper", (DL_FUNC) &multiMatchUpper, 3},
	{"multiMatchCharNotNA", (DL_FUNC) &multiMatchCharNotNA, 1},
	{"replaceChars", (DL_FUNC) &replaceChars, 3},
	{"replaceChar", (DL_FUNC) &replaceChar, 3},
	{"trimChar", (DL_FUNC) &trimChar, 2},
	{"intMatch", (DL_FUNC) &intMatch, 3},
	{"firstMatchUpper", (DL_FUNC) &firstMatchUpper, 3},
	{"terminalMismatch", (DL_FUNC) &terminalMismatch, 5},
	{"NNLS", (DL_FUNC) &NNLS, 10},
	{"sparseMult", (DL_FUNC) &sparseMult, 6},
	{"calculateDeltaG", (DL_FUNC) &calculateDeltaG, 3},
	{"calculateFISH", (DL_FUNC) &calculateFISH, 2},
	{"alignProfiles", (DL_FUNC) &alignProfiles, 11},
	{"alignProfilesAA", (DL_FUNC) &alignProfilesAA, 11},
	{"consensusProfile", (DL_FUNC) &consensusProfile, 2},
	{"consensusProfileAA", (DL_FUNC) &consensusProfileAA, 2},
	{"matchLists", (DL_FUNC) &matchLists, 4},
	{"adjustHeights", (DL_FUNC) &adjustHeights, 1},
	{"enumerateSequence", (DL_FUNC) &enumerateSequence, 2},
	{"enumerateSequenceAA", (DL_FUNC) &enumerateSequenceAA, 2},
	{"enumerateSequenceAA", (DL_FUNC) &enumerateSequenceReducedAA, 3},
	{"enumerateGappedSequence", (DL_FUNC) &enumerateGappedSequence, 2},
	{"enumerateGappedSequenceAA", (DL_FUNC) &enumerateGappedSequenceAA, 2},
	{"matchOrder", (DL_FUNC) &matchOrder, 4},
	{"matchRanges", (DL_FUNC) &matchRanges, 5},
	{"firstSeqsEqual", (DL_FUNC) &firstSeqsEqual, 6},
	{"matchOrderDual", (DL_FUNC) &matchOrderDual, 3},
	{"boundedMatches", (DL_FUNC) &boundedMatches, 3},
	{"gcContent", (DL_FUNC) &gcContent, 3},
	{"intDist", (DL_FUNC) &intDist, 6},
	{"meltPolymer", (DL_FUNC) &meltPolymer, 4},
	{NULL, NULL, 0}
};

void R_init_DECIPHER(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
