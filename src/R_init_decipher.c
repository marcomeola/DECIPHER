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
	{"clusterNJ", (DL_FUNC) &clusterNJ, 4},
	{"reclusterNJ", (DL_FUNC) &reclusterNJ, 2},
	{"clusterUPGMA", (DL_FUNC) &clusterUPGMA, 5},
	{"reclusterUPGMA", (DL_FUNC) &reclusterUPGMA, 2},
	{"clusterML", (DL_FUNC) &clusterML, 5},
	{"distMatrix", (DL_FUNC) &distMatrix, 6},
	{"gaps", (DL_FUNC) &gaps, 1},
	{"commonGaps", (DL_FUNC) &commonGaps, 1},
	{"replaceChars", (DL_FUNC) &replaceChars, 2},
	{"replaceChar", (DL_FUNC) &replaceChar, 3},
	{"trimChar", (DL_FUNC) &trimChar, 2},
	{NULL, NULL, 0}
};

void R_init_DECIPHER(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
