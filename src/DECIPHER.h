// ConsensusSequence.c

SEXP consensusSequence(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps);

// DistanceMatrix.c

SEXP distMatrix(SEXP x, SEXP terminalGaps, SEXP penalizeGapGaps, SEXP penalizeGapLetters, SEXP verbose, SEXP pBar);

SEXP gaps(SEXP x);

// R_init_DECIPHER.c

void R_init_DECIPHER(DllInfo *info);

// ClusterNJ.c

SEXP clusterNJ(SEXP x, SEXP cutoff, SEXP verbose, SEXP pBar);

SEXP reclusterNJ(SEXP ans, SEXP cutoff);

// ClusterUPGMA.c

SEXP clusterUPGMA(SEXP x, SEXP cutoff, SEXP method, SEXP verbose, SEXP pBar);

SEXP reclusterUPGMA(SEXP ans, SEXP cutoff);

// ClusterML.c

SEXP clusterML(SEXP x, SEXP y, SEXP cutoff, SEXP verbose, SEXP pBar);

// CommonGaps.c

SEXP commonGaps(SEXP x);

// ReplaceChars.c

SEXP replaceChars(SEXP x, SEXP r);

SEXP replaceChar(SEXP x, SEXP c, SEXP r);

SEXP trimChar(SEXP x, SEXP y);

// CalculateDeltaG.c

SEXP calculateDeltaG(SEXP p, SEXP t, SEXP deltaGrules);
