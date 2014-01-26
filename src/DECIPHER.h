// R_init_DECIPHER.c

void R_init_DECIPHER(DllInfo *info);

// ConsensusSequence.c

SEXP consensusSequence(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps);

SEXP consensusSequenceAA(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps);

SEXP consensusProfile(SEXP x, SEXP weight);

SEXP consensusProfileAA(SEXP x, SEXP weight);

// DistanceMatrix.c

SEXP distMatrix(SEXP x, SEXP t, SEXP terminalGaps, SEXP penalizeGapGaps, SEXP penalizeGapLetters, SEXP fullMatrix, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP gaps(SEXP x, SEXP t);

SEXP firstSeqsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y);

// ClusterNJ.c

SEXP clusterNJ(SEXP x, SEXP cutoff, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP reclusterNJ(SEXP ans, SEXP cutoff);

SEXP adjustHeights(SEXP x);

// ClusterUPGMA.c

SEXP clusterUPGMA(SEXP x, SEXP cutoff, SEXP method, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP reclusterUPGMA(SEXP ans, SEXP cutoff);

// ClusterML.c

SEXP clusterML(SEXP x, SEXP y, SEXP model, SEXP branches, SEXP lengths, SEXP nThreads);

// DesignProbes.c

SEXP designProbes(SEXP x, SEXP max_pl, SEXP min_pl, SEXP max_c, SEXP numMMs, SEXP numPs, SEXP st, SEXP en, SEXP max_ov, SEXP h_percent, SEXP min_f, SEXP max_f, SEXP minS, SEXP verbose, SEXP pBar, SEXP nThreads);

// CommonGaps.c

SEXP commonGaps(SEXP x);

// MultiMatch.c

SEXP multiMatch(SEXP x, SEXP y, SEXP z);

SEXP multiMatchUpper(SEXP x, SEXP y, SEXP z);

SEXP multiMatchCharNotNA(SEXP x);

SEXP intMatch(SEXP x, SEXP y, SEXP nThreads);

SEXP firstMatchUpper(SEXP x, SEXP y, SEXP nThreads);

SEXP matchLists(SEXP x, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP matchOrder(SEXP x, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP matchOrderDual(SEXP x, SEXP y, SEXP nThreads);

SEXP matchRanges(SEXP x, SEXP y, SEXP wordSize, SEXP maxLength, SEXP threshold);

// ReplaceChars.c

SEXP replaceChars(SEXP x, SEXP r, SEXP t);

SEXP replaceChar(SEXP x, SEXP c, SEXP r);

SEXP trimChar(SEXP x, SEXP y);

// TerminalMismatch.c

SEXP terminalMismatch(SEXP p, SEXP t, SEXP cutoff, SEXP mGaps, SEXP nThreads);

// NNLS.c

SEXP NNLS(SEXP row, SEXP col, SEXP value, SEXP nrows, SEXP ncols, SEXP b, SEXP tol, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP sparseMult(SEXP row, SEXP col, SEXP value, SEXP nrows, SEXP ncols, SEXP b);

// CalculateDeltaG.c

SEXP calculateDeltaG(SEXP p, SEXP t, SEXP deltaGrules);

// CalculateFISH.c

SEXP calculateFISH(SEXP probes, SEXP targets);

// AlignProfiles.c

SEXP alignProfiles(SEXP p, SEXP s, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads);

SEXP alignProfilesAA(SEXP p, SEXP s, SEXP subMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads);

// EnumerateSequence.c

SEXP enumerateSequence(SEXP x, SEXP wordSize);

SEXP enumerateSequenceAA(SEXP x, SEXP wordSize);

SEXP enumerateGappedSequence(SEXP x, SEXP wordSize);

SEXP enumerateGappedSequenceAA(SEXP x, SEXP wordSize);
