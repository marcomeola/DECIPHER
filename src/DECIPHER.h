// R_init_DECIPHER.c

void R_init_DECIPHER(DllInfo *info);

// ConsensusSequence.c

SEXP consensusSequence(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps);

SEXP consensusSequenceAA(SEXP x, SEXP threshold, SEXP ambiguity, SEXP minInformation, SEXP ignoreNonLetters, SEXP terminalGaps);

SEXP consensusProfile(SEXP x, SEXP weight, SEXP structs);

SEXP consensusProfileAA(SEXP x, SEXP weight, SEXP structs);

SEXP colScores(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP weights, SEXP structs, SEXP dbnMatrix);

SEXP colScoresAA(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP weights, SEXP structs, SEXP hecMatrix);

SEXP shiftGaps(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP gl, SEXP sc, SEXP thresh, SEXP weights);

SEXP shiftGapsAA(SEXP x, SEXP subMatrix, SEXP go, SEXP ge, SEXP gl, SEXP sc, SEXP thresh, SEXP weights);

// DistanceMatrix.c

SEXP distMatrix(SEXP x, SEXP t, SEXP terminalGaps, SEXP penalizeGapGaps, SEXP penalizeGapLetters, SEXP fullMatrix, SEXP verbose, SEXP pBar, SEXP nThreads);

SEXP gaps(SEXP x, SEXP t);

SEXP firstSeqsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y);

SEXP firstSeqsGapsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP t);

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

SEXP boundedMatches(SEXP x, SEXP bl, SEXP bu);

SEXP intMatchOnce(SEXP x, SEXP y);

SEXP matchListsDual(SEXP x, SEXP y, SEXP verbose, SEXP pBar, SEXP nThreads);

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

SEXP alignProfiles(SEXP p, SEXP s, SEXP subMatrix, SEXP dbnMatrix, SEXP pm, SEXP mm, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads);

SEXP alignProfilesAA(SEXP p, SEXP s, SEXP subMatrix, SEXP hecMatrix, SEXP go, SEXP ge, SEXP exp, SEXP power, SEXP endGapPenaltyLeft, SEXP endGapPenaltyRight, SEXP boundary, SEXP nThreads);

// EnumerateSequence.c

SEXP enumerateSequence(SEXP x, SEXP wordSize);

SEXP enumerateSequenceAA(SEXP x, SEXP wordSize);

SEXP enumerateGappedSequence(SEXP x, SEXP wordSize, SEXP ordering);

SEXP enumerateGappedSequenceAA(SEXP x, SEXP wordSize, SEXP ordering);

SEXP enumerateSequenceReducedAA(SEXP x, SEXP wordSize, SEXP alphabet);

// Compositions.c

SEXP gcContent(SEXP x, SEXP begins, SEXP ends);

SEXP composition(SEXP x);

SEXP positionWeightMatrix(SEXP x, SEXP begins, SEXP ends, SEXP width);

// IntDist.c

SEXP intDist(SEXP x, SEXP levels, SEXP bins, SEXP maxBins, SEXP numRows, SEXP totRows);

// MeltPolymer.c

SEXP meltPolymer(SEXP x, SEXP temps, SEXP ions, SEXP output);

// InsertGaps.c

SEXP insertGaps(SEXP x, SEXP positions, SEXP lengths, SEXP type, SEXP nThreads);

// ExpandAmbiguities.c

SEXP expandAmbiguities(SEXP x, SEXP c);

// RemoveGaps.c

SEXP removeCommonGaps(SEXP x, SEXP type, SEXP nThreads);

// SubsetXStringSet.c

SEXP subsetXStringSet(SEXP x, SEXP subset, SEXP type, SEXP nThreads);

// AppendXStringSets.c

SEXP appendXStringSets(SEXP x, SEXP y, SEXP type, SEXP nThreads);

// PredictHEC.c

SEXP predictHEC(SEXP x, SEXP windowSize, SEXP background, SEXP HEC_MI1, SEXP HEC_MI2, SEXP output);

// AssignIndels.c

SEXP clearIns(SEXP x);

SEXP all(SEXP x);

SEXP any(SEXP x);

// ConsolidateGaps.c

SEXP consolidateGaps(SEXP x, SEXP type);

// FindFrameshifts.c

SEXP findFrameshifts(SEXP t, SEXP l, SEXP f, SEXP index, SEXP maxComp, SEXP go, SEXP ge, SEXP fs, SEXP minD, SEXP maxD, SEXP subMatrix, SEXP verbose, SEXP pBar);

// Order.c

SEXP radixOrder(SEXP x);

// ChainSegments.c

SEXP fillOverlaps(SEXP m, SEXP n);

SEXP indexByContig(SEXP starts, SEXP ends, SEXP order, SEXP index, SEXP widths);

SEXP chainSegments(SEXP x_s, SEXP x_e, SEXP x_i, SEXP x_f, SEXP y_s, SEXP y_e, SEXP y_i, SEXP y_f, SEXP weights, SEXP sepCost, SEXP gapCost, SEXP shiftCost, SEXP codingCost, SEXP maxSep, SEXP maxGap, SEXP ordering, SEXP minScore);

SEXP extendSegments(SEXP X, SEXP W1, SEXP W2, SEXP S1, SEXP S2, SEXP O1P, SEXP O1N, SEXP O2P, SEXP O2N, SEXP S, SEXP maxDrop, SEXP INDEX1, SEXP INDEX2);

SEXP maskRepeats(SEXP e, SEXP size, SEXP minL, SEXP maxL, SEXP totL);

// Translate.c

SEXP basicTranslate(SEXP x, SEXP code, SEXP starts);
