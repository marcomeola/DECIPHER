/****************************************************************************
 *                        Quickly Replace Characters                        *
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

//ans_start <- .Call("replaceChars", sequences, replaceChar, PACKAGE="DECIPHER")
SEXP replaceChars(SEXP x, SEXP r)
{
	int i, j, l, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	const char *repChar = CHAR(STRING_ELT(r, 0));
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char s[longest + 1]; // each sequence
	
	// write new character vector
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		count = 0;
		for (j = 0; j < l; j++) {
			if (seq[j]!='U' && seq[j]!='u') {
				switch (seq[j]) {
					case '-':
					case 'A':
					case 'a':
					case 'C':
					case 'c':
					case 'G':
					case 'g':
					case 'T':
					case 't':
					case 'N':
					case 'n':
					case 'M':
					case 'm':
					case 'R':
					case 'r':
					case 'W':
					case 'w':
					case 'S':
					case 's':
					case 'Y':
					case 'y':
					case 'K':
					case 'k':
					case 'V':
					case 'v':
					case 'H':
					case 'h':
					case 'D':
					case 'd':
					case 'B':
					case 'b':
					case '+':
						s[count] = seq[j];
						count++;
						break;
					default:
						if (repChar[0] != '\0') {
							s[count] = repChar[0];
							count++;
						}
						break;
				}
			} else {
				s[count] = 'T';
				count++;
			}
		}
		s[count] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(s));
	}
	
	UNPROTECT(1);
	
	return seqs;
}

//ans_start <- .Call("replaceChar", sequences, charReplace, replaceChar, PACKAGE="DECIPHER")
SEXP replaceChar(SEXP x, SEXP c, SEXP r)
{
	int i, j, l, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	const char *repChar = CHAR(STRING_ELT(r, 0));
	const char *charRep = CHAR(STRING_ELT(c, 0));
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char s[longest + 1]; // each sequence
	
	// write new character vector
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		count = 0;
		for (j = 0; j < l; j++) {
			if (seq[j]==charRep[0]) {
				if (repChar[0] != '\0') {
					s[count] = repChar[0];
					count++;
				}
			} else {
				s[count] = seq[j];
				count++;
			}
		}
		s[count] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(s));
	}
	
	UNPROTECT(1);
	
	return seqs;
}

//ans_start <- .Call("trimChar", sequences, numChar, PACKAGE="DECIPHER")
SEXP trimChar(SEXP x, SEXP y)
{
	int i, j, l;
	int num = asInteger(y);
	int n = length(x);
	int longest = 0;
	const char *seq;
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char s[longest + 1 - num]; // each sequence
	
	// write new character vector
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		for (j = 0; j < (l - num); j++)
			s[j] = seq[j];
		s[j] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(s));
	}
	
	UNPROTECT(1);
	
	return seqs;
}
