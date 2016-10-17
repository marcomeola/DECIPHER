/****************************************************************************
 *                       Helper Functions for Import                        *
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

// DECIPHER header file
#include "DECIPHER.h"

// concatenate x between index1 and index2 indices
SEXP collapse(SEXP x, SEXP index1, SEXP index2)
{
	int i, j, k, l, tot, count;
	char *s;
	const char *x_i;
	int n = length(index1);
	int *i1 = INTEGER(index1);
	int *i2 = INTEGER(index2);
	
	SEXP ans;
	PROTECT(ans = allocVector(STRSXP, n));
	
	// write new character vector
	for (i = 0; i < n; i++) {
		tot = 1; // 1 for null-termination
		for (j = i1[i] - 1; j <  i2[i]; j++) {
			tot += length(STRING_ELT(x, j));
		}
		
		s = Calloc(tot, char);
		count = 0;
		for (j = i1[i] - 1; j <  i2[i]; j++) {
			l = length(STRING_ELT(x, j));
			x_i = CHAR(STRING_ELT(x, j));
			for (k = 0; k < l; k++) {
				s[count] = x_i[k];
				count++;
			}
		}
		s[count] = '\0'; // null-terminate
		SET_STRING_ELT(ans, i, mkChar(s));
		Free(s);
	}
	
	UNPROTECT(1);
	
	return ans;
}

// extract GenBank fields between sequence starts
SEXP extractFields(SEXP x, SEXP fields, SEXP starts, SEXP ends)
{
	int i, j, k, l1, l2, o, p, hit, curlen;
	int n = length(fields);
	int l = length(starts);
	const char *str, *field;
	int maxlen = 1000;
	char *s = Calloc(maxlen, char);
	int *ss = INTEGER(starts);
	int *se = INTEGER(ends);
	
	SEXP ret, ans;
	PROTECT(ret = allocVector(VECSXP, n));
	
	for (i = 0; i < n; i++) {
		field = CHAR(STRING_ELT(fields, i));
		l1 = length(STRING_ELT(fields, i));
		PROTECT(ans = allocVector(STRSXP, l));
		
		for (j = 0; j < l; j++) {
			curlen = 0;
			
			for (k = ss[j] - 1; k < se[j]; k++) {
				str = CHAR(STRING_ELT(x, k));
				l2 = length(STRING_ELT(x, k));
				
				hit = 0;
				o = 0;
				p = 0;
				
				// find leading spaces
				while (p < l2 && p < 12) {
					if (str[p] != ' ')
						break;
					p++;
				}
				
				// search for hit
				while (o < l1 && p < l2 && p < 12) {
					if (str[p] != field[o]) {
						break;
					} else {
						p++;
						if (o==(l1 - 1)) {
							hit = 1;
							break;
						}
						o++;
					}
				}
				
				// find trailing spaces
				if (hit==1) {
					while (p < l2 && p < 12) {
						if (str[p] != ' ') {
							hit = 0;
							break;
						}
						p++;
					}
					
					if (hit==1) {
						if (curlen != 0) {
							if (curlen==maxlen) {
								maxlen *= 2;
								s = Realloc(s, maxlen, char);
							}
							s[curlen] = '\n';
							curlen++;
							
							if (curlen==maxlen) {
								maxlen *= 2;
								s = Realloc(s, maxlen, char);
							}
							s[curlen] = '\n';
							curlen++;
						}
						
						while (p < l2) {
							if (curlen==maxlen) {
								maxlen *= 2;
								s = Realloc(s, maxlen, char);
							}
							s[curlen] = str[p];
							curlen++;
							p++;
						}
						
						k++;
						while (k < se[j]) {
							str = CHAR(STRING_ELT(x, k));
							l2 = length(STRING_ELT(x, k));
							
							p = 0;
							while (p < l2 && p < 12) {
								if (str[p] != ' ') {
									hit = 0;
									break;
								}
								p++;
							}
							
							if (hit==1) {
								if (curlen==maxlen) {
									maxlen *= 2;
									s = Realloc(s, maxlen, char);
								}
								s[curlen] = '\n';
								curlen++;
								
								while (p < l2) {
									if (curlen==maxlen) {
										maxlen *= 2;
										s = Realloc(s, maxlen, char);
									}
									s[curlen] = str[p];
									curlen++;
									p++;
								}
							} else {
								k--;
								break;
							}
							k++;
						}
					}
				}
			}
			
			if (curlen==maxlen) {
				maxlen *= 2;
				s = Realloc(s, maxlen, char);
			}
			s[curlen] = '\0';
			
			SET_STRING_ELT(ans, j, mkChar(s));
		}
		
		SET_VECTOR_ELT(ret, i, ans);
		UNPROTECT(1);
	}
	
	Free(s);
	UNPROTECT(1); // ret
	
	return ret;
}
