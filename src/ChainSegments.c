/****************************************************************************
 *                         Chains Graph of Segments                         *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

// in-place fills missing sequence in the input
SEXP fillOverlaps(SEXP m, SEXP n)
{	
	int i;
	int *x = INTEGER(m);
	int l = length(m);
	int y = asInteger(n);
	
	int k = y - 1; // current position
	int last = k; // last non-incrementing position
	int p = 0; // last position
	
	while (k < l) {
		if (x[k]!=NA_INTEGER &&
			x[p]!=NA_INTEGER &&
			x[k]==(x[p] + y - 1)) {
			if (last > p) {
				for (i = p + 1; i < k; i++)
					x[i] = x[i - 1] + 1;
				last = p;
			}
		} else {
			last = k;
		}
		k++;
		p++;
	}
	
	return R_NilValue;
}

// in-place adjusts starts and ends in order to be within widths
SEXP indexByContig(SEXP starts, SEXP ends, SEXP order, SEXP index, SEXP widths)
{	
	int j, k, p;
	int *s = INTEGER(starts);
	int *e = INTEGER(ends);
	int *o = INTEGER(order);
	int *w = INTEGER(widths);
	int *i = INTEGER(index);
	int l = length(starts);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l));
	int *rans = INTEGER(ans);
	
	// fill initial values
	for (j = 0; j < l; j++) {
		p = o[j] - 1;
		if (s[p] > w[0])
			break;
		rans[p] = i[0];
	}
	
	// index sequences by contig
	k = 1;
	for (; j < l; j++) { // j = j
		p = o[j] - 1;
		while (s[p] > w[k]) {
			k++;
		}
		s[p] = s[p] - w[k - 1];
		e[p] = e[p] - w[k - 1];
		rans[p] = i[k];
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP chainSegments(SEXP x_s, SEXP x_e, SEXP x_i, SEXP x_f, SEXP y_s, SEXP y_e, SEXP y_i, SEXP y_f, SEXP weights, SEXP sepCost, SEXP gapCost, SEXP shiftCost, SEXP codingCost, SEXP maxSep, SEXP maxGap, SEXP ordering, SEXP minScore)
{	
	int *xs = INTEGER(x_s);
	int *xe = INTEGER(x_e);
	int *xi = INTEGER(x_i);
	int *xf = INTEGER(x_f);
	int *ys = INTEGER(y_s);
	int *ye = INTEGER(y_e);
	int *yi = INTEGER(y_i);
	int *yf = INTEGER(y_f);
	double *we = REAL(weights);
	double sepC = asReal(sepCost);
	double gapC = asReal(gapCost);
	double shiC = asReal(shiftCost);
	double codC = asReal(codingCost);
	double maxS = asReal(maxSep);
	double maxG = asReal(maxGap);
	int *xo = INTEGER(ordering);
	double minS = asReal(minScore);
	
	int l = length(x_s);
	int i = 0;
	int j = -1;
	int prev = 0;
	int k, max, dy, dx, sep, gap, xok, xoj;
	double temp, score;
	
	// initialize an array of activities
	int *A = Calloc(l, int); // initialized to zero
	// initialize an array of scores
	double *S = Calloc(l, double); // initialized to zero
	// initialize an array of return indicies
	int *R = Calloc(l, int); // initialized to zero
	// initialize an array of origin indices
	int *O = Calloc(l, int); // initialized to zero
	
	while (i < l || j < (l - 1)) {
		if (i < l &&
			xs[i] <= xe[xo[j + 1]] &&
			xi[i]==xi[xo[j + 1]]) {
			// use i, left of rectangle
			if (i > 0 &&
				xi[i] != xi[i - 1])
				prev = 0; // deactivate left
			
			while (prev > 0 &&
				   (xs[i] - xe[xo[j - prev + 1]]) > maxS) {
				prev--;
			}
			
			// find the highest scoring rectangle below the start
			max = -1;
			score = 0;
			for (k = j - prev + 1; k <= j; k++) {
				xok = xo[k];
				
				if (A[xok] != 1 ||
					yi[i] != yi[xok])
					continue;
				
				dy = ys[i] - ye[xok] - 1;
				if (dy <= 0)
					continue;
				
				dx = xs[i] - xe[xok] - 1;
				if (dx > dy) {
					sep = dy;
					gap = dx - dy;
				} else {
					sep = dx;
					gap = dy - dx;
				}
				if (sep > maxS ||
					gap > maxG)
					continue;
				
				// add cost for shifting reading frames
				if (xf[i]==xf[xok] &&
					(xf[i]==0 || (dx % 3)==0)) {
					temp = 0;
				} else {
					if (xf[i]==0 || xf[xok]==0) {
						temp = codC;
					} else {
						temp = shiC;
					}
				}
				if (!(yf[i]==yf[xok] &&
					(yf[i]==0 || (dy % 3)==0))) {
					if (yf[i]==0 || yf[xok]==0) {
						temp = codC;
					} else {
						temp = shiC;
					}
				}
				
				temp += sep*sepC;
				temp += gap*gapC;
				temp += S[xok];
				//temp += we[i];
				if (sep > 0 && gap > 0) {
					temp -= log((double)(sep*gap));
				} else if (sep > 0) {
					temp -= log((double)sep);
				} else if (gap > 0) {
					temp -= log((double)gap);
				}
				
				if (temp > score) {
					max = xok;
					score = temp;
				}
			}
			//Rprintf("\ni = %d score = %1.2f xo[i] = %d we[i] = %1.2f", i, score, xo[i], we[i]);
			if (score > 0) {
				S[i] = score + we[i];
				R[i] = max;// + 1;
				O[i] = O[max];// + 1;
			} else {
				S[i] = we[i];
				O[i] = i;// + 1;
				R[i] = -1;//0;
			}
			
			i++;
		} else {
			// use j, right of rectangle
			j++;
			xoj = xo[j];
			
			// find the nearest rectangle below the end
			max = -1;
			score = 1e53;
			for (k = j - prev; k < j; k++) {
				xok = xo[k];
				if (A[xok] != 1 ||
					yi[xoj] != yi[xok])
					continue;
				
				dy = ye[xoj] - ye[xok];
				if (dy >= 0 &&
					dy < score) {
					max = xok;
					score = dy;
				}
			}
			
			if (max != -1 &&
				S[xoj] > S[max]) {
				// find rectangles above the end with lower score and the same origin
				for (k = j - prev; k < j; k++) {
					xok = xo[k];
					if (ye[xok] >= ye[xoj] &&
						S[xok] < S[xoj] &&
						yi[xok]==yi[xoj] &&
						O[xok]==O[xoj])
						A[xok] = 0; // deactivate
				}
			}
			
			A[xoj] = 1;
			prev++;
		}
	}
	
	int count = 0;
	int *p;
	int n, min, overlap;
	int size = 1000;
	int **ptrs = Calloc(size, int *); // chains
	int *lens = Calloc(size, int); // length of each chain
	double *scores = Calloc(size, double); // score of each chain
	// start and end of the rectangle encompassing each chain
	int *rectXS = Calloc(size, int);
	int *rectXE = Calloc(size, int);
	int *rectYS = Calloc(size, int);
	int *rectYE = Calloc(size, int);
	int *rectXI = Calloc(size, int);
	int *rectYI = Calloc(size, int);
	
	for (i = 0; i < l; i++)
		A[i] = 1; // re-activate all
	
	while (1) {
		// find the highest active score
		max = 0;
		score = 0;
		for (i = 0; i < l; i++) {
			if (A[i] == 1 &&
				S[i] > score) {
				score = S[i];
				max = i;
			}
		}
		if (score < minS)
			break;
		
		j = 0; // length of chain
		// find the maximum length chain
		// with the same origin as max
		// and passing through max score
		for (i = l - 1; i >= O[max]; i--) {
			if (O[i]==O[max] &&
				A[i]==1) {
				n = 1;
				k = i;
				if (S[k]==score) {
					overlap = 1;
				} else {
					overlap = 0;
				}
				A[i] = 0; // deactivate
				while (R[k] != -1) {
					k = R[k];
					if (A[k] != 0)
						A[k] = 0; // deactivate
					if (S[k]==score)
						overlap = 1;
					n++;
				}
				if (overlap==1 &&
					n > j) {
					max = i;
					j = n;
				}
			}
		}
		
		// score must be increasing
		while (R[max] != -1 &&
			S[max] < S[R[max]]) {
			max = R[max];
		}
		
		// pull back overlapping end
		while (max != -1 &&
			max >= O[max]) {
			overlap = 0;
			for (i = 0; i < count; i++) {
				if (xi[max]==rectXI[i] &&
					yi[max]==rectYI[i] &&
					((xe[max] >= rectXS[i] &&
					xe[max] <= rectXE[i]) ||
					(ye[max] >= rectYS[i] &&
					ye[max] <= rectYE[i]))) {
					overlap = 1;
					break;
				}
			}
			if (overlap==0)
				break;
			max = R[max];
		}
		
		if (max==-1) // no chain
			continue;
		
		// pull back overlapping start
		min = max;
		n = min;
		j = 0;
		while (n >= O[max]) {
			overlap = 0;
			for (i = 0; i < count; i++) {
				if (xi[n]==rectXI[i] &&
					yi[n]==rectYI[i] &&
					((xs[n] >= rectXS[i] &&
					xs[n] <= rectXE[i]) ||
					(ys[n] >= rectYS[i] &&
					ys[n] <= rectYE[i]))) {
					overlap = 1;
					break;
				}
			}
			if (overlap==1)
				break;
			j++;
			min = n; // prior value of min
			n = R[min];
		}
		
		if (j==0)
			continue; // fully overlapping
		
		// shorten to reach minScore
		score = S[max] - S[min] + we[min];
		while (score < minS && max > min) {
			max = R[max]; // shorten from end
			j--;
			score = S[max] - S[min] + we[min];
		}
		
		if (score >= minS) {
			// find the nearest rectangle
			int minDx = 2e9, minDy = 2e9, minX = -1, minY = -2, merge = 0, upX, upY;
			for (i = 0; i < count; i++) {
				if (xi[max]==rectXI[i] &&
					yi[max]==rectYI[i]) {
					dx = xs[min] - rectXE[i];
					if (dx < 0) {
						dx = rectXS[i] - xe[max];
					}
					if (dx < minDx) {
						minDx = dx;
						minX = i;
						upX = (xs[min] > rectXE[i]) ? 1 : 0;
					}
					dy = ys[min] - rectYE[i];
					if (dy < 0) {
						dy = rectYS[i] - ye[max];
					}
					if (dy < minDy) {
						minDy = dy;
						minY = i;
						upY = (ys[min] > rectYE[i]) ? 1 : 0;
					}
				}
			}
			if (minDx < 0 || minDy < 0) // fully overlapping
				continue;
			
			sep = 1e9;
			gap = 1e9;
			if (minX == minY && upX == upY) {
				if (minDx < minDy) {
					sep = minDx;
					gap = minDy - minDx;
				} else {
					sep = minDy;
					gap = minDx - minDy;
				}
				if (gap <= maxG && sep <= maxS) {
					temp = scores[minX] + score + sepC*sep + gapC*gap;
					if (temp >= minS)
						merge = 1;
				}
			}
			
			if (merge) {
				ptrs[minX] = Realloc(ptrs[minX], j + lens[minX], int);
				p = ptrs[minX];
				
				if (upX) { // new chain is last
					j = j + lens[minX];
					lens[minX] = j;
					
					// add new chain at the end
					i = max;
					p[--j] = i;
					while (R[i] >= min) {
						i = R[i];
						p[--j] = i;
					}
					
					rectXE[minX] = xe[max];
					rectYE[minX] = ye[max];
				} else { // new chain is first
					// shift old chain to the end
					for (i = lens[minX] - 1; i >= 0; i--)
						p[i + j] = p[i];
					
					lens[minX] = j + lens[minX];
					
					// add new chain at the beginning
					i = max;
					p[--j] = i;
					while (R[i] >= min) {
						i = R[i];
						p[--j] = i;
					}
					
					rectXS[minX] = xs[min];
					rectYS[minX] = ys[min];
				}
				
				scores[count] = temp;
				
				continue;
			}
			
			if (count >= size) {
				size += 1000;
				ptrs = Realloc(ptrs, size, int *);
				lens = Realloc(lens, size, int);
				scores = Realloc(scores, size, double);
				rectXS = Realloc(rectXS, size, int);
				rectXE = Realloc(rectXE, size, int);
				rectYS = Realloc(rectYS, size, int);
				rectYE = Realloc(rectYE, size, int);
				rectXI = Realloc(rectXI, size, int);
				rectYI = Realloc(rectYI, size, int);
			}
			
			ptrs[count] = Calloc(j, int);
			lens[count] = j;
			p = ptrs[count];
			
			i = max;
			p[--j] = i;
			while (R[i] >= min) {
				i = R[i];
				p[--j] = i;
			}
			
			scores[count] = score;
			
			rectXI[count] = xi[max];
			rectYI[count] = yi[max];
			rectXS[count] = xs[min];
			rectXE[count] = xe[max];
			rectYS[count] = ys[min];
			rectYE[count] = ye[max];
			
			count++;
		}
	}
	
	Free(A);
	Free(S);
	Free(R);
	Free(O);
	Free(rectXS);
	Free(rectXE);
	Free(rectYS);
	Free(rectYE);
	Free(rectXI);
	Free(rectYI);
	
	SEXP ret, chains, chain, cs;
	int *pchain;
	double *pcs;
	PROTECT(cs = allocVector(REALSXP, count));
	pcs = REAL(cs);
	PROTECT(chains = allocVector(VECSXP, count));
	for (i = 0; i < count; i++) {
		PROTECT(chain = allocVector(INTSXP, lens[i]));
		pchain = INTEGER(chain);
		
		p = ptrs[i];
		for (j = 0; j < lens[i]; j++)
			pchain[j] = p[j] + 1;
		
		Free(p);
		
		SET_VECTOR_ELT(chains, i, chain);
		UNPROTECT(1);
		
		pcs[i] = scores[i];
	}
	
	Free(ptrs);
	Free(lens);
	Free(scores);
	
	PROTECT(ret = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret, 0, chains);
	SET_VECTOR_ELT(ret, 1, cs);
	
	UNPROTECT(3);
	
	return ret;
}

int extend(const Chars_holder *S1, const Chars_holder *S2, int *s1_start, int *s2_start, int s1_dir, int s2_dir, int bound, int rev, double drop) {
	int i, v1, v2, score = 0, max = 0, max_i = 0;
	const char *s1, *s2;
	
	for (i = 1, s1 = (S1->ptr + *s1_start + s1_dir - 1), s2 = (S2->ptr + *s2_start + s2_dir - 1);
		i <= bound;
		i++, s1 += s1_dir, s2 += s2_dir) {
		if ((*s1) & 0x01) {
			v1 = 0;
		} else if ((*s1) & 0x02) {
			v1 = 1;
		} else if ((*s1) & 0x04) {
			v1 = 2;
		} else { // assume 'T'
			v1 = 3;
		}
		if (rev) { // reverse complement
			if ((*s2) & 0x01) {
				v2 = 3;
			} else if ((*s2) & 0x02) {
				v2 = 2;
			} else if ((*s2) & 0x04) {
				v2 = 1;
			} else { // assume 'A'
				v2 = 0;
			}
		} else {
			if ((*s2) & 0x01) {
				v2 = 0;
			} else if ((*s2) & 0x02) {
				v2 = 1;
			} else if ((*s2) & 0x04) {
				v2 = 2;
			} else { // assume 'T'
				v2 = 3;
			}
		}
		
		if (v1==v2) { // perfect match
			score += 2;
		} else if ((v1==0 && v2==2) || // A/G
			(v1==2 && v2==0) || // G/A
			(v1==1 && v2==3) || // C/T
			(v1==3 && v2==1)) { // T/C
			// transition
			score -= 1;
		} else { // transversion
			score -= 2;
		}
		
		if (score > max) {
			max_i = i;
			max = score;
		} else if ((double)score < ((double)max + drop)) {
			break;
		}
	}
	
	//Rprintf("\nextend = %d score = %d max = %d rev = %d", max_i, score, max, rev);
	*s1_start = *s1_start + s1_dir*max_i;
	*s2_start = *s2_start + s2_dir*max_i;
	
	return(max);
}

SEXP extendSegments(SEXP X, SEXP W1, SEXP W2, SEXP S1, SEXP S2, SEXP O1P, SEXP O1N, SEXP O2P, SEXP O2N, SEXP S, SEXP maxDrop, SEXP INDEX1, SEXP INDEX2)
{
	int j, b, b1, b2, d, strand, *begin, score;
	int *x = INTEGER(X); // results matrix
	int *w1 = INTEGER(W1); // widths of sequence set 1
	int *w2 = INTEGER(W2); // widths of sequence set 2
	int *o1p = INTEGER(O1P); // order of previous starts
	int *o1n = INTEGER(O1N); // order of next starts
	int *o2p = INTEGER(O2P); // order of previous starts
	int *o2n = INTEGER(O2N); // order of next starts
	int *s = INTEGER(S); // subset (rows) of results
	int *index1 = INTEGER(INDEX1); // indices in sequence set 1
	int *index2 = INTEGER(INDEX2); // indices in sequence set 2
	int mD = asReal(maxDrop); // max drop score
	int l = length(S);
	SEXP dim = getAttrib(X, R_DimSymbol);
	int n = INTEGER(dim)[0]; // total rows in results
	
	XStringSet_holder s1_set, s2_set;
	Chars_holder s1, s2;
	s1_set = hold_XStringSet(S1);
	s2_set = hold_XStringSet(S2);
	
	for (j = 0; j < l; j++) {
		//Rprintf("\ns[j] = %d index1 = %d index2 = %d", s[j], index1[s[j]], index2[s[j]]);
		s1 = get_elt_from_XStringSet_holder(&s1_set, index1[s[j]]);
		s2 = get_elt_from_XStringSet_holder(&s2_set, index2[s[j]]);
		
		strand = x[2*n + s[j]];
		
		// extend left of start1
		if (o1p[j]==NA_INTEGER) {
			b1 = x[4*n + s[j]] - 1; // start1 - 1
		} else {
			b1 = x[4*n + s[j]] - x[6*n + s[o1p[j]]] - 1; // start - end - 1
		}
		if (strand) { // strand==1
			d = 1;
			begin = x + 7*n + s[j];
			// extend right of end2
			if (o2n[j]==NA_INTEGER) {
				b2 = w2[x[n + s[j]] - 1] - *begin; // width - end
			} else {
				b2 = x[5*n + s[o2n[j]]] - *begin - 1; // start - end - 1
			}
		} else {
			d = -1;
			begin = x + 5*n + s[j];
			// extend left of start2
			if (o2p[j]==NA_INTEGER) {
				b2 = *begin - 1; // start - 1
			} else {
				b2 = *begin - x[7*n + s[o2p[j]]] - 1; // start - end - 1
			}
		}
		b = b1 < b2 ? b1 : b2; // max allowable extension
		if (b > 0) {
			score = extend(&s1, &s2, x + 4*n + s[j], begin, -1, d, b, strand, mD);
			x[3*n + s[j]] += score;
		}
		//Rprintf("\nj = %d b1 = %d b2 = %d start1 = %d start2 = %d", j, b1, b2, x[4*n + s[j]], *begin);
		
		// extend right of end1
		if (o1n[j]==NA_INTEGER) {
			b1 = w1[x[s[j]] - 1] - x[6*n + s[j]]; // width - end
		} else {
			b1 = x[4*n + s[o1n[j]]] - x[6*n + s[j]] - 1; // start - end - 1
		}
		if (strand) { // strand==1
			d = -1;
			begin = x + 5*n + s[j];
			// extend left of end2
			if (o2p[j]==NA_INTEGER) {
				b2 = *begin - 1; // start - 1
			} else {
				b2 = *begin - x[7*n + s[o2p[j]]] - 1; // start - end - 1
			}
		} else {
			d = 1;
			begin = x + 7*n + s[j];
			// extend right of start2
			if (o2n[j]==NA_INTEGER) {
				b2 = w2[x[n + s[j]] - 1] - *begin; // width - end
			} else {
				b2 = x[5*n + s[o2n[j]]] - *begin - 1; // start - end - 1
			}
		}
		b = b1 < b2 ? b1 : b2; // max allowable extension
		if (b > 0) {
			score = extend(&s1, &s2, x + 6*n + s[j], begin, 1, d, b, strand, mD);
			x[3*n + s[j]] += score;
		}
		//Rprintf("\nj = %d b1 = %d b2 = %d start1 = %d start2 = %d", j, b1, b2, x[6*n + s[j]], *begin);
	}
	
	return R_NilValue;
}

// changes repeat regions to NAs
SEXP maskRepeats(SEXP e, SEXP size, SEXP minL, SEXP maxL, SEXP totL)
{	
	int i, p, j, k;
	int *x = INTEGER(e); // enumerated sequence
	int l = length(e);
	int n = asInteger(size);
	int l1 = asInteger(minL); // min period
	int l2 = asInteger(maxL); // max period
	int l3 = asInteger(totL); // min length of repeat
	
	i = 0; // current position
	while (i < (l - l2)) {
		if (x[i]!=NA_INTEGER) {
			for (p = l1; p <= l2; p++) { // periodicity
				if (x[i]==x[i + p]) { // repeat
					j = i + 1;
					
					while (j < (l - p)) {
						if (x[j]!=x[j + p])
							break;
						j++;
					}
					
					if ((j - i + n) > p && // continuous repeat
						(j + p - i + n) > l3) {
						for (k = i; k <= (j + p - 1); k++)
							x[k] = NA_INTEGER;
						i = k - 1;
						break;
					}
				}
			}
		}
		i++;
	}
	
	return R_NilValue;
}
