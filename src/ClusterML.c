/****************************************************************************
 *                        Cluster Maximum Likelihood                        *
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

static void L_known(const char *p, double *Ls, int *length)
{
	switch (*p) {
		case 1: // A
			*(Ls) = 1;
			break;
		case 2: // C
			*(Ls + *length) = 1;
			break;
		case 3: // M
			*(Ls) = .5; *(Ls + *length) = .5; // AC
			break;
		case 4: // G
			*(Ls + 2*(*length)) = 1;
			break;
		case 5: // R
			*(Ls) = .5; *(Ls + 2*(*length)) = .5; // AG
			break;
		case 6: // S
			*(Ls + *length) = .5; *(Ls + 2*(*length)) = .5; // CG
			break;
		case 7: // V
			*(Ls) = (double)1/3; *(Ls + *length) = (double)1/3; *(Ls + 2*(*length)) = (double)1/3; // ACG
			break;
		case 8: // T
			*(Ls + 3*(*length)) = 1;
			break;
		case 9: // W
			*(Ls) = .5; *(Ls + 3*(*length)) = .5; // AT
			break;
		case 10: // Y
			*(Ls + *length) = .5; *(Ls + 3*(*length)) = .5; // CT
			break;
		case 11: // H
			*(Ls) = (double)1/3; *(Ls + *length) = (double)1/3; *(Ls + 3*(*length)) = (double)1/3; // ACT
			break;
		case 12: // K
			*(Ls + 2*(*length)) = .5; *(Ls + 3*(*length)) = .5; // GT
			break;
		case 13: // D
			*(Ls) = (double)1/3; *(Ls + 2*(*length)) = (double)1/3; *(Ls + 3*(*length)) = (double)1/3; // AGT
			break;
		case 14: // B
			*(Ls + *length) = (double)1/3; *(Ls + 2*(*length)) = (double)1/3; *(Ls + 3*(*length)) = (double)1/3; // CGT
			break;
		case 15: // N
			//*(Ls) = .25; *(Ls + *length) = .25; *(Ls + 2*(*length)) = .25; *(Ls + 3*(*length)) = .25; // ACGT
			break;
		case 16: // -
			//*(Ls) = .25; *(Ls + *length) = .25; *(Ls + 2*(*length)) = .25; *(Ls + 3*(*length)) = .25; // all possible
			break;
		case 32: // +
			//*(Ls) = .25; *(Ls + *length) = .25; *(Ls + 2*(*length)) = .25; *(Ls + 3*(*length)) = .25; // all possible
			break;
		default:
			error("not DNA!");
			break;
	}
}

static void L_unknown(double *Ls, int offset, int j, int row, int length, double *P1, double *P2)
{
	double L1, L2;
	
	if (*(Ls + 0*length + row)!=0 ||
		*(Ls + 1*length + row)!=0 ||
		*(Ls + 2*length + row)!=0 ||
		*(Ls + 3*length + row)!=0) {
		if (*(Ls + 4*length + row)!=0 ||
			*(Ls + 5*length + row)!=0 ||
			*(Ls + 6*length + row)!=0 ||
			*(Ls + 7*length + row)!=0) {
			// neither branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(*(Ls + 0*length + row)); // Laa
			L1 += *(P1 + 10)*(*(Ls + 1*length + row)); // Lca
			L1 += *(P1 + 11)*(*(Ls + 2*length + row)); // Lga
			L1 += *(P1 + 12)*(*(Ls + 3*length + row)); // Lta
			L2 = *(P2 + 0)*(*(Ls + 4*length + row)); // Laa
			L2 += *(P2 + 10)*(*(Ls + 5*length + row)); // Lca
			L2 += *(P2 + 11)*(*(Ls + 6*length + row)); // Lga
			L2 += *(P2 + 12)*(*(Ls + 7*length + row)); // Lta
			*(Ls + offset*length + j) = L1*L2;
			
			// L(C)
			L1 = *(P1 + 4)*(*(Ls + 1*length + row)); // Lcc
			L1 += *(P1 + 1)*(*(Ls + 0*length + row)); // Lac
			L1 += *(P1 + 13)*(*(Ls + 2*length + row)); // Lgc
			L1 += *(P1 + 14)*(*(Ls + 3*length + row)); // Ltc
			L2 = *(P2 + 4)*(*(Ls + 5*length + row)); // Lcc
			L2 += *(P2 + 1)*(*(Ls + 4*length + row)); // Lac
			L2 += *(P2 + 13)*(*(Ls + 6*length + row)); // Lgc
			L2 += *(P2 + 14)*(*(Ls + 7*length + row)); // Ltc
			*(Ls + (offset + 1)*length + j) = L1*L2;
			
			// L(G)
			L1 = *(P1 + 7)*(*(Ls + 2*length + row)); // Lgg
			L1 += *(P1 + 2)*(*(Ls + 0*length + row)); // Lag
			L1 += *(P1 + 5)*(*(Ls + 1*length + row)); // Lcg
			L1 += *(P1 + 15)*(*(Ls + 3*length + row)); // Ltg
			L2 = *(P2 + 7)*(*(Ls + 6*length + row)); // Lgg
			L2 += *(P2 + 2)*(*(Ls + 4*length + row)); // Lag
			L2 += *(P2 + 5)*(*(Ls + 5*length + row)); // Lcg
			L2 += *(P2 + 15)*(*(Ls + 7*length + row)); // Ltg
			*(Ls + (offset + 2)*length + j) = L1*L2;
			
			// L(T)
			L1 = *(P1 + 9)*(*(Ls + 3*length + row)); // Ltt
			L1 += *(P1 + 3)*(*(Ls + 0*length + row)); // Lat
			L1 += *(P1 + 6)*(*(Ls + 1*length + row)); // Lct
			L1 += *(P1 + 8)*(*(Ls + 2*length + row)); // Lgt
			L2 = *(P2 + 9)*(*(Ls + 7*length + row)); // Ltt
			L2 += *(P2 + 3)*(*(Ls + 4*length + row)); // Lat
			L2 += *(P2 + 6)*(*(Ls + 5*length + row)); // Lct
			L2 += *(P2 + 8)*(*(Ls + 6*length + row)); // Lgt
			*(Ls + (offset + 3)*length + j) = L1*L2;
		} else {
			// second branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(*(Ls + 0*length + row)); // Laa
			L1 += *(P1 + 10)*(*(Ls + 1*length + row)); // Lca
			L1 += *(P1 + 11)*(*(Ls + 2*length + row)); // Lga
			L1 += *(P1 + 12)*(*(Ls + 3*length + row)); // Lta
			*(Ls + offset*length + j) = L1;
			
			// L(C)
			L1 = *(P1 + 4)*(*(Ls + 1*length + row)); // Lcc
			L1 += *(P1 + 1)*(*(Ls + 0*length + row)); // Lac
			L1 += *(P1 + 13)*(*(Ls + 2*length + row)); // Lgc
			L1 += *(P1 + 14)*(*(Ls + 3*length + row)); // Ltc
			*(Ls + (offset + 1)*length + j) = L1;
			
			// L(G)
			L1 = *(P1 + 7)*(*(Ls + 2*length + row)); // Lgg
			L1 += *(P1 + 2)*(*(Ls + 0*length + row)); // Lag
			L1 += *(P1 + 5)*(*(Ls + 1*length + row)); // Lcg
			L1 += *(P1 + 15)*(*(Ls + 3*length + row)); // Ltg
			*(Ls + (offset + 2)*length + j) = L1;
			
			// L(T)
			L1 = *(P1 + 9)*(*(Ls + 3*length + row)); // Ltt
			L1 += *(P1 + 3)*(*(Ls + 0*length + row)); // Lat
			L1 += *(P1 + 6)*(*(Ls + 1*length + row)); // Lct
			L1 += *(P1 + 8)*(*(Ls + 2*length + row)); // Lgt
			*(Ls + (offset + 3)*length + j) = L1;
		}
	} else if (*(Ls + 0*length + row)==0 &&
			   *(Ls + 1*length + row)==0 &&
			   *(Ls + 2*length + row)==0 &&
			   *(Ls + 3*length + row)==0) {
		if (*(Ls + 4*length + row)!=0 ||
			*(Ls + 5*length + row)!=0 ||
			*(Ls + 6*length + row)!=0 ||
			*(Ls + 7*length + row)!=0) {
			// first branch can be disregarded
			
			// L(A)
			L2 = *(P2 + 0)*(*(Ls + 4*length + row)); // Laa
			L2 += *(P2 + 10)*(*(Ls + 5*length + row)); // Lca
			L2 += *(P2 + 11)*(*(Ls + 6*length + row)); // Lga
			L2 += *(P2 + 12)*(*(Ls + 7*length + row)); // Lta
			*(Ls + offset*length + j) = L2;
			
			// L(C)
			L2 = *(P2 + 4)*(*(Ls + 5*length + row)); // Lcc
			L2 += *(P2 + 1)*(*(Ls + 4*length + row)); // Lac
			L2 += *(P2 + 13)*(*(Ls + 6*length + row)); // Lgc
			L2 += *(P2 + 14)*(*(Ls + 7*length + row)); // Ltc
			*(Ls + (offset + 1)*length + j) = L2;
			
			// L(G)
			L2 = *(P2 + 7)*(*(Ls + 6*length + row)); // Lgg
			L2 += *(P2 + 2)*(*(Ls + 4*length + row)); // Lag
			L2 += *(P2 + 5)*(*(Ls + 5*length + row)); // Lcg
			L2 += *(P2 + 15)*(*(Ls + 7*length + row)); // Ltg
			*(Ls + (offset + 2)*length + j) = L2;
			
			// L(T)
			L2 = *(P2 + 9)*(*(Ls + 7*length + row)); // Ltt
			L2 += *(P2 + 3)*(*(Ls + 4*length + row)); // Lat
			L2 += *(P2 + 6)*(*(Ls + 5*length + row)); // Lct
			L2 += *(P2 + 8)*(*(Ls + 6*length + row)); // Lgt
			*(Ls + (offset + 3)*length + j) = L2;
		}
	}

}

static void ProbChange(double *m, double *P, double v)
{
	/*
	+-                                                         -+
	|       G                           G                       |
	|  #6 + -- + #9,      #12,     #4 - -- + #10,      #14      |
	|       #2                          #2                      |
	|                                                           |
	|                     T                                 T   |
	|       #13,     #5 + -- + #8,      #11,      #3 + #7 - --  |
	|                     #1                                #1  |
	|                                                           |
	|       A                           A                       |
	|  #6 - -- + #9,      #12,     #4 + -- + #10,      #14      |
	|       #2                          #2                      |
	|                                                           |
	|                     C                                 C   |
	|       #13,     #5 - -- + #8,      #11,      #3 + #7 + --  |
	|                     #1                                #1  |
	+-                                                         -+
	
	where
	
	#1 = exp(A t u + G t u + T k2 t u + C k2 t u) (C + T)
		
	#2 = exp(C t u + T t u + A k1 t u + G k1 t u) (A + G)
		
			T
	#3 = -------------
		A + C + G + T
		
			G
	#4 = -------------
		A + C + G + T
		
			C
	#5 = -------------
		A + C + G + T
		
			A
	#6 = -------------
		A + C + G + T
		
		A T + G T
	#7 = ---------
			#15
		
		C (A + G)
	#8 = ---------
			#15
		
		A (C + T)
	#9 = ---------
			#16
		
		G (C + T)
	#10 = ---------
			#16
		
			G #17
	#11 = - -------------
			A + C + G + T
		
			C #17
	#12 = - -------------
			A + C + G + T
		
			A #17
	#13 = - -------------
			A + C + G + T
		
			T #17
	#14 = - -------------
			A + C + G + T
		
				2    2
	#15 = #18 (C  + T  + A C + C G + A T + 2 C T + G T)
		
	#16 = #18 (A + G) (A + C + G + T)
		
								1
	#17 = ------------------------------------------- - 1
		exp(A t u) exp(C t u) exp(G t u) exp(T t u)
		
	#18 = exp(A t u + C t u + G t u + T t u)
	*/
	
	double e1, e2, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18;
	double A = *(m), C = *(m + 1), G = *(m + 2), T = *(m + 3), k1 = *(m + 4), k2 = *(m + 5);
	e1 = exp((A + G + T*k2 + C*k2)*v)*(C + T);
	e2 = exp((C + T + A*k1 + G*k1)*v)*(A + G);
	e18 = exp(v);
	e15 = e18*(C*C + T*T + A*C + C*G + A*T + 2*C*T + G*T);
	e7 = (A*T + G*T)/e15;
	e8 = C*(A + G)/e15;
	e16 = e18*(A + G);
	e9 = A*(C + T)/e16;
	e10 = G*(C + T)/e16;
	e17 = 1/(exp(A*v)*exp(C*v)*exp(G*v)*exp(T*v)) - 1;
	e11 = -1*G*e17;
	e12 = -1*C*e17;
	e13 = -1*A*e17;
	e14 = -1*T*e17;
	
	*P = A + G/e2 + e9; // Paa
	*(P + 1) = e12; // Pac
	*(P + 2) = G - G/e2 + e10; // Pag
	*(P + 3) = e14; // Pat
	*(P + 4) = C + T/e1 + e8; // Pcc
	*(P + 5) = e11; // Pcg
	*(P + 6) = T + e7 - T/e1; // Pct
	*(P + 7) = G + A/e2 + e10; // Pgg
	*(P + 8) = e14; // Pgt
	*(P + 9) = T + e7 + C/e1; // Ptt
	*(P + 10) = e13; // Pca
	*(P + 11) = A - A/e2 + e9; // Pga
	*(P + 12) = e13; // Pta
	*(P + 13) = e12; // Pgc
	*(P + 14) = C - C/e1 + e8; // Ptc
	*(P + 15) = e11; // Ptg
}

SEXP clusterML(SEXP x, SEXP y, SEXP model, SEXP verbose, SEXP pBar)
{
	// input is the output tree from clusterNJ.c
	
	// step 1:  figure out widths of each sequence
	// step 2:  for each position calculate the likelihood(L)
	// step 3:  for each node calculate L at the node
	//      3a: if node is external (leaf) then L(base) is known
	//          if base is ambigous then L(base) is split up
	//          gaps and N are disregarded as unknowns L() = 0
	//      3b: if node is internal then L must be calculated
	//          L = probaility of changing over that branch length
	//              times L of previous node
	//              multiplied together for both nodes
	//          if branch length is zero then L is same as previous
	//          if L of previous node is zero then it is disregarded
	// step 4:  L of the whole tree must be calculated
	//          L = probability at final node for each base
	//              multipled by that bases' frequency
	//          since L is a very small number -Ln(L) is used
	//          -LnL = -1*sum(Ln(base)*frequency)
	
	// initialize variables
	cachedXStringSet y_set;
	cachedCharSeq y_i;
	y_set = cache_XStringSet(y);
	int length = get_cachedXStringSet_length(&y_set);
	int i, j, k, numRates, row;
	double *T = REAL(x); // Tree Topology
	double *m = REAL(model); // Substitution Model [%A %C %G %T k1 k2]
	int *widths = (int *) R_alloc(length, sizeof(int));
	int v = asLogical(verbose);
	
	// calculate a vector of sequence lengths
	int maxWidth = 0;
	for (i = 0; i < length; i++) {
		y_i = get_cachedXStringSet_elt(&y_set, i);
		*(widths + i) = y_i.length;
		if (*(widths + i) > maxWidth) {
			maxWidth = *(widths + i);
		}
	}
	
	double *sumL = Calloc(maxWidth, double);
	numRates = (length(model) - 6)/2; // number of bins for the gamma distribution
	for (k = 0; k < numRates; k++) { // for each bin of the gamma distribution determined by alfa
		// P = [Paa Pac Pag Pat Pcc Pcg Pct Pgg Pgt Ptt Pca Pga Pta Pgc Ptc Ptg]
		double *P = Calloc((length - 1)*32, double); // initialized to zero
		#pragma omp parallel for private(i) schedule(guided)
		for (i = 0; i < (length - 1); i++) {
			ProbChange(m, (P + i*32), T[3*(length - 1) + i] * *(m + k + 6));
			ProbChange(m, (P + i*32 + 16), T[4*(length - 1) + i] * *(m + k + 6));
		}
		
		#pragma omp parallel for private(i,j,y_i,row,sumL) schedule(guided)
		for (i = 0; i < maxWidth; i++) { // for each position
			//double *Ls = (double *) R_alloc(length*8, sizeof(double));
			double *Ls = Calloc(length*8, double); // initialized to zero
			
			for (j = 0; j < (length - 1); j++) { // for each node
				// if first branch is leaf then its base L is 1
				if ((int)T[6*(length - 1) + j] < 0) { // first branch is a leaf
					if (*(widths + (-1*(int)T[6*(length - 1) + j] - 1)) > i) { // position exist in this sequence
						y_i = get_cachedXStringSet_elt(&y_set, (-1*(int)T[6*(length - 1) + j] - 1));
						L_known(&y_i.seq[i], (Ls + 0*length + j), &length);
					}
				}
				// if second branch is leaf then its base L is 1
				if ((int)T[7*(length - 1) + j] < 0) { // second branch is a leaf
					if (*(widths + (-1*(int)T[7*(length - 1) + j] - 1)) > i) { // position exist in this sequence
						y_i = get_cachedXStringSet_elt(&y_set, (-1*(int)T[7*(length - 1) + j] - 1));
						L_known(&y_i.seq[i], (Ls + 4*length + j), &length);
					}
				}
				// if the first branch is a node then L must be calculated
				if ((int)T[6*(length - 1) + j] > 0) { // first branch is a node
					// L is probability(branch lengths) * L(previous nodes)
					row = (int)T[6*(length - 1) + j] - 1;
					
					L_unknown(Ls, 0, j, row, length, (P + row*32), (P + row*32 + 16));
				}
				// if the second branch is a node then L must be calculated
				if ((int)T[7*(length - 1) + j] > 0) { // first branch is a node
					// L is probability(branch lengths) * L(previous nodes)
					row = (int)T[7*(length - 1) + j] - 1;
					
					L_unknown(Ls, 4, j, row, length, (P + row*32), (P + row*32 + 16));
				}
			}
			
			// calculate Likelihood of each base at final node
			row = length - 2;  // final node: j = length - 1
			
			L_unknown(Ls, 0, j, row, length, (P + row*32), (P + row*32 + 16));
			/* // Prints last common ancestor:
			if (*(Ls + 0*length + j)==0 &&
				*(Ls + 1*length + j)==0 &&
				*(Ls + 2*length + j)==0 &&
				*(Ls + 3*length + j)==0) {
				if (i==0) Rprintf("\n");
				Rprintf("-");
				continue;
			} else if (*(Ls + 0*length + j) > *(Ls + 1*length + j) &&
				*(Ls + 0*length + j) > *(Ls + 2*length + j) &&
				*(Ls + 0*length + j) > *(Ls + 3*length + j)) {
				Rprintf("A");
			} else if (*(Ls + 1*length + j) > *(Ls + 2*length + j) &&
					   *(Ls + 1*length + row) > *(Ls + 3*length + j)) {
				Rprintf("C");
			} else if (*(Ls + 2*length + j) > *(Ls + 3*length + j)) {
				Rprintf("G");
			} else {
				Rprintf("T");
			}
			*/
			
			// sum likelihoods of tree for every position
			*(sumL + i) += (*(m) * *(Ls + 0*length + j) +
				 *(m + 1) * *(Ls + 1*length + j) +
				 *(m + 2) * *(Ls + 2*length + j) +
				*(m + 3) * *(Ls + 3*length + j)) * *(m + numRates + k + 6);
			Free(Ls);
		}
	}
	
	double LnL = 0;
	for (i = 0; i < maxWidth; i++) {
		if (*(sumL + i) > 0) {
			LnL -= log(*(sumL + i));
		}
	}
	Free(sumL);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 1));
	double *rans = REAL(ans);
	rans[0] = LnL;
	UNPROTECT(1);
	
	return ans;
}
