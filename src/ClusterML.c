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
#include <omp.h>

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

static void L_unknown(double *Ls, int offset, int j, int row, int length, double Pii1, double Pij1, double Pii2, double Pij2)
{
	double L1, L2;
	
	if (Pij1==0) {
		// first branch is zero length
		*(Ls + offset*length + j) = *(Ls + 0*length + row);
		*(Ls + (offset + 1)*length + j) = *(Ls + 1*length + row);
		*(Ls + (offset + 2)*length + j) = *(Ls + 2*length + row);
		*(Ls + (offset + 3)*length + j) = *(Ls + 3*length + row);
	} else if (Pij2==0) {
		// second branch is zero length
		*(Ls + offset*length + j) = *(Ls + 4*length + row);
		*(Ls + (offset + 1)*length + j) = *(Ls + 5*length + row);
		*(Ls + (offset + 2)*length + j) = *(Ls + 6*length + row);
		*(Ls + (offset + 3)*length + j) = *(Ls + 7*length + row);
	} else if (*(Ls + 0*length + row)!=0 ||
		*(Ls + 1*length + row)!=0 ||
		*(Ls + 2*length + row)!=0 ||
		*(Ls + 3*length + row)!=0) {
		if (*(Ls + 4*length + row)!=0 ||
			*(Ls + 5*length + row)!=0 ||
			*(Ls + 6*length + row)!=0 ||
			*(Ls + 7*length + row)!=0) {
			// neither branch can be disregarded
			
			// L(A)
			L1 = Pii1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 1*length + row));
			L1 += Pij1*(*(Ls + 2*length + row));
			L1 += Pij1*(*(Ls + 3*length + row));
			L2 = Pii2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 5*length + row));
			L2 += Pij2*(*(Ls + 6*length + row));
			L2 += Pij2*(*(Ls + 7*length + row));
			*(Ls + offset*length + j) = L1*L2;
			
			// L(C)
			L1 = Pii1*(*(Ls + 1*length + row));
			L1 += Pij1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 2*length + row));
			L1 += Pij1*(*(Ls + 3*length + row));
			L2 = Pii2*(*(Ls + 1*length + row));
			L2 += Pij2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 6*length + row));
			L2 += Pij2*(*(Ls + 7*length + row));
			*(Ls + (offset + 1)*length + j) = L1*L2;
			
			// L(G)
			L1 = Pii1*(*(Ls + 2*length + row));
			L1 += Pij1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 1*length + row));
			L1 += Pij1*(*(Ls + 3*length + row));
			L2 = Pii2*(*(Ls + 6*length + row));
			L2 += Pij2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 5*length + row));
			L2 += Pij2*(*(Ls + 7*length + row));
			*(Ls + (offset + 2)*length + j) = L1*L2;
			
			// L(T)
			L1 = Pii1*(*(Ls + 3*length + row));
			L1 += Pij1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 1*length + row));	
			L1 += Pij1*(*(Ls + 2*length + row));
			L2 = Pii2*(*(Ls + 7*length + row));
			L2 += Pij2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 5*length + row));
			L2 += Pij2*(*(Ls + 6*length + row));
			*(Ls + (offset + 3)*length + j) = L1*L2;
		} else {
			// second branch can be disregarded
			
			// L(A)
			L1 = Pii1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 1*length + row));
			L1 += Pij1*(*(Ls + 2*length + row));
			L1 += Pij1*(*(Ls + 3*length + row));
			*(Ls + offset*length + j) = L1;
			
			// L(C)
			L1 = Pii1*(*(Ls + 1*length + row));
			L1 += Pij1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 2*length + row));
			L1 += Pij1*(*(Ls + 3*length + row));
			*(Ls + (offset + 1)*length + j) = L1;
			
			// L(G)
			L1 = Pii1*(*(Ls + 2*length + row));
			L1 += Pij1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 1*length + row));
			L1 += Pij1*(*(Ls + 3*length + row));
			*(Ls + (offset + 2)*length + j) = L1;
			
			// L(T)
			L1 = Pii1*(*(Ls + 3*length + row));
			L1 += Pij1*(*(Ls + 0*length + row));
			L1 += Pij1*(*(Ls + 1*length + row));	
			L1 += Pij1*(*(Ls + 2*length + row));
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
			L2 = Pii2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 5*length + row));
			L2 += Pij2*(*(Ls + 6*length + row));
			L2 += Pij2*(*(Ls + 7*length + row));
			*(Ls + offset*length + j) = L2;
			
			// L(C)
			L2 = Pii2*(*(Ls + 1*length + row));
			L2 += Pij2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 6*length + row));
			L2 += Pij2*(*(Ls + 7*length + row));
			*(Ls + (offset + 1)*length + j) = L2;
			
			// L(G)
			L2 = Pii2*(*(Ls + 6*length + row));
			L2 += Pij2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 5*length + row));
			L2 += Pij2*(*(Ls + 7*length + row));
			*(Ls + (offset + 2)*length + j) = L2;
			
			// L(T)
			L2 = Pii2*(*(Ls + 7*length + row));
			L2 += Pij2*(*(Ls + 4*length + row));
			L2 += Pij2*(*(Ls + 5*length + row));
			L2 += Pij2*(*(Ls + 6*length + row));
			*(Ls + (offset + 3)*length + j) = L2;
		}
	}

}

SEXP clusterML(SEXP x, SEXP y, SEXP cutoff, SEXP verbose, SEXP pBar)
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
	int i, j, o, p, row;
	double *T = REAL(x); // Tree Topology
	double Pii1, Pij1, Pii2, Pij2, LnL = 0;
	double *Ls = (double *) R_alloc(length*8, sizeof(double));
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
	
	for (i = 0; i < maxWidth; i++) { // for each position
		// zero Ls matrix
		for (o = 0; o < length; o++) {
			for (p = 0; p < 8; p++) {
				*(Ls + p*length + o) = 0;
			}
		}
		
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
				
				Pii1 = .25 + .75*exp(-1*(T[3*(length - 1) + row])); // probability same base
				Pij1 = .25 - .25*exp(-1*(T[3*(length - 1) + row])); // probability diff base
				Pii2 = .25 + .75*exp(-1*(T[4*(length - 1) + row])); // probability same base
				Pij2 = .25 - .25*exp(-1*(T[4*(length - 1) + row])); // probability diff base
				
				L_unknown(Ls, 0, j, row, length, Pii1, Pij1, Pii2, Pij2);
			}
			// if the second branch is a node then L must be calculated
			if ((int)T[7*(length - 1) + j] > 0) { // first branch is a node
				// L is probability(branch lengths) * L(previous nodes)
				row = (int)T[7*(length - 1) + j] - 1;
				
				Pii1 = .25 + .75*exp(-1*(T[3*(length - 1) + row])); // probability same base
				Pij1 = .25 - .25*exp(-1*(T[3*(length - 1) + row])); // probability diff base
				Pii2 = .25 + .75*exp(-1*(T[4*(length - 1) + row])); // probability same base
				Pij2 = .25 - .25*exp(-1*(T[4*(length - 1) + row])); // probability diff base
				
				L_unknown(Ls, 4, j, row, length, Pii1, Pij1, Pii2, Pij2);
			}
		}
		
		// calculate Likelihood of each base at final node
		row = length - 2;  // final node: j = length - 1
		
		Pii1 = .25 + .75*exp(-1*(T[3*(length - 1) + row])); // probability same base
		Pij1 = .25 - .25*exp(-1*(T[3*(length - 1) + row])); // probability diff base
		Pii2 = .25 + .75*exp(-1*(T[4*(length - 1) + row])); // probability same base
		Pij2 = .25 - .25*exp(-1*(T[4*(length - 1) + row])); // probability diff base
		
		L_unknown(Ls, 0, j, row, length, Pii1, Pij1, Pii2, Pij2);
		
		/*
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
		
		// sum -LnL of tree for every position
		LnL -= log(.25*(*(Ls + 0*length + j) +
						*(Ls + 1*length + j) +
						*(Ls + 2*length + j) +
						*(Ls + 3*length + j)));
	}
	
	if (v)
		Rprintf("\n-LnL = %d", (int)LnL);
	
	return(x);
}
