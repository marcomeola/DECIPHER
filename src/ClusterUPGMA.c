/****************************************************************************
 *                           Cluster Ultrametric                            *
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


void binUPGMA(double *rans, int i, int clusterNumber, double maxHeight, int length) {
	// NOTE:  cluster number is going to be the index i + 1
	if (rans[8*(length - 1) + i] == 0 || rans[9*(length - 1) + i] == 0) { // cluster number unassigned
		// then assign the node a cluster number
		if (rans[8*(length - 1) + i] == 0) { // the first fork is unassigned
			if (rans[6*(length - 1) + i] < 0) { // the first fork is a leaf
				rans[8*(length - 1) + i] = clusterNumber;
			} else { // the first fork is a branch
				rans[8*(length - 1) + i] = -1;
			}
		}
		if (rans[9*(length - 1) + i] == 0) { // the second fork is unassigned
			if (rans[7*(length - 1) + i] < 0) { // the second fork is a leaf
				rans[9*(length - 1) + i] = clusterNumber;
			} else { // the second fork is a branch
				rans[9*(length - 1) + i] = -1;
			}
		}
		
		// keep going up the tree
		for (int j = i + 1; j < length - 1; j++) {
			if (rans[6*(length - 1) + j] == (i + 1) || rans[7*(length - 1) + j] == (i + 1)) { // node is merged again
				if (rans[5*(length - 1) + j] <= maxHeight) { // if the node is within reach
					// then assign the same clusterNumber to this node
					binUPGMA(rans, j, clusterNumber, maxHeight, length);
					break;
				}
			}
		}
	}
	
	// follow the branches down the tree
	if (rans[6*(length - 1) + i] > 0) { // if first fork is a branch
		binUPGMA(rans, (int)(rans[6*(length - 1) + i] - 1), clusterNumber, maxHeight, length);
	}
	if (rans[7*(length - 1) + i] > 0) { // if second fork is a branch
		binUPGMA(rans, (int)(rans[7*(length - 1) + i] - 1), clusterNumber, maxHeight, length);
	}
}

//ans_start <- .Call("cluster", myDistMatrix, verbose, pBar, PACKAGE="DECIPHER")
SEXP clusterUPGMA(SEXP x, SEXP cutoff, SEXP method, SEXP verbose, SEXP pBar)
{	
	/*
	 * **** Input ****
	 *
	 *    Dist Matrix
	 *   A B C D E F G
	 * A 0
	 * B x 0
	 * C x x 0
	 * D x x x 0
	 * E x x x x 0
	 * F x x x x x 0
	 * G x x x x x x 0
	 */
	
	/*
	 * **********   Output (ans)   **********
	 *
	 * rans[0*(length - 1) + k] // 1st col = row merged
	 * rans[1*(length - 1) + k] // 2nd col = col merged
	 * rans[2*(length - 1) + k] // 3rd col = cluster number
	 * rans[3*(length - 1) + k] // 4th col = row branch length
	 * rans[4*(length - 1) + k] // 5th col = col branch length
	 * rans[5*(length - 1) + k] // 6th col = height of merger
	 * rans[6*(length - 1) + k] // 7th col = alternative numbering
	 * rans[7*(length - 1) + k] // 8th col = alternative numbering
	 * rans[8*(length - 1) + k] // 9th col = cutoff numbering
	 * rans[9*(length - 1) + k] // 10th col = cutoff numbering
	 */
	
	// distanceMatrix is a pointer to x
	// clusterNum hold the number of clusters
	// dMatrix is a manipulatable copy of distance matrix
	// *rans is a pointer to the output ans
	// cluster number is stored in the outer col, row of dMatrix
	// single leaves are negative, clusters are positive
	// Cluster numbering notes:
	// starting at lowest leaf, climb the tree until over maxHeight
	// maxHeight is merge height plus cutoff minus longest branch
	// do not go down the tree past floorHeight
	
	int i, j, k, clusterNum, size, minRow, minCol, index, minR, minC, met;
	int before, v, *rPercentComplete;
	double soFar, total, minHeight, *cut, *rans, *distanceMatrix, minH;
	SEXP ans, percentComplete, utilsPackage;
	
	// initialize variables
	clusterNum = 0; // increments with each new cluster
	size = sqrt(length(x)); // square distance matrix dimension
	const int length = size; // doesn't decrease in size
	PROTECT(ans = allocMatrix(REALSXP, (size - 1), 10));
	rans = REAL(ans);
	distanceMatrix = REAL(x);
	double *dMatrix = (double *) R_alloc(size*size, sizeof(double)); // final row & col contain cluster numbers
	double dTemp[size - 2], dist1, dist2;
	double cumHeight[size - 2];
	int clusterNums[size - 2];
	cut = REAL(cutoff);
	met = asInteger(method);
	v = asLogical(verbose);
	
	if (v) { // initialize progress variables
		soFar = 0;
		before = 0;
		total = (2*pow(length - 1,3) + 3*pow(length - 1,2) + length - 1)/6;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	// copy the lower triangle of the distance matrix (x)
	// into the lower triangle of the dMatrix
	for (i = 0; i < (size - 1); i++) { // length is the same as size but it doesn't change
		for (j = i; j < (size - 1); j++) {
			*(dMatrix + j*length + i) = *(distanceMatrix + i*size + j + 1);
		}
	}
	
	// initialize the order of groupings in the matrix
	for (i = 0; i < (size - 1); i++) {
		*(dMatrix + i*length + (size - 1)) = -(i + 2); // final column of dMatrix has cluster numbers
	}
	for (j = 0; j < (size - 1); j++) {
		*(dMatrix + (size - 1)*length + j) = -(j + 1); // final row of dMatrix has cluster numbers
	}
	
	// initialize rans to zero representing incomplete answer
	for (i = 0; i < (size - 1); i++) {
		for (j = 0; j < 10; j++) {
			*(rans + j*(length - 1) + i) = 0;
		}
	}
	
	// start the loop that goes from tree leaf to root
	for (k = 0; k < (length - 1); k++) {
		// calculate the Q matrix & find the smallest element in the Q matrix
		minRow = 0;
		minCol = 0;
		minHeight = 1e50;
		#pragma omp parallel for private(i,j,minR,minC,minH) schedule(guided)
		for (i = (size - 2); i >= 0; i--) {
			minH = minHeight;
			for (j = 0; j <= i; j++) {
				if (*(dMatrix + i*length + j) < minH) {
					minH = *(dMatrix + i*length + j);
					minR = i;
					minC = j;
				}
			}
			#pragma omp critical
			if (minH < minHeight) {
				minHeight = minH;
				minRow = minR;
				minCol = minC;
			}
		}
		
		// merge into a cluster
		rans[0*(length - 1) + k] = *(dMatrix + minRow*length + length - 1); // row merged
		rans[1*(length - 1) + k] = *(dMatrix + (length - 1)*length + minCol); // column merged
		
		// cluster
		if (((rans[0*(length - 1) + k] < 0) && (rans[1*(length - 1) + k] < 0)) ||
			((rans[0*(length - 1) + k] > 0) && (rans[1*(length - 1) + k] > 0))) {
			// merge into a new cluster
			clusterNum++;
			rans[2*(length - 1) + k] = clusterNum; // cluster formed
			// calculate both branch lengths
			rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2; // col
			rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2; // row
			
			// add longest path to cumulative height
			if ((rans[0*(length - 1) + k] < 0) && (rans[1*(length - 1) + k] < 0)) {
				cumHeight[clusterNum - 1] = rans[4*(length - 1) + k];
				
				// alternative cluster numbering
				rans[6*(length - 1) + k] = rans[0*(length - 1) + k];
				rans[7*(length - 1) + k] = rans[1*(length - 1) + k];
			} else {
				cumHeight[clusterNum - 1] = rans[3*(length - 1) + k];
				rans[4*(length - 1) + k] -= cumHeight[(int)rans[0*(length - 1) + k] - 1]; // col
				rans[3*(length - 1) + k] -= cumHeight[(int)rans[1*(length - 1) + k] - 1]; // row
				
				// alternative cluster numbering
				rans[6*(length - 1) + k] = clusterNums[(int)rans[0*(length - 1) + k] - 1];
				rans[7*(length - 1) + k] = clusterNums[(int)rans[1*(length - 1) + k] - 1];
			}
			clusterNums[clusterNum - 1] = k + 1;
			rans[5*(length - 1) + k] = cumHeight[clusterNum - 1];
		} else if (rans[0*(length - 1) + k] > 0) {
			// row is a cluster from before
			rans[2*(length - 1) + k] = rans[0*(length - 1) + k]; // merge with previous cluster
			// calculate both branch lengths
			rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2; // col
			rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2 - cumHeight[(int)rans[0*(length - 1) + k] - 1]; // row
			
			// add longest path to cumulative height
			cumHeight[(int)rans[0*(length - 1) + k] - 1] += rans[3*(length - 1) + k];
			if (rans[4*(length - 1) + k] > cumHeight[(int)rans[0*(length - 1) + k] - 1]) {
				cumHeight[(int)rans[0*(length - 1) + k] - 1] = rans[4*(length - 1) + k];
			}
			rans[5*(length - 1) + k] = cumHeight[(int)rans[0*(length - 1) + k] - 1];
			
			// alternative cluster numbering
			rans[6*(length - 1) + k] = clusterNums[(int)rans[0*(length - 1) + k] - 1];
			rans[7*(length - 1) + k] = rans[1*(length - 1) + k];
			clusterNums[(int)rans[0*(length - 1) + k] - 1] = k + 1;
		} else if (rans[1*(length - 1) + k] > 0) {
			rans[2*(length - 1) + k] = rans[1*(length - 1) + k]; // merge with previous cluster
			// calculate both branch lengths
			rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2 - cumHeight[(int)rans[1*(length - 1) + k] - 1]; // col
			rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2; // row
			
			// add longest path to cumulative height
			cumHeight[(int)rans[1*(length - 1) + k] - 1] += rans[4*(length - 1) + k];
			if (rans[3*(length - 1) + k] > cumHeight[(int)rans[1*(length - 1) + k] - 1]) {
				cumHeight[(int)rans[1*(length - 1) + k] - 1] = rans[3*(length - 1) + k];
			}
			rans[5*(length - 1) + k] = cumHeight[(int)rans[1*(length - 1) + k] - 1];
			
			// alternative cluster numbering
			rans[7*(length - 1) + k] = clusterNums[(int)rans[1*(length - 1) + k] - 1];
			rans[6*(length - 1) + k] = rans[0*(length - 1) + k];
			clusterNums[(int)rans[1*(length - 1) + k] - 1] = k + 1;
		}
		
		// calculate distances to the new node/cluster
		if (met==4) { // complete
			index = 0;
			for (i = -1; i < (size - 1); i++) {
				if (!(i==minRow) && !(i==minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dist1 = *(dMatrix + minRow*length + (i + 1));
					} else {
						dist1 = *(dMatrix + i*length + (minRow + 1));
					}
					if (i >= minCol) {
						dist2 = *(dMatrix + i*length + minCol);
					} else {
						dist2 = *(dMatrix + (minCol - 1)*length + (i + 1));
					}
					if (dist1 > dist2) { // pick max distance
						dTemp[index] = dist1;
					} else {
						dTemp[index] = dist2;
					}
					index++;
				}
			}
		} else if (met==5) { // single
			index = 0;
			for (i = -1; i < (size - 1); i++) {
				if (!(i==minRow) && !(i==minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dist1 = *(dMatrix + minRow*length + (i + 1));
					} else {
						dist1 = *(dMatrix + i*length + (minRow + 1));
					}
					if (i >= minCol) {
						dist2 = *(dMatrix + i*length + minCol);
					} else {
						dist2 = *(dMatrix + (minCol - 1)*length + (i + 1));
					}
					if (dist1 < dist2) { // pick min distance
						dTemp[index] = dist1;
					} else {
						dTemp[index] = dist2;
					}
					index++;
				}
			}
		} else { // UPGMA/average
			index = 0;
			for (i = -1; i < (size - 1); i++) {
				if (!(i==minRow) && !(i==minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dTemp[index] += *(dMatrix + minRow*length + (i + 1));
					} else {
						dTemp[index] += *(dMatrix + i*length + (minRow + 1));
					}
					if (i >= minCol) {
						dTemp[index] += *(dMatrix + i*length + minCol);
					} else {
						dTemp[index] += *(dMatrix + (minCol - 1)*length + (i + 1));
					}
					dTemp[index] /= 2; // average distance
					index++;
				}
			}
		}
		
		// make the new distance matrix
		// (note: minRow is always greater than minCol)
		#pragma omp parallel for private(i,j) schedule(guided)
		for (j = 0; j < (size - 1); j++) { // for each column
			// move each row up by one after minRow
			int start = minRow + 1;
			if (j >= start) {
				start = j;
			}
			for (i = start; i < (size - 1); i++) { // for each row
				dMatrix[(i - 1)*length + j] = dMatrix[i*length + j];
			}
		}
		
		#pragma omp parallel for private(i,j) schedule(guided)
		for (i = minRow; i < (size - 1); i++) { // for each row
			// move each column left after minRow
			for (j = minRow + 1; j <= i; j++) { // for each column
				dMatrix[i*length + j] = dMatrix[i*length + j + 1];
			}
		}
		
		// move the cluster numbers
		for (i = minRow + 1; i < (size - 1); i++) {
			dMatrix[(i - 1)*length + (length - 1)] = dMatrix[i*length + (length - 1)];
		}
		for (j = minRow + 2; j < (size - 1); j++) {
			dMatrix[(length - 1)*length + (j - 1)] = dMatrix[(length - 1)*length + j];
		}
		
		// give the cluster its new number
		*(dMatrix + (length - 1)*length + minCol) = rans[2*(length - 1) + k];
		if ((minCol - 1) >= 0)
			*(dMatrix + (minCol - 1)*length + (length - 1)) = rans[2*(length - 1) + k];
		
		// decrement size of the matricies; length remains constant
		size--;
		
		// put new distances into the new cluster at minCol
		index = 0;
		for (j = 0; j < minCol; j++) {
			*(dMatrix + (minCol - 1)*length + j) = dTemp[index];
			index++;
		}
		for (i = minCol; i < (size - 1); i++) {
			*(dMatrix + i*length + minCol) = dTemp[index];
			index++;
		}
		
		if (v) {
			// print the percent completed so far
			soFar += pow(length - 1 - k, 2);
			*rPercentComplete = floor(100*soFar/total);
			
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			void R_CheckUserInterrupt(void);
		}
	}
	
	// bin sequences
	// assign cluster numbers to nodes within cutoff percent different
	int clusterNumber = 1;
	for (i = 0; i < length - 1; i++) {
		if (rans[5*(length - 1) + i] > *cut/2 &&
			rans[8*(length - 1) + i] == 0 && // first fork is unassigned
			rans[6*(length - 1) + i] < 0 && // first fork is a leaf
			rans[9*(length - 1) + i] == 0 && // second fork is unassigned
			rans[7*(length - 1) + i] < 0) { // second fork is a leaf
			rans[8*(length - 1) + i] = clusterNumber;
			clusterNumber++;
			rans[9*(length - 1) + i] = clusterNumber;
			clusterNumber++;
		} else {
			if (rans[8*(length - 1) + i] == 0 && // first fork is unassigned
				rans[6*(length - 1) + i] < 0) { // first fork is a leaf
				binUPGMA(rans, i, clusterNumber, *cut/2, length);
				clusterNumber++;
			}
			if (rans[9*(length - 1) + i] == 0 && // second fork is unassigned
				rans[7*(length - 1) + i] < 0) { // second fork is a leaf
				binUPGMA(rans, i, clusterNumber, *cut/2, length);
				clusterNumber++;
			}
		}
	}
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}

SEXP reclusterUPGMA(SEXP ans, SEXP cutoff)
{
	int i;
	double *cut, *rans;
	cut = REAL(cutoff);
	rans = REAL(ans);
	const int length = length(ans)/10 + 1; // number of rows
	
	// zero out previous clusters
	for (i = 0; i < length - 1; i++) {
		rans[8*(length - 1) + i] = 0;
		rans[9*(length - 1) + i] = 0;
	}
	
	// rebin sequences
	// assign cluster numbers to nodes within cutoff percent different
	int clusterNumber = 1;
	for (i = 0; i < length - 1; i++) {
		if (rans[5*(length - 1) + i] > *cut/2 &&
			rans[8*(length - 1) + i] == 0 && // first fork is unassigned
			rans[6*(length - 1) + i] < 0 && // first fork is a leaf
			rans[9*(length - 1) + i] == 0 && // second fork is unassigned
			rans[7*(length - 1) + i] < 0) { // second fork is a leaf
			rans[8*(length - 1) + i] = clusterNumber;
			clusterNumber++;
			rans[9*(length - 1) + i] = clusterNumber;
			clusterNumber++;
		} else {
			if (rans[8*(length - 1) + i] == 0 && // first fork is unassigned
				rans[6*(length - 1) + i] < 0) { // first fork is a leaf
				binUPGMA(rans, i, clusterNumber, *cut/2, length);
				clusterNumber++;
			}
			if (rans[9*(length - 1) + i] == 0 && // second fork is unassigned
				rans[7*(length - 1) + i] < 0) { // second fork is a leaf
				binUPGMA(rans, i, clusterNumber, *cut/2, length);
				clusterNumber++;
			}
		}
	}
	
	return ans;
}
