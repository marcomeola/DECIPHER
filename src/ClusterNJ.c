/****************************************************************************
 *                         Cluster Neighbor Joining                         *
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

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

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

// DECIPHER	header file
#include "DECIPHER.h"

void FollowBranch(double *rans, int i, double *branchLength, int length) {
	// add the longest branch length
	if (rans[8*(length - 1) + i] == 0) { // cluster number unassigned
		double alternative;
		if (rans[6*(length - 1) + i] < 0 && // first fork is a leaf and
			rans[7*(length - 1) + i] < 0) { // second fork is a leaf
			// add the longest one to branch length
			if (rans[3*(length - 1) + i] < rans[4*(length - 1) + i] && // second leaf is longest and
				rans[9*(length - 1) + i]==0) { // second leaf is not assigned to a cluster
				// add second leaf to branch length
				*branchLength += rans[4*(length - 1) + i];
			} else {
				// add first leaf to branch length
				*branchLength += rans[3*(length - 1) + i];
			}
		} else if (rans[6*(length - 1) + i] > 0) { // first fork is a branch
			alternative = *branchLength + rans[4*(length - 1) + i];
			*branchLength += rans[3*(length - 1) + i];
			FollowBranch(rans, (int)(rans[6*(length - 1) + i] - 1), branchLength, length);
			if (*branchLength < alternative) {
				*branchLength = alternative;
			}
		} else if (rans[7*(length - 1) + i] > 0) { // second fork is a branch
			alternative = *branchLength + rans[3*(length - 1) + i];
			*branchLength += rans[4*(length - 1) + i];
			FollowBranch(rans, (int)(rans[7*(length - 1) + i] - 1), branchLength, length);
			if (*branchLength < alternative) {
				*branchLength = alternative;
			}
		}
	}

}

void assignNumber(double *rans, int i, int clusterNumber, double maxHeight, double floorHeight, int length) {
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
		double branchLength;
		for (int j = i + 1; j < length - 1; j++) {
			if (rans[6*(length - 1) + j] == (i + 1) || rans[7*(length - 1) + j] == (i + 1)) { // node is merged again
				branchLength = 0;
				FollowBranch(rans, j, &branchLength, length);
				if ((rans[5*(length - 1) + j] + branchLength) <= maxHeight) { // if the node is within reach
					// then assign the same clusterNumber to this node
					assignNumber(rans, j, clusterNumber, maxHeight, floorHeight, length);
					break;
				}
			}
		}
	}
	// follow the branches down the tree
	if (rans[6*(length - 1) + i] > 0) { // if first fork is a branch
		if (rans[5*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] >= floorHeight) { // condition to stop following branches
			// then recursively number the branch
			if (((rans[5*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] -
				  rans[3*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] >= floorHeight) || // first fork is within reach or
				 rans[8*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] != 0) && // the fork is already assigned and
				((rans[5*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] -
				  rans[4*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] >= floorHeight) || // second fork is within reach or
				 rans[9*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] != 0)) { // the fork is already assigned
				assignNumber(rans, (int)(rans[6*(length - 1) + i] - 1), clusterNumber, maxHeight, floorHeight, length);
			}
		}
	}
	if (rans[7*(length - 1) + i] > 0) { // if second fork is a branch
		if (rans[5*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] >= floorHeight) { // condition to stop following branches
			// then recursively number the branch
			if (((rans[5*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] -
				  rans[3*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] >= floorHeight) || // first fork is within reach or
				 rans[8*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] != 0) && // the fork is already assigned and
				((rans[5*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] -
				  rans[4*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] >= floorHeight) || // second fork is within reach or
				 rans[9*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] != 0)) { // the fork is already assigned
					assignNumber(rans, (int)(rans[7*(length - 1) + i] - 1), clusterNumber, maxHeight, floorHeight, length);
			}
		}
	}
}

void Offset(int i, double *rans, double *offset, int length) {
	for (int j = i + 1; j < length - 1; j++) {
		if (rans[6*(length - 1) + j] == (i + 1)) {
			*offset = *offset + rans[5*(length - 1) + j] - rans[5*(length - 1) + i] - rans[3*(length - 1) + j];
			Offset(j, rans, offset, length);
			break;
		}
		if (rans[7*(length - 1) + j] == (i + 1)) {
			*offset = *offset + rans[5*(length - 1) + j] - rans[5*(length - 1) + i] - rans[4*(length - 1) + j];
			Offset(j, rans, offset, length);
			break;
		}
	}
}

//ans_start <- .Call("cluster", myDistMatrix, verbose, pBar, PACKAGE="DECIPHER")
SEXP clusterNJ(SEXP x, SEXP cutoff, SEXP verbose, SEXP pBar, SEXP nThreads)
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
	// nDiv is the net divergence
	// cluster number is stored in the outer col, row of dMatrix
	// single leaves are negative, clusters are positive
	// Cluster numbering notes:
	// starting at lowest leaf, climb the tree until over maxHeight
	// maxHeight is merge height plus cutoff minus longest branch
	// do not go down the tree past floorHeight
	
	int i, j, k, clusterNum, size, minRow, minCol, index, minR, minC;
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
	double nDiv[size];
	double *dTemp = (double *) R_alloc(size - 1, sizeof(double));
	double *cumHeight = (double *) R_alloc(size - 1, sizeof(double));
	int *clusterNums = (int *) R_alloc(size - 1, sizeof(int));
	cut = REAL(cutoff);
	v = asLogical(verbose);
	int nthreads = asInteger(nThreads);
	
	if (v) { // initialize progress variables
		soFar = 0;
		before = 0;
		total = length*(length - 1);
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
		// zero out the net divergence
		for (i = 0; i < size; i++) {
			nDiv[i] = 0;
		}
		// calculate the net divergence
		for (i = 0; i < (size - 1); i++) {
			for (j = 0; j <= i; j++) {
				nDiv[j] += *(dMatrix + i*length + j); // col sums
				nDiv[i + 1] += *(dMatrix + i*length + j); // row sums
			}
		}
		
		// calculate the Q matrix & find the smallest element in the Q matrix
		minRow = 0;
		minCol = 0;
		minHeight = 1e50;
		#pragma omp parallel for private(i,j,minR,minC,minH) schedule(guided) num_threads(nthreads)
		for (i = (size - 2); i >= 0; i--) {
			minH = minHeight;
			for (j = 0; j <= i; j++) {
				if (*(dMatrix + i*length + j) - (nDiv[i + 1] + nDiv[j])/(size - 2) < minH) {
					minH = *(dMatrix + i*length + j) - (nDiv[i + 1] + nDiv[j])/(size - 2);
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
			if ((size - 2)==0) { // case of (0/0 == NaN)
				rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2; // col
				rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2; // row
			} else {
				rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2 + (nDiv[minCol] - nDiv[minRow + 1])/(2*(size-2)); // col
				rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol) - rans[4*(length - 1) + k]; // row
			}
			
			// zero negative branch lengths
			if (rans[3*(length - 1) + k] < 0 && rans[4*(length - 1) + k] < 0) {
				rans[3*(length - 1) + k] = 0;
				rans[4*(length - 1) + k] = 0;
			} else if (rans[4*(length - 1) + k] < 0) {
				rans[3*(length - 1) + k] = -1*rans[4*(length - 1) + k]; // add difference to other branch
				rans[4*(length - 1) + k] = 0;
			} else if (rans[3*(length - 1) + k] < 0) {
				rans[4*(length - 1) + k] = -1*rans[3*(length - 1) + k];
				rans[3*(length - 1) + k] = 0;
			}
			
			// add longest path to cumulative height
			if ((rans[0*(length - 1) + k] < 0) && (rans[1*(length - 1) + k] < 0)) {
				if (rans[4*(length - 1) + k] > rans[3*(length - 1) + k]) {
					cumHeight[clusterNum - 1] = rans[4*(length - 1) + k];
				} else {
					cumHeight[clusterNum - 1] = rans[3*(length - 1) + k];
				}
				
				// alternative cluster numbering
				rans[6*(length - 1) + k] = rans[0*(length - 1) + k];
				rans[7*(length - 1) + k] = rans[1*(length - 1) + k];
			} else {
				if ((cumHeight[(int)rans[0*(length - 1) + k] - 1] + rans[3*(length - 1) + k]) >
					(cumHeight[(int)rans[1*(length - 1) + k] - 1] + rans[4*(length - 1) + k])) {
					cumHeight[clusterNum - 1] = cumHeight[(int)rans[0*(length - 1) + k] - 1] + rans[3*(length - 1) + k];
				} else {
					cumHeight[clusterNum - 1] = cumHeight[(int)rans[1*(length - 1) + k] - 1] + rans[4*(length - 1) + k];
				}
				
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
			if ((size - 2)==0) { // case of (0/0 == NaN)
				rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol); // col
				rans[3*(length - 1) + k] = 0; // row
			} else {
				rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2 + (nDiv[minCol] - nDiv[minRow + 1])/(2*(size-2)); // col
				rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol) - rans[4*(length - 1) + k]; // row
			}
			
			// zero negative branch lengths
			if (rans[3*(length - 1) + k] < 0 && rans[4*(length - 1) + k] < 0) {
				rans[3*(length - 1) + k] = 0;
				rans[4*(length - 1) + k] = 0;
			} else if (rans[4*(length - 1) + k] < 0) {
				rans[3*(length - 1) + k] = -1*rans[4*(length - 1) + k]; // add difference to other branch
				rans[4*(length - 1) + k] = 0;
			} else if (rans[3*(length - 1) + k] < 0) {
				rans[4*(length - 1) + k] = -1*rans[3*(length - 1) + k];
				rans[3*(length - 1) + k] = 0;
			}
			
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
			if ((size - 2)==0) { // case of (0/0 == NaN)
				// unclear what to do in this case - splitting branch (edge) lengths as compromise
				rans[4*(length - 1) + k] = 0; // col
				rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol); // row
			} else {
				rans[4*(length - 1) + k] = *(dMatrix + minRow*length + minCol)/2 + (nDiv[minCol] - nDiv[minRow + 1])/(2*(size-2)); // col
				rans[3*(length - 1) + k] = *(dMatrix + minRow*length + minCol) - rans[4*(length - 1) + k]; // row
			}
			
			// zero negative branch lengths
			if (rans[3*(length - 1) + k] < 0 && rans[4*(length - 1) + k] < 0) {
				rans[3*(length - 1) + k] = 0;
				rans[4*(length - 1) + k] = 0;
			} else if (rans[4*(length - 1) + k] < 0) {
				rans[3*(length - 1) + k] = -1*rans[4*(length - 1) + k]; // add difference to other branch
				rans[4*(length - 1) + k] = 0;
			} else if (rans[3*(length - 1) + k] < 0) {
				rans[4*(length - 1) + k] = -1*rans[3*(length - 1) + k];
				rans[3*(length - 1) + k] = 0;
			}
			
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
				dTemp[index] -= *(dMatrix + minRow*length + minCol);
				dTemp[index] /= 2;
				index++;
			}
		}
		
		// make the new distance matrix
		// (note: minRow is always greater than minCol)
		#pragma omp parallel for private(i,j) schedule(guided) num_threads(nthreads)
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
		
		#pragma omp parallel for private(i,j) schedule(guided) num_threads(nthreads)
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
			soFar = (2*length - 2 - k)*(k + 1);
			*rPercentComplete = floor(100*soFar/total);
			
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	
	// insure that nodes are at the correct heights
	double offset;
	for (i = 0; i < length - 1; i++) {
		offset = 0;
		Offset(i, rans, &offset, length);
		rans[5*(length - 1) + i] = rans[5*(length - 1) + i] + offset;
	}
	// bin sequences
	// assign cluster numbers to nodes within cutoff percent different
	int clusterNumber = 1;
	double maxHeight, lowHeight, tempHeight;
	double longLeaf, longestLeaf, floorHeight;
	double tempSwap; // for swapping order of leaves
	int lowestBranch;
	for (i = 0; i < length; i++) {
		// find the lowest unassigned leaf
		lowestBranch = -1;
		lowHeight = 1e50;
		for (j = 0; j < length - 1; j++) {
			if (rans[8*(length - 1) + j] == 0 && // cluster number unassigned
				(rans[6*(length - 1) + j] < 0 || // first fork is a leaf or
				rans[7*(length - 1) + j] < 0)) { // second fork is a leaf
				if (rans[6*(length - 1) + j] < 0 && // first fork is a leaf and
					rans[7*(length - 1) + j] < 0) { // second fork is a leaf
					if (rans[3*(length - 1) + j] < rans[4*(length - 1) + j] && // second leaf is longest and
						rans[9*(length - 1) + j] == 0) { // second cluster number is unassigned
						tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
						longLeaf = rans[4*(length - 1) + j];
					} else {
						tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
						longLeaf = rans[3*(length - 1) + j];
					}
				} else if (rans[6*(length - 1) + j] < 0) { // first fork is a leaf
					tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
					longLeaf = rans[3*(length - 1) + j];
				} else if (rans[7*(length - 1) + j] < 0) { // second fork is a leaf
					tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
					longLeaf = rans[4*(length - 1) + j];
				}
				if (tempHeight < lowHeight) { // found a lower leaf
					lowestBranch = j;
					lowHeight = tempHeight;
					longestLeaf = longLeaf;
				}
			}
		}
		// assign a number to the lowest leaf
		if (lowestBranch >= 0) {
			// handle the special case where the two leaves are far apart
			if (rans[6*(length - 1) + lowestBranch] < 0 && // first fork is a leaf and
				rans[7*(length - 1) + lowestBranch] < 0 && // second fork is a leaf and
				rans[9*(length - 1) + lowestBranch] == 0 && // second cluster number is unassigned and
				((rans[3*(length - 1) + lowestBranch] + rans[4*(length - 1) + lowestBranch]) > *cut)) { // leaves farther apart then cutoff
				// then assign the longest leaf its own cluster number
				rans[9*(length - 1) + lowestBranch] = clusterNumber;
				if (rans[3*(length - 1) + lowestBranch] > rans[4*(length - 1) + lowestBranch]) { // first leaf is longest
					// then swap leaves
					tempSwap = rans[4*(length - 1) + lowestBranch];
					rans[4*(length - 1) + lowestBranch] = rans[3*(length - 1) + lowestBranch];
					rans[3*(length - 1) + lowestBranch] = tempSwap;
					tempSwap = rans[1*(length - 1) + lowestBranch];
					rans[1*(length - 1) + lowestBranch] = rans[0*(length - 1) + lowestBranch];
					rans[0*(length - 1) + lowestBranch] = tempSwap;
					tempSwap = rans[7*(length - 1) + lowestBranch];
					rans[7*(length - 1) + lowestBranch] = rans[6*(length - 1) + lowestBranch];
					rans[6*(length - 1) + lowestBranch] = tempSwap;
				}
			} else {
				// assign clusters going down the tree
				floorHeight = lowHeight - *cut + 2*longestLeaf;
				maxHeight = lowHeight + *cut;
				assignNumber(rans, lowestBranch, clusterNumber, maxHeight, floorHeight, length);
			}
			clusterNumber++;
		} else { // no leaves left to assign
			break;
		}

	}
		
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}

SEXP reclusterNJ(SEXP x, SEXP cutoff)
{
	// bin sequences
	// assign cluster numbers to nodes within cutoff percent different
	int clusterNumber = 1;
	double maxHeight, lowHeight, tempHeight;
	double longLeaf, longestLeaf, floorHeight;
	double tempSwap; // for swapping order of leaves
	int lowestBranch;
	
	int i, j;
	double *cut, *rans;
	cut = REAL(cutoff);
	SEXP ans;
	PROTECT(ans = duplicate(x));
	rans = REAL(ans);
	const int length = length(ans)/10 + 1; // number of rows
	
	// zero out previous clusters
	for (i = 0; i < length - 1; i++) {
		rans[8*(length - 1) + i] = 0;
		rans[9*(length - 1) + i] = 0;
	}
	
	for (i = 0; i < length; i++) {
		// find the lowest unassigned leaf
		lowestBranch = -1;
		lowHeight = 1e50;
		for (j = 0; j < length - 1; j++) {
			if (rans[8*(length - 1) + j] == 0 && // cluster number unassigned
				(rans[6*(length - 1) + j] < 0 || // first fork is a leaf or
				 rans[7*(length - 1) + j] < 0)) { // second fork is a leaf
					if (rans[6*(length - 1) + j] < 0 && // first fork is a leaf and
						rans[7*(length - 1) + j] < 0) { // second fork is a leaf
						if (rans[3*(length - 1) + j] < rans[4*(length - 1) + j] && // second leaf is longest and
							rans[9*(length - 1) + j] == 0) { // second cluster number is unassigned
							tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
							longLeaf = rans[4*(length - 1) + j];
						} else {
							tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
							longLeaf = rans[3*(length - 1) + j];
						}
					} else if (rans[6*(length - 1) + j] < 0) { // first fork is a leaf
						tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
						longLeaf = rans[3*(length - 1) + j];
					} else if (rans[7*(length - 1) + j] < 0) { // second fork is a leaf
						tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
						longLeaf = rans[4*(length - 1) + j];
					}
					if (tempHeight < lowHeight) { // found a lower leaf
						lowestBranch = j;
						lowHeight = tempHeight;
						longestLeaf = longLeaf;
					}
				}
		}
		// assign a number to the lowest leaf
		if (lowestBranch >= 0) {
			// handle the special case where the two leaves are far apart
			if (rans[6*(length - 1) + lowestBranch] < 0 && // first fork is a leaf and
				rans[7*(length - 1) + lowestBranch] < 0 && // second fork is a leaf and
				rans[9*(length - 1) + lowestBranch] == 0 && // second cluster number is unassigned and
				((rans[3*(length - 1) + lowestBranch] + rans[4*(length - 1) + lowestBranch]) > *cut)) { // leaves farther apart then cutoff
				// then assign the longest leaf its own cluster number
				rans[9*(length - 1) + lowestBranch] = clusterNumber;
				if (rans[3*(length - 1) + lowestBranch] > rans[4*(length - 1) + lowestBranch]) { // first leaf is longest
					// then swap leaves
					tempSwap = rans[4*(length - 1) + lowestBranch];
					rans[4*(length - 1) + lowestBranch] = rans[3*(length - 1) + lowestBranch];
					rans[3*(length - 1) + lowestBranch] = tempSwap;
					tempSwap = rans[1*(length - 1) + lowestBranch];
					rans[1*(length - 1) + lowestBranch] = rans[0*(length - 1) + lowestBranch];
					rans[0*(length - 1) + lowestBranch] = tempSwap;
					tempSwap = rans[7*(length - 1) + lowestBranch];
					rans[7*(length - 1) + lowestBranch] = rans[6*(length - 1) + lowestBranch];
					rans[6*(length - 1) + lowestBranch] = tempSwap;
				}
			} else {
				// assign clusters going down the tree
				floorHeight = lowHeight - *cut + 2*longestLeaf;
				maxHeight = lowHeight + *cut;
				assignNumber(rans, lowestBranch, clusterNumber, maxHeight, floorHeight, length);
			}
			clusterNumber++;
		} else { // no leaves left to assign
			break;
		}
		
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP adjustHeights(SEXP x)
{
	// insure that nodes are at the correct heights
	double offset, *rans;
	int length = length(x)/10 + 1;
	rans = REAL(x);
	
	for (int i = 0; i < length - 1; i++) {
		offset = 0;
		Offset(i, rans, &offset, length);
		rans[5*(length - 1) + i] = rans[5*(length - 1) + i] + offset;
	}
	
	return x;
}
