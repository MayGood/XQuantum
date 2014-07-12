/*
 * lin.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "lin.h"
#include "lapacke.h"


int gen_sym_eigen(
	int nDim,                           // array dimension
	const float *A, const float *S,   // input A and S
	float *e, float *C){              // output value and vector

	int i, j, k;   // looping

	//
	// LAPACK dspgvd_ variables
	//
	long int ITYPE = 1;  // A*x = (lambda)*B*x
	char     JOBZ = 'V'; // compute eigenvalues and vectors
	char     UPLO = 'L'; // upper triangle of A and B are stored
	long int N;        // dimension of matrices A and B
	float  *AP;       // upper triangle matrix A
	float  *BP;       // upper triangle matrix B
	float  *W;        // output eigen value
	float  *Z;        // output eigen vectors
	long int LDZ;      // dimension of matrix
	float  *WORK;     // working array size LWORK
	long int LWORK;    // should be >= 1 + 6N + 2N^2
	long int *IWORK;   // array dimension LIWORK
	long int LIWORK;   // shouble be >= 3 + 5N
	long int INFO;     // exit status

	// sanity check
	if (nDim < 1){
		printf("gen_sym_eigen: Error - Invalid dimension range\n");
		exit(-1);
	}

	// allocate memory
	N = nDim;
	LDZ = nDim;
	AP = calloc(N*(N + 1) / 2, sizeof(float));
	BP = calloc(N*(N + 1) / 2, sizeof(float));
	W = calloc(N, sizeof(float));
	Z = calloc(LDZ*N, sizeof(float));

	// computation of LWORK is from DSPGVD code
	i = (int)(log(N) / log(2.0));
	if (pow(2.0, i) < (float)N) i++;
	if (pow(2.0, i) < (float)N) i++;
	LWORK = 1 + 5 * N + 2 * N*N + 2 * N*i;

	WORK = calloc(LWORK, sizeof(float));
	LIWORK = 3 + 5 * N;
	IWORK = calloc(LIWORK, sizeof(long int));
	INFO = 0;

	// read lower triangle array A
	// Important!!!!
	// According to LAPACK document, it should be
	// upper. But from emprical experiment, you really
	// need to enter the matrix as lower triangle
	// - Teepanis Chachiyo Dec 23, 2007
	k = 0;

	for (i = 0; i < N; i++) {
		for (j = i; j < N; j++){
			AP[k] = A[i*N + j];
			k++;
		}
	}


	// read lower triangle array B
	k = 0;

	for (i = 0; i < N; i++) {
		for (j = i; j < N; j++){
			BP[k] = S[i*N + j];
			k++;
		}
	}

	sspgvd_(&ITYPE, &JOBZ, &UPLO, &N, AP, BP, W, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK,	&INFO);

	if (INFO != 0){
		printf("gen_sym_eigen: Error - LAPACK subroutine returns error %ld\n"
			, INFO);
		exit(EXIT_FAILURE);
	}
	// print eigen values
	//for(i=0; i < N; i++){
	//	fprintf(outFile, "%20.8lE", W[i]);
	//}
	//fprintf(outFile, "\n");
	for (i = 0; i < N; i++)
		e[i] = W[i];

	// print eigen vector per row
	//for(i=0; i < N; i++){
	//	for(j=0; j < N; j++){
	//		fprintf(outFile, "%20.8lE", Z[i*N+j]);
	//	}
	//	fprintf(outFile, "\n");
	//}
	for (i = 0; i < N; i++)
	for (j = 0; j < N; j++)
		C[i*N + j] = Z[i*N + j];

	// clean up memory
	free(AP);
	free(BP);
	free(W);
	free(Z);
	free(WORK);
	free(IWORK);

	// return TRUE on success
	return 1;
}

int linear_solver(int nRow,        // number of row
	int nCol,        // number column
	const float *A, // input matrix A
	float *B){      // input matrix B, output matrix X

	int i, j;

	//
	// LAPACK dgesv_ variables
	//
	float *LA;    // LAPACK routine matrix
	float *LB;    // LAPACK routine matrix
	long int N;    // number of rows
	long int NRHS; // number of columns
	long int LDA;  // number of rows
	long int *IPIV;// internal usage
	long int LDB;  // number of rows
	long int INFO; // return code

	// validate range
	if (nRow < 1 || nCol < 1){
		printf("linear_solver: Error - invalid dimension range\n");
		exit(-1);
	}

	// setup LAPACK variables
	N = nRow;
	NRHS = nCol;
	LDA = nRow;
	LDB = nRow;
	INFO = 0;

	// allocate memory
	LA = calloc(nRow*nRow, sizeof(float));
	LB = calloc(nCol*nRow, sizeof(float));
	IPIV = calloc(nRow, sizeof(long int));
	if (IPIV == NULL || LA == NULL || LB == NULL){
		printf("linear_solver: Error - cannnot allocate memeory\n");
		exit(-1);
	}

	// read matrix values using transpose due to LAPACK's FORTRAN notation
	for (i = 0; i < nRow; i++)
	for (j = 0; j < nRow; j++)
		LA[i*nRow + j] = A[j*nRow + i];

	for (i = 0; i < nCol; i++)
	for (j = 0; j < nRow; j++)
		LB[i*nRow + j] = B[j*nCol + i];
	sgesv_(&N, &NRHS, LA, &LDA, IPIV, LB, &LDB, &INFO);
	// check for error
	if (INFO != 0){
		printf("linear_solver: Error - LAPACK returns error code %ld\n", INFO);
		exit(-1);
	}

	// output
	for (i = 0; i < nRow; i++)
	for (j = 0; j < nCol; j++)
		B[i*nCol + j] = LB[j*nRow + i];

	// clean memory
	free(IPIV);
	free(LA);
	free(LB);

	// return TRUE on success
	return 1;
}
