/*
 * lin.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef LIN_H_
#define LIN_H_

int gen_sym_eigen(
	int nDim,                           // array dimension
	const float *A, const float *S,   // input A and S
	float *e, float *C);              // output value and vector

int linear_solver(int nRow,         // number of row
	int nCol,         // number column
	const float *A,  // input matrix A
	float *B);       // input matrix B, output matrix X

#endif /* LIN_H_ */
