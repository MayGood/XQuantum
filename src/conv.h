/*
 * conv.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef CONV_H_
#define CONV_H_

float conv_damping(int nBasis,     // number of basis
	float alpha,   // drag coefficient between 0 and 1
	float *PA,     // alpha density matrix
	float *PB);    // beta density matrix

float conv_diis2(int nBasis,      // number of basis
	float alpha,    // drag coefficient between 0 and 1
	float *PA,      // alpha density matrix
	float *PB);     // beta density matrix

float conv_diis3(int nBasis,      // number of basis
	float alpha,    // drag coefficient between 0 and 1
	float *PA,      // alpha density matrix
	float *PB);     // beta density matrix

float conv_diis4(int nBasis,      // number of basis
	float alpha,    // drag coefficient between 0 and 1
	float *PA,      // output alpha density matrix
	float *PB);     // output beta  density matrix

#endif /* CONV_H_ */
