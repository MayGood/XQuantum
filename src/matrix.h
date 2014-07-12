/*
 * matrix.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "basis.h"
#include "mol.h"
#include "int.h"
#include "option.h"

// Schwarz inequality cut-off for 2e integral
#define SCHWARZ_CUTOFF 1.0E-12

float GTO_moment(int i,                            // ith basis
	int j,                            // jth basis
	const struct GTOBasis_t *gto,     // basis database
	int mx, int my, int mz,           // moment order
	float xc, float yc, float zc); // moment center

float GTO_overlap(int i,                        // ith basis
	int j,                        // jth basis
	const struct GTOBasis_t *gto);// basis database

float GTO_kinetic(int i,                        // ith basis
	int j,                        // jth basis
	const struct GTOBasis_t *gto);// basis database

float GTO_nuclei(int i,                        // ith basis
	int j,                        // jth basis
	const struct GTOBasis_t *gto, // basis database
	const struct Molecule_t *mol);// molecule database

float * create_Schwarz(
	int nBasis,                     // number of basis functions
	const struct GTOBasis_t *gto);  // basis set info

void GTO_JK_Matrix_Quartet(
	int nBasis,                    // number of basis functions
	const float *PA,              // density matrix for spin up
	const float *PB,              // density matrix for spin down
	const struct GTOBasis_t *gto,  // basis set info
	const float *schwarz_basis,   // pointer to schwarz matrix
	float fixedCutoff,            // cutoff to ignore
	float *GA,                    // return G for spin up
	float *GB,                    // return G for spin down
	struct option_t *opt);         // global option

void GTO_JK_Matrix_Quartet_Parallel(
	int childID,                   // child id number
	int nBasis,                    // number of basis functions
	const float *PA,              // density matrix for spin up
	const float *PB,              // density matrix for spin down
	const struct GTOBasis_t *gto,  // basis set info
	const float *schwarz_basis,   // pointer to schwarz matrix
	float fixedCutoff,            // cutoff to ignore
	float *GA,                    // return G for spin up
	float *GB,                    // return G for spin down
	struct option_t *opt);         // global option

void GTO_JK_Matrix(
	int nBasis,                  // number of basis functions
	float *P,                   // density matrix
	struct GTOBasis_t *gto,      // basis set info
	float *Schwarz,             // pointer to schwarz matrix
	float cutoff,               // cutoff to ignore
	float *G);                  // return matrix

void GTO_JK_Matrix1(
	int nBasis,                  // number of basis functions
	float *P,                   // density matrix
	struct GTOBasis_t *gto,      // basis set info
	float *Schwarz,             // pointer to schwarz matrix
	float cutoff,               // cutoff to ignore
	float *G);                  // return matrix

int isSameShell(
	int nBasis,                    // number of basis functions
	int i, int j,                  // which basis to test
	const struct GTOBasis_t *gto); // basis set info

#endif /* MATRIX_H_ */
