/*
 * uhf.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef UHF_H_
#define UHF_H_

#include <math.h>
#include "mol.h"
#include "basis.h"
#include "option.h"

float uhf(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nEA,                 // total number of spin up electrons
	int nEB,                 // total number of spin down electrons
	float *CA,              // returned molecular alpha spin orbital
	float *CB,              // returned molecular beta spin orbital
	float *eA,              // returned eigen values
	float *eB,              // returned eigen values
	struct option_t *opt);   // global option

#endif /* UHF_H_ */
