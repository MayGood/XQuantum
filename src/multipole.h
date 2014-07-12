/*
 * multipole.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef MULTIPOLE_H_
#define MULTIPOLE_H_

#include "basis.h"

#define MULTIPOLE_CHARGE_CUTOFF 1.0E-18
#define MULTIPOLE_RADII2_CUTOFF 1.0E-00

// multipole expansion of a product between two basis functions
struct Multipole_t{
	float q;                    // total charge
	float x, y, z;                // center
};

struct Multipole_t *genMultipole(
	int nBasis,                     // the number of basis functions
	const struct GTOBasis_t *gto);  // pointer to basis function information

#endif /* MULTIPOLE_H_ */
