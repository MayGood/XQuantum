/*
 * multipole.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "multipole.h"
#include "basis.h"
#include "matrix.h"

//
// genMultipole : generate multipole moments for a product of basis function.
// This can be used to approximate coulomb integral.
//
// 23 Dec, 2012 - Teepanis Chachiyo
//       Initial implementation and testing
//
struct Multipole_t *genMultipole(
	int nBasis,                     // the number of basis functions
	const struct GTOBasis_t *gto){  // pointer to basis function information

	struct Multipole_t *mpole;  // pointer to multipole data structure
	float q;                   // total charge
	float px, py, pz;            // dipole moment
	float x, y, z;               // center
	int i, j;                    // basis function index

	// allocate memory
	mpole = (struct Multipole_t *)
		calloc(nBasis*nBasis, sizeof(struct Multipole_t));
	if (mpole == NULL){
		printf("genMultipole: Error - cannot allocate memory\n");
		exit(-1);
	}

	// loop thru all basis functions
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){

		// compute total charge
		q = GTO_overlap(i, j, gto);

		// truncate total charge
		if (fabs(q) < MULTIPOLE_CHARGE_CUTOFF) q = 0.0;

		// compute charge center
		x = (gto[i].x0 + gto[j].x0) / 2.0;
		y = (gto[i].y0 + gto[j].y0) / 2.0;
		z = (gto[i].z0 + gto[j].z0) / 2.0;

		// compute dipole moment
		px = GTO_moment(i, j, gto, 1, 0, 0, x, y, z);
		py = GTO_moment(i, j, gto, 0, 1, 0, x, y, z);
		pz = GTO_moment(i, j, gto, 0, 0, 1, x, y, z);

		// re-evaluate new center so that dipole becomes zero
		if (q != 0.0){
			x = x + px / q;
			y = y + py / q;
			z = z + pz / q;
		}

		// store values
		mpole[i*nBasis + j].q = q;
		mpole[i*nBasis + j].x = x;
		mpole[i*nBasis + j].y = y;
		mpole[i*nBasis + j].z = z;
	}

	return mpole;
}
