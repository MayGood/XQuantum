/*
 * basis.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef BASIS_H_
#define BASIS_H_

#include "mol.h"

#define MAXCONTRACT 16
struct GTOBasis_t {
	int nContract;
	int l;
	int m;
	int n;
	float x0;
	float y0;
	float z0;
	float coef[MAXCONTRACT];
	float exp[MAXCONTRACT];
	float norm[MAXCONTRACT];
};

#define MAX_BASIS 10000
#define BASIS_NAME 256

struct GTOBasisSet_t {
	char name[BASIS_NAME];
	int Z;
	int nContract;
	int l;
	int m;
	int n;
	float *exp;
	float *coef;
};

struct GTOBasis_t *cleanGTOBasis(struct GTOBasis_t *gtoBasis, int nBasis);
struct GTOBasisSet_t *read_GAMESS_BasisSet(FILE *inFile, char *name, int *nBasisSet);

struct GTOBasis_t *genBasis(struct Molecule_t *mol, int *nBasis, int dbSize, const struct GTOBasisSet_t *basisDB);

int get_nPrim(int nBasis, const struct GTOBasis_t *gto);

#endif /* BASIS_H_ */
