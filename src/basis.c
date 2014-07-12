/*
 * basis.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basis.h"
#include "util.h"
#include "int.h"

#define MAX_BASIS_PER_FILE 20000
struct GTOBasisSet_t *read_GAMESS_BasisSet(FILE *inFile, char *name, int *nBasisSet) {
	int i;
	int Z = 0;
	int nContract;

	// problem is that when reading file we do not know how
	// many basis set is in there before hand.
	struct GTOBasisSet_t * basisDB;     // at returning stage
	int nDBSet = 0;                       // number of db items
	char keyword[1024];                 // current keyword
	float bas_exp[1024];               // buffer for exponent
	float bas_coef_a[1024];            // buffer for coefficient
	float bas_coef_b[1024];            // 2nd buffer if needed

	// allocate memeory at maximum
	basisDB = calloc(MAX_BASIS_PER_FILE, sizeof(struct GTOBasisSet_t));

	// search for word "$DATA"
	if (findf(inFile, 1, "$DATA") == EOF) {
		printf("read_GAMESS_BasisSet - "
			"Error: Cannot find $DATA keyword\n");
		exit(-1);
	}

	do {
		// read keyword
		if (fscanf(inFile, "%s", keyword) != 1) {
			printf("read_GAMESS_BasisSet - "
				"Error: File corrupted\n");
			exit(-1);
		}

		// keyword could be atom name
		if (sym2Z(keyword, SYMB_LONGNAME) > 0) {
			Z = sym2Z(keyword, SYMB_LONGNAME);
			continue;
		}

		// keyword could be orbital index
		if (strcmp(keyword, "S") == 0 ||
			strcmp(keyword, "L") == 0 ||
			strcmp(keyword, "P") == 0 ||
			strcmp(keyword, "D") == 0 ||
			strcmp(keyword, "F") == 0) {

			// read number of contracted function
			if (fscanf(inFile, "%d", &nContract) != 1) {
				printf("read_GAMESS_BasisSet - "
					"Error: No contracted function number\n");
				exit(-1);
			}

			// read in
			for (i = 0; i < nContract; i++) {
				// special for L: both coefficients at the same time
				if (strcmp(keyword, "L") == 0) {
					if (fscanf(inFile, "%*d %f %f %f", bas_exp + i,
						bas_coef_a + i,
						bas_coef_b + i) != 3){
						printf("read_GAMESS_BasisSet - "
							"Error: Reading coefficient\n");
						exit(-1);
					}
				}
				// process normal orbitals
				else if (fscanf(inFile, "%*d %f %f", bas_exp + i,
					bas_coef_a + i) != 2) {
					printf("read_GAMESS_BasisSet - "
						"Error: Reading coefficient\n");
					exit(-1);
				}
			}

#define ALLOCATE												\
	basisDB[nDBSet].nContract = nContract;						\
	basisDB[nDBSet].coef = calloc(nContract, sizeof(float));	\
	basisDB[nDBSet].exp = calloc(nContract, sizeof(float))

#define SETINFO(li, lj, lk, ee, cc)								\
	strcpy(basisDB[nDBSet].name, name);							\
	basisDB[nDBSet].Z = Z;										\
	basisDB[nDBSet].l = li;										\
	basisDB[nDBSet].m = lj;										\
	basisDB[nDBSet].n = lk;										\
			for (i = 0; i < nContract; i++) {					\
				basisDB[nDBSet].coef[i] = cc[i];				\
				basisDB[nDBSet].exp[i] = ee[i];					\
			}

			// process orbital
			if (strcmp(keyword, "S") == 0) {
				ALLOCATE; SETINFO(0, 0, 0, bas_exp, bas_coef_a); nDBSet++;
			}

			// process l orbital
			if (strcmp(keyword, "L") == 0) {
				ALLOCATE; SETINFO(0, 0, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1, 0, 0, bas_exp, bas_coef_b); nDBSet++;
				ALLOCATE; SETINFO(0, 1, 0, bas_exp, bas_coef_b); nDBSet++;
				ALLOCATE; SETINFO(0, 0, 1, bas_exp, bas_coef_b); nDBSet++;
			}

			// process p orbital
			if (strcmp(keyword, "P") == 0) {
				ALLOCATE; SETINFO(1, 0, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 1, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 0, 1, bas_exp, bas_coef_a); nDBSet++;
			}

			// process d orbital
			if (strcmp(keyword, "D") == 0) {
				ALLOCATE; SETINFO(2, 0, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 2, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 0, 2, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1, 1, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1, 0, 1, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 1, 1, bas_exp, bas_coef_a); nDBSet++;
			}

			// process f orbital
			if (strcmp(keyword, "F") == 0) {
				ALLOCATE; SETINFO(3, 0, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 3, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 0, 3, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(2, 1, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1, 2, 0, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(2, 0, 1, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1, 0, 2, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 2, 1, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0, 1, 2, bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1, 1, 1, bas_exp, bas_coef_a); nDBSet++;
			}
#undef SETINFO
#undef ALLOCATE
			continue;
		}
	} while (strcmp(keyword, "$END") != 0);

	// truckate memory
	basisDB = realloc(basisDB, sizeof(struct GTOBasisSet_t) * nDBSet);

	// output and return value
	*nBasisSet = nDBSet;
	return basisDB;
}

struct GTOBasis_t * genBasis(
	struct Molecule_t *mol,               // pointer to molecular structure
	int *nBasis,                          // return the number of basis created
	int dbSize,                           // number of record in basis set database
	const struct GTOBasisSet_t *basisDB){ // pointer to basis set database

	int i, b, j, n = 0;
	struct GTOBasis_t *gto;
	float norm;

	// allocate memory
	gto = calloc(MAX_BASIS, sizeof(struct GTOBasis_t));
	if (gto == NULL){
		printf("genBasis: Error - Cannot allocate memory\n");
		exit(EXIT_FAILURE);
	}

	// loop thru all atoms
	for (i = 0; i < mol->nAtom; i++){
		// make sure we have at least one matching
		for (b = 0; b < dbSize; b++) {
			if (basisDB[b].Z == mol->Z[i])
				break;
		}
		if (b == dbSize){
			printf("genBasis: Error - Cannot find atom %d in basis file.\n",
				mol->Z[i]);
			exit(EXIT_FAILURE);
		}
		// find all matched atomic number
		for (b = 0; b < dbSize; b++){
			if (basisDB[b].Z == mol->Z[i]){
				// copy values
				gto[n].nContract = basisDB[b].nContract;
				gto[n].l = basisDB[b].l;
				gto[n].m = basisDB[b].m;
				gto[n].n = basisDB[b].n;
				gto[n].x0 = mol->x[i];
				gto[n].y0 = mol->y[i];
				gto[n].z0 = mol->z[i];
				// allocate memory
//				gto[n].exp = calloc(gto[n].nContract, sizeof(float));
//				gto[n].coef = calloc(gto[n].nContract, sizeof(float));
//				gto[n].norm = calloc(gto[n].nContract, sizeof(float));
//				if (gto[n].exp == NULL ||
//					gto[n].coef == NULL ||
//					gto[n].norm == NULL){
//					printf("genBasis: Error - Cannot allocate memory\n");
//					exit(EXIT_FAILURE);
//				}
				// copy values
				for (j = 0; j < gto[n].nContract; j++){
					gto[n].exp[j] = basisDB[b].exp[j];
					gto[n].coef[j] = basisDB[b].coef[j];
				}
				n++;
			}
		}
	}

	// shrink memory
	gto = realloc(gto, sizeof(struct GTOBasis_t) * n);
	*nBasis = n;

	//
	// normalize contracted coefficient
	//
	// loop over basis set
	for (i = 0; i < *nBasis; i++){
		// loop over contracted functions
		for (j = 0; j < gto[i].nContract; j++){
			norm = overlap(gto[i].exp[j],
				gto[i].l,
				gto[i].m,
				gto[i].n,
				gto[i].x0,
				gto[i].y0,
				gto[i].z0,

				gto[i].exp[j],
				gto[i].l,
				gto[i].m,
				gto[i].n,
				gto[i].x0,
				gto[i].y0,
				gto[i].z0);

			gto[i].norm[j] = 1.0 / sqrt(norm);
		}
	}

	// report stat
	printf("There are %d Cartesian-Gaussian contracted basis functions\n", n);
	i = get_nPrim(*nBasis, gto);
	printf("There are %d Primitive GTO functions\n", i);

	return gto;
}

int get_nPrim(int nBasis, const struct GTOBasis_t *gto){
	int i;
	int nPrim = 0;

	for (i = 0; i < nBasis; i++)
		nPrim += gto[i].nContract;
	return nPrim;
}

