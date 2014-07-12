/*
 * mol.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef MOL_H_
#define MOL_H_

struct Molecule_t {
	int nAtom;	// number of nuclei
	int Q;		// total molecular charge
	int *Z;		// atomic number array
	float *x;		// Cartesian x-coordinate
	float *y;		// Cartesian y-coordinate
	float *z;		// Cartesian z-coordinate
};

void printMolecule_XYZ(const struct Molecule_t *mol, FILE *fd);
struct Molecule_t *readMolecule_XYZ(FILE *inFile);
struct Molecule_t *cleanMolecule(struct Molecule_t *mol);
int get_nElectron(struct Molecule_t *mol);
float nuclei_coulomb(struct Molecule_t *mol);
int get_nEA(struct Molecule_t *mol, int M);
int get_nEB(struct Molecule_t *mol, int M);

#endif /* MOL_H_ */
