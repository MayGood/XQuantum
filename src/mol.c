/*
 * mol.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <stdio.h>
#include <math.h>
#include "mol.h"
#include "util.h"

int get_nElectron(struct Molecule_t *mol){
	int i;
	int Q = 0;
	for (i = 0; i < mol->nAtom; i++)
		Q += mol->Z[i];
	return Q - mol->Q;
}

struct Molecule_t *readMolecule_XYZ(FILE *inFile) {
	int i;
	int nItem;
	struct Molecule_t *mol;
	char str[256];

	// allocate memory
	mol = (struct Molecule_t *) calloc(1, sizeof(struct Molecule_t));
	if (mol == NULL){
		printf("readMolecule_XYZ: Error - Cannot allocate memory\n");
		exit(-1);
	}

	// reading number of atoms in unit cell
	if (fgets(str, 256, inFile) != str){
		printf("readMolecule_XYZ: Error - Reading number of atoms\n");
		exit(-1);
	}
	nItem = sscanf(str, "%d", &(mol->nAtom));
	if (nItem != 1){
		printf("readMolecule_XYZ: Error - Reading number of atoms\n");
		exit(-1);
	}

	printf("Reading molecule in XYZ file format\n");
	printf("Detected %d atoms\n", mol->nAtom);

	// allocate atomic number and coordinates
	mol->Z = calloc(mol->nAtom, sizeof(int));
	if (mol->Z == NULL){
		printf("readMolecule_XYZ: Error - Cannot allocate atomic number\n");
		exit(-1);
	}
	mol->x = calloc(mol->nAtom, sizeof(float));
	mol->y = calloc(mol->nAtom, sizeof(float));
	mol->z = calloc(mol->nAtom, sizeof(float));
	if (mol->x == NULL || mol->y == NULL || mol->z == NULL){
		printf("readMolecule_XYZ: Error - Cannot allocate coordinates\n");
		exit(-1);
	}

	// reading comments to buffer
	if (fgets(str, 256, inFile) != str){
		printf("readMolecule_XYZ: Error - Reading comment from XYZ file\n");
		exit(-1);
	}

	// reading atomic number and coordinates
	for (i = 0; i < mol->nAtom; i++){
		nItem = fscanf(inFile, "%s %f %f %f", str,
			mol->x + i,
			mol->y + i,
			mol->z + i);
		if (nItem != 4){
			printf("readMolecule_XYZ: Error - Reading coordinate\n");
			exit(-1);
		}
		// parse atom name in short format
		mol->Z[i] = sym2Z(str, SYMB_SHORTNAME);

		// convert angstrom to bohr
		mol->x[i] = (mol->x[i]) * ANGSTROM2BOHR;
		mol->y[i] = (mol->y[i]) * ANGSTROM2BOHR;
		mol->z[i] = (mol->z[i]) * ANGSTROM2BOHR;
	}

	// print molecular structure
	printMolecule_XYZ(mol, stdout);

	// set the total molecular charge to zero by default
	mol->Q = 0;

	return mol;
}

void printMolecule_XYZ(
	const struct Molecule_t * mol,// molecule info
	FILE *fd){                    // output file pointer

	int A;        // atom index;
	char str[64]; // atom name in short format

	fprintf(fd,
		"-------------------------------------------------------------\n"
		"                     Atom Coordinates (Angstroms)            \n"
		"     Atom         X              Y             Z             \n"
		"-------------------------------------------------------------\n");
#define CAPVALUE(a) (fabs(a)<1.0E-6)?0.0:a
	for (A = 0; A < mol->nAtom; A++){
		Z2SymShort(mol->Z[A], str);
		fprintf(fd, "   %5s %15.6f%15.6f%15.6f\n",
			str,
			CAPVALUE(mol->x[A])*BOHR2ANGSTROM,
			CAPVALUE(mol->y[A])*BOHR2ANGSTROM,
			CAPVALUE(mol->z[A])*BOHR2ANGSTROM);
	}
	fprintf(fd,
		"-------------------------------------------------------------\n");
#undef CAPVALUE
}

int get_nEA(struct Molecule_t *mol, int M){
	int N;
	N = get_nElectron(mol);
	if ((N + M - 1) % 2 == 1){
		printf("get_nEA - error number of electron does not match multiplicity\n");
		exit(-1);
	}
	return (N + M - 1) / 2;
}
int get_nEB(struct Molecule_t *mol, int M){
	int N;
	N = get_nElectron(mol);
	if ((N - M + 1) % 2 == 1){
		printf("get_nEB - error number of electron does not match multiplicity\n");
		exit(-1);
	}
	if ((N - M + 1) < 0){
		printf("get_nEB - error multiplicity cannot larger than number of electron\n");
		exit(-1);
	}
	return (N - M + 1) / 2;
}

float nuclei_coulomb(struct Molecule_t *mol){
	float U;
	int i, j;

	U = 0.0;
	for (i = 0; i < mol->nAtom; i++)
	for (j = 0; j < i; j++)
		U += mol->Z[i] * mol->Z[j] /
		sqrt((mol->x[i] - mol->x[j])*(mol->x[i] - mol->x[j]) +
		(mol->y[i] - mol->y[j])*(mol->y[i] - mol->y[j]) +
		(mol->z[i] - mol->z[j])*(mol->z[i] - mol->z[j]));
	return U;
}
