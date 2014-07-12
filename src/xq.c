/*
 * qc.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "basis.h"
#include "mol.h"
#include "option.h"
#include "util.h"
#include "uhf.h"

static void greeting(FILE *outFile) {
	char time[256];

	time_str(256, time);
	fprintf(outFile, "Thanks to SQ - XQuantum\n"
		"Begin: %s\n", time);
}

int main(int argc, char *argv[]) {
	struct Molecule_t *mol = NULL;			// pointer to molecular structure info
	struct GTOBasis_t *gto = NULL;			// pointer to basis function
	struct option_t opt;					// store user specifed options
	struct GTOBasisSet_t *dbSet = NULL;		// pointer to basis set information
	float *CA, *eA, *CB, *eB;				// molecular orbitals and eigen values
	int nBasis;				// number of basis function
	int dbSetItem;			// number of basis set item
	FILE *inFile;			// xyz file pointer
	FILE *basisFile;		// basis function file pointer
	float Etot;			// total energy
	char time[256];
	struct timeval startTime;
	struct timeval endTime;

	gettimeofday(&startTime, NULL);

	greeting(stdout);

	if (argc < 3) {
		exit(EXIT_FAILURE);
	}

	parse_option(&opt, argc, argv);

	inFile = fopen(argv[1], "r");
	if (inFile == NULL) {
		printf("Cannot open file %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	mol = readMolecule_XYZ(inFile);

	////////////////////////////////////
	// handle default for multiplicity
	////////////////////////////////////
	// even number of electron set multiplicity to one by default
	if (opt.multiplicity == 0 && get_nElectron(mol) % 2 == 0)
		opt.multiplicity = 1;
	// odd number of electron set multiplicity to two by default
	if (opt.multiplicity == 0 && get_nElectron(mol) % 2 == 1)
		opt.multiplicity = 2;

	/////////////////////////////////////
	// handle default for RHF and UHF
	/////////////////////////////////////
	if (opt.RHF == 0 && opt.UHF == 0) {
		if (opt.multiplicity == 1) opt.RHF = 1; else opt.UHF = 1;
	}

	basisFile = fopen(argv[2], "r");
	if (basisFile == NULL) {
		printf("Cannot open file %s\n", argv[2]);
		exit(EXIT_FAILURE);
	}
	dbSet = read_GAMESS_BasisSet(basisFile, argv[2], &dbSetItem);

	gto = genBasis(mol, &nBasis, dbSetItem, dbSet);

	// allocate meory for molecular orbitals and their eigen values
	CA = calloc(nBasis * nBasis, sizeof(float));
	eA = calloc(nBasis, sizeof(float));
	CB = calloc(nBasis * nBasis, sizeof(float));
	eB = calloc(nBasis, sizeof(float));
	if (CA == NULL || eA == NULL || CB == NULL || eB == NULL){
		printf("main: error - cannot allocate memory\n");
		exit(-1);
	}

	Etot = uhf(nBasis, gto, mol,
		get_nEA(mol, opt.multiplicity), get_nEB(mol, opt.multiplicity),
		CA, CB, eA, eB, &opt);

	// successful calculations
	gettimeofday(&endTime, NULL);
	printf("Job cost time: %ld us.\n", 1000000 * (endTime.tv_sec - startTime.tv_sec) + endTime.tv_usec - startTime.tv_usec);
	return 0;
}
