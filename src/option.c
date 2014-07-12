/*
 * option.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "option.h"

void parse_option(struct option_t *opt, int argc, char *argv[]) {
	int i;

	// set default
	opt->molCharge = 0;
	opt->multiplicity = 0;
	opt->RHF = 0;
	opt->MP2 = 0;
	opt->UHF = 0;
	opt->force = 0;
	opt->opt = 0;
	opt->outVolumeType = VOLUME_NONE;
	opt->outWhichMO = 0;
	opt->outFormat = VOLUME_XSF;
	opt->outVolumeCut = 1.0E-4;
	opt->outGAUSSIAN = 0;
	opt->SCFGuess = SCFGUESS_DIAG;
	opt->SCFConv = 1.0E-6;
	opt->SCFCutoff = 1.0E-15;
	opt->SCFDrag = 0.25;
	opt->SCFMax = 80;
	opt->convMethod = CONVMETHOD_DIIS4;
	opt->SCFAccuracy = SCFACCURACY_3STEP;
	opt->maxMem = 250;
	opt->loadDMatrix = 0;
	opt->saveDMatrix = 0;
	strcpy(opt->DMatrixFile, "dmatrix.txt");
	opt->saveCheck = 0;
	opt->saveCheckAll = 0;
	opt->loadCheck = 0;
	strcpy(opt->CheckFile, "checkpoint.txt");
	opt->nCPU = 1;
	strcpy(opt->prefixStr, "SQ");
	opt->MECP = 0;
	opt->mecpMax = 30;
	opt->mecpMA = 0;
	opt->mecpMB = 0;
	strcpy(opt->DMatrixFileA, "dmatrixA.txt");
	strcpy(opt->DMatrixFileB, "dmatrixB.txt");
	strcpy(opt->CheckFileA, "checkpointA.txt");
	strcpy(opt->CheckFileB, "checkpointB.txt");
	strcpy(opt->gaussEXE, "\0");
	strcpy(opt->gaussINA, "\0");
	strcpy(opt->gaussINB, "\0");
	opt->optMax = 30;

	// loop throu all options
	for (i = 3; i < argc; i++) {
		if (strncmp(argv[i], "-Q=", 3) == 0) {
			opt->molCharge = atoi(argv[i] + 3);
			continue;
		}
		if (strncmp(argv[i], "-M=", 3) == 0) {
			opt->multiplicity = atoi(argv[i] + 3);
			continue;
		}
		if (strcmp(argv[i], "-RHF") == 0) {
			opt->RHF = 1;
			continue;
		}
		if (strcmp(argv[i], "-MP2") == 0) {
			opt->MP2 = 1;
			continue;
		}
		if (strcmp(argv[i], "-UHF") == 0) {
			opt->UHF = 1;
			continue;
		}
		if (strcmp(argv[i], "-FORCE") == 0) {
			opt->force = 1;
			continue;
		}
		if (strcmp(argv[i], "-OPT") == 0) {
			opt->opt = 1;
			continue;
		}
		if (strncmp(argv[i], "-OPTMAX=", 8) == 0) {
			opt->optMax = atoi(argv[i] + 8);
			continue;
		}
		if (strcmp(argv[i], "-DENSITY") == 0) {
			if (opt->outVolumeType != VOLUME_NONE) {
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_DENSITY_TOTAL;
			continue;
		}
		if (strcmp(argv[i], "-POTENTIAL") == 0) {
			if (opt->outVolumeType != VOLUME_NONE) {
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_POTENTIAL;
			continue;
		}
		if (strncmp(argv[i], "-MOALPHA=", 9) == 0) {
			if (opt->outVolumeType != VOLUME_NONE) {
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_MO_ALPHA;
			opt->outWhichMO = atoi(argv[i] + 9);
			continue;
		}
		if (strncmp(argv[i], "-MOBETA=", 8) == 0) {
			if (opt->outVolumeType != VOLUME_NONE) {
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_MO_BETA;
			opt->outWhichMO = atoi(argv[i] + 8);
			continue;
		}
		if (strncmp(argv[i], "-VOLCUT=", 8) == 0) {
			opt->outVolumeCut = atof(argv[i] + 8);
			continue;
		}
		if (strcmp(argv[i], "-GAUSSIAN") == 0) {
			opt->outGAUSSIAN = 1;
			continue;
		}
		if (strcmp(argv[i], "-XSF") == 0) {
			opt->outFormat = VOLUME_XSF;
			continue;
		}
		if (strcmp(argv[i], "-CUBE") == 0) {
			opt->outFormat = VOLUME_CUBE;
			continue;
		}
		if (strcmp(argv[i], "-GUESS=DIAG") == 0) {
			opt->SCFGuess = SCFGUESS_DIAG;
			continue;
		}
		if (strcmp(argv[i], "-GUESS=CORE") == 0) {
			opt->SCFGuess = SCFGUESS_CORE;
			continue;
		}
		if (strcmp(argv[i], "-GUESS=CHECK") == 0) {
			opt->SCFGuess = SCFGUESS_CHECK;
			continue;
		}
		if (strcmp(argv[i], "-SCFACC=3STEP") == 0) {
			opt->SCFAccuracy = SCFACCURACY_3STEP;
			continue;
		}
		if (strcmp(argv[i], "-SCFACC=1STEP") == 0) {
			opt->SCFAccuracy = SCFACCURACY_1STEP;
			continue;
		}
		if (strcmp(argv[i], "-SCFDIIS") == 0) {
			opt->convMethod = CONVMETHOD_DIIS4;
			continue;
		}
		if (strcmp(argv[i], "-SCFDIIS3") == 0) {
			opt->convMethod = CONVMETHOD_DIIS3;
			continue;
		}
		if (strcmp(argv[i], "-SCFDIIS2") == 0) {
			opt->convMethod = CONVMETHOD_DIIS2;
			continue;
		}
		if (strcmp(argv[i], "-SCFDAMP") == 0) {
			opt->convMethod = CONVMETHOD_DAMPING;
			continue;
		}
		if (strncmp(argv[i], "-MAXMEM=", 8) == 0) {
			opt->maxMem = atoi(argv[i] + 8);
			continue;
		}
		if (strncmp(argv[i], "-SCFCONV=", 9) == 0) {
			opt->SCFConv = atof(argv[i] + 9);
			continue;
		}
		if (strncmp(argv[i], "-SCFDRAG=", 9) == 0) {
			opt->SCFDrag = atof(argv[i] + 9);
			continue;
		}
		if (strncmp(argv[i], "-SCFMAX=", 8) == 0) {
			opt->SCFMax = atoi(argv[i] + 8);
			continue;
		}
		if (strcmp(argv[i], "-LDMATRIX") == 0) {
			opt->loadDMatrix = 1;
			continue;
		}
		if (strcmp(argv[i], "-SDMATRIX") == 0) {
			opt->saveDMatrix = 1;
			continue;
		}
		if (strncmp(argv[i], "-FDMATRIX=", 10) == 0) {
			strcpy(opt->DMatrixFile, argv[i] + 10);
			continue;
		}
		if (strcmp(argv[i], "-SCHECK") == 0) {
			opt->saveCheck = 1;
			continue;
		}
		if (strcmp(argv[i], "-SCHECK=ALL") == 0) {
			opt->saveCheckAll = 1;
			continue;
		}
		if (strcmp(argv[i], "-LCHECK") == 0) {
			opt->loadCheck = 1;
			continue;
		}
		if (strncmp(argv[i], "-FCHECK=", 8) == 0) {
			strcpy(opt->CheckFile, argv[i] + 8);
			continue;
		}
		if (strncmp(argv[i], "-NCPU=", 6) == 0) {
			opt->nCPU = atoi(argv[i] + 6);
			continue;
		}
		if (strncmp(argv[i], "-PREFIX=", 8) == 0) {
			strcpy(opt->prefixStr, argv[i] + 8);
			continue;
		}
		if (strncmp(argv[i], "-MECP=", 6) == 0) {
			if (sscanf(argv[i] + 6, "%d,%d", &opt->mecpMA, &opt->mecpMB) != 2) {
				printf("parse_option - error cannot recognize option %s\n", argv[i]);
				exit(-1);
			}
			opt->MECP = 1;
			continue;
		}
		if (strncmp(argv[i], "-MECPMAX=", 9) == 0) {
			opt->mecpMax = atoi(argv[i] + 9);
			continue;
		}
		if (strncmp(argv[i], "-FDMATRIXA=", 11) == 0) {
			strcpy(opt->DMatrixFileA, argv[i] + 11);
			continue;
		}
		if (strncmp(argv[i], "-FDMATRIXB=", 11) == 0) {
			strcpy(opt->DMatrixFileB, argv[i] + 11);
			continue;
		}
		if (strncmp(argv[i], "-FCHECKA=", 9) == 0) {
			strcpy(opt->CheckFileA, argv[i] + 9);
			continue;
		}
		if (strncmp(argv[i], "-FCHECKB=", 9) == 0) {
			strcpy(opt->CheckFileB, argv[i] + 9);
			continue;
		}
		if (strncmp(argv[i], "-GAUSSEXE=", 10) == 0) {
			strcpy(opt->gaussEXE, argv[i] + 10);
			continue;
		}
		if (strncmp(argv[i], "-GAUSSINA=", 10) == 0) {
			strcpy(opt->gaussINA, argv[i] + 10);
			continue;
		}
		if (strncmp(argv[i], "-GAUSSINB=", 10) == 0) {
			strcpy(opt->gaussINB, argv[i] + 10);
			continue;
		}

		// cannot recongnize parameter
		printf("parse_option - error cannot recognize option %s\n", argv[i]);
		exit(-1);
	}

	// validate
	if (opt->outVolumeType == VOLUME_MO_ALPHA || opt->outVolumeType == VOLUME_MO_BETA) {
		if (opt->outWhichMO <= 0) {
			printf("parse_option - error molecular orbital index should be greater than zero\n");
			exit(-1);
		}
	}

	// validate SCFConv
	if (opt->SCFConv <= 0.0) {
		printf("parse_option - error invalid SCFConv range\n");
		exit(-1);
	}

	// validate volumeCut
	if (opt->outVolumeCut <= 0.0) {
		printf("parse_option - error invalid volumeCut range\n");
		exit(-1);
	}

	// validate SCFDrag
	if (opt->SCFDrag <= 0.0 || opt->SCFDrag > 1) {
		printf("parse_option - error invalid SCFDrag range\n");
		exit(-1);
	}

	// validate SCFMax
	if (opt->SCFMax < 0) {
		printf("parse_option - error invalid SCFMax range\n");
		exit(-1);
	}

	// validate maxMem
	if (opt->maxMem < 0) {
		printf("parse_option - error invalid MAXMEM range\n");
		exit(-1);
	}

	// validate RHF and UHF choice
	if (opt->RHF + opt->UHF > 1) {
		printf("parse_option - error can not choose both RHF and UHF\n");
		exit(-1);
	}

	// validate checkpoint related options
	if (opt->loadCheck && opt->opt) {
		printf("parse_option - error OPT cannot be used with LCHECK\n");
		exit(-1);
	}
	if (opt->loadCheck && (opt->saveCheck || opt->saveCheckAll)) {
		printf("parse_option - error LCHECK cannot be used with saving checkpoints\n");
		exit(-1);
	}
	if (opt->saveCheck && opt->saveCheckAll) {
		printf("parse_option - error SCHECK cannot be used with SCHECK=ALL");
		exit(-1);
	}

	// optimization does not support mp2 yet
	if (opt->opt && opt->MP2) {
		printf("parse_option - error OPT does not support MP2 at the moment\n");
		exit(-1);
	}

	// validate nCPU
	if (opt->nCPU <= 0) {
		printf("parse_option - error invalid number of cpus\n");
		exit(-1);
	}

	// validate MECP
	if (opt->MECP) {
		if (opt->opt){
			printf("parse_option - error MECP cannot be used with OPT\n");
			exit(-1);
		}
		if (opt->MP2){
			printf("parse_option - error MECP cannot be used with MP2\n");
			exit(-1);
		}
	}
}

