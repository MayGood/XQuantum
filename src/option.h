/*
 * option.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef OPTION_H_
#define OPTION_H_

struct option_t {
	int molCharge;
	int multiplicity;

	int outVolumeType;
#define VOLUME_NONE				0
#define VOLUME_DENSITY_TOTAL	1
#define VOLUME_DENSITY_SPIN		2
#define VOLUME_MO_ALPHA			3
#define VOLUME_MO_BETA			4
#define VOLUME_POTENTIAL		5

	int outWhichMO;
	int outFormat;
#define VOLUME_CUBE				1
#define VOLUME_XSF				2

	float outVolumeCut;
	int outGAUSSIAN;

	float SCFConv;
	float SCFDrag;
	float SCFCutoff;
	int SCFMax;
#define CONVMETHOD_DAMPING		0
#define CONVMETHOD_DIIS2		1
#define CONVMETHOD_DIIS3		2
#define CONVMETHOD_DIIS4		3
	int convMethod;
	int SCFAccuracy;
#define SCFACCURACY_3STEP		0
#define SCFACCURACY_1STEP		1

#define SCFACCURACY_3STEP_A		1.0E-3
#define SCFACCURACY_3STEP_B		1.0E-6
#define SCFACCURACY_3STEP_C		1.0E-9

	int SCFGuess;
#define SCFGUESS_DIAG			0
#define SCFGUESS_CORE			1
#define SCFGUESS_CHECK			2

	int maxMem;
	int loadDMatrix;
	int saveDMatrix;
	int saveCheck;
	int saveCheckAll;
	int loadCheck;
	int opt;
	int optMax;
	int RHF;
	int MP2;
	int UHF;
	int force;
	char DMatrixFile[256];
	char CheckFile[256];
	int nCPU;
	char prefixStr[256];
	int MECP;
	int mecpMax;
	int mecpMA;
	int mecpMB;
	char DMatrixFileA[256];
	char DMatrixFileB[256];
	char CheckFileA[256];
	char CheckFileB[256];
	char gaussEXE[256];
	char gaussINA[256];
	char gaussINB[256];
};
void parse_option(struct option_t *opt, int argc, char *argv[]);

#endif /* OPTION_H_ */
