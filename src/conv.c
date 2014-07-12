/*
 * conv.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "conv.h"
#include "lin.h"

float conv_damping(int nBasis,     // number of basis
	float alpha,   // drag coefficient between 0 and 1
	float *PA,     // alpha density matrix
	float *PB){    // beta density matrix

	static float *prevPA = NULL;  // previous iteration alpha density matrix
	static float *prevPB = NULL;  // previous iteration beta density matrix
	static int nIter = 0;          // number of iteration
	int i, j;                     // generic loop indexes
	float avgdP;                // average deviation

	// reset if requested
	if (nBasis == 0 && alpha == 0.0 && PA == NULL && PB == NULL){
		nIter = 0;
		if (prevPA != NULL){ free(prevPA); prevPA = NULL; }
		if (prevPB != NULL){ free(prevPB); prevPB = NULL; }
		return 0.0;
	}

	// check range of alpha
	if (alpha < 0.0 || alpha > 1.0){
		printf("conv_damping : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix
	////////////////////////////////////////////////////////
	if (nIter == 0){

		// check if array is NULL
		if (prevPA != NULL || prevPB != NULL){
			printf("conv_damping : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		prevPA = calloc(nBasis*nBasis, sizeof(float));
		prevPB = calloc(nBasis*nBasis, sizeof(float));
		if (prevPA == NULL || prevPB == NULL){
			printf("conv_damping : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			prevPA[i*nBasis + j] = PA[i*nBasis + j];
			prevPB[i*nBasis + j] = PB[i*nBasis + j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	//////////////////////////////////////////////////////////////////
	/// for subsequent iteration we first compute average deviation
	//////////////////////////////////////////////////////////////////
	avgdP = 0.0;
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		avgdP += (PA[i*nBasis + j] - prevPA[i*nBasis + j])
			*(PA[i*nBasis + j] - prevPA[i*nBasis + j]);
		avgdP += (PB[i*nBasis + j] - prevPB[i*nBasis + j])
			*(PB[i*nBasis + j] - prevPB[i*nBasis + j]);
	}
	avgdP = sqrt(avgdP / nBasis / nBasis / 2.0);

	// update
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		PA[i*nBasis + j] = prevPA[i*nBasis + j] +
			alpha *(PA[i*nBasis + j] - prevPA[i*nBasis + j]);

		PB[i*nBasis + j] = prevPB[i*nBasis + j] +
			alpha *(PB[i*nBasis + j] - prevPB[i*nBasis + j]);
	}

	// store previous density matrix
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		prevPA[i*nBasis + j] = PA[i*nBasis + j];
		prevPB[i*nBasis + j] = PB[i*nBasis + j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}

float conv_diis2(int nBasis,      // number of basis
	float alpha,    // drag coefficient between 0 and 1
	float *PA,      // output alpha density matrix
	float *PB){     // output beta  density matrix

	static float *PA0 = NULL;     // input alpha density matrix
	static float *PB0 = NULL;     // input beta  density matrix
	static float *PA1 = NULL;     // input alpha density matrix
	static float *PB1 = NULL;     // input beta  density matrix
	static float *dPA0 = NULL;    // error alpha density matrix
	static float *dPB0 = NULL;    // error beta  density matrix
	static float *dPA1 = NULL;    // error alpha density matrix
	static float *dPB1 = NULL;    // error beta  density matrix
	static int nIter = 0;          // number of iteration
	int i, j;                     // generic loop indexes
	float avgdP;                // average deviation
	float c1, c0, B00, B01, B11;    // diis variables

	// reset if requested
	if (nBasis == 0 && alpha == 0.0 && PA == NULL && PB == NULL){
		nIter = 0;
		if (PA0 != NULL){ free(PA0); PA0 = NULL; }
		if (PB0 != NULL){ free(PB0); PB0 = NULL; }
		if (PA1 != NULL){ free(PA1); PA1 = NULL; }
		if (PB1 != NULL){ free(PB1); PB1 = NULL; }

		if (dPA0 != NULL){ free(dPA0); dPA0 = NULL; }
		if (dPB0 != NULL){ free(dPB0); dPB0 = NULL; }
		if (dPA1 != NULL){ free(dPA1); dPA1 = NULL; }
		if (dPB1 != NULL){ free(dPB1); dPB1 = NULL; }
		return 0.0;
	}

	// check range of alpha
	if (alpha < 0.0 || alpha > 1.0){
		printf("conv_diis2 : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix
	////////////////////////////////////////////////////////
	if (nIter == 0){

		// check if array is NULL
		if (PA0 != NULL || PB0 != NULL ||
			PA1 != NULL || PB1 != NULL ||
			dPA0 != NULL || dPB0 != NULL ||
			dPA1 != NULL || dPB1 != NULL){
			printf("conv_diis2 : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		PA0 = calloc(nBasis*nBasis, sizeof(float));
		PB0 = calloc(nBasis*nBasis, sizeof(float));
		PA1 = calloc(nBasis*nBasis, sizeof(float));
		PB1 = calloc(nBasis*nBasis, sizeof(float));
		dPA0 = calloc(nBasis*nBasis, sizeof(float));
		dPB0 = calloc(nBasis*nBasis, sizeof(float));
		dPA1 = calloc(nBasis*nBasis, sizeof(float));
		dPB1 = calloc(nBasis*nBasis, sizeof(float));
		if (PA0 == NULL || PB0 == NULL ||
			PA1 == NULL || PB1 == NULL ||
			dPA0 == NULL || dPB0 == NULL ||
			dPA1 == NULL || dPB1 == NULL){
			printf("conv_diis2 : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA1[i*nBasis + j] = PA[i*nBasis + j];
			PB1[i*nBasis + j] = PB[i*nBasis + j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	////////////////////////////////////////////////////////
	/// subsequent iteration
	////////////////////////////////////////////////////////

	// compute rms deviation
	avgdP = 0.0;
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		avgdP += (PA[i*nBasis + j] - PA1[i*nBasis + j])
			*(PA[i*nBasis + j] - PA1[i*nBasis + j]);
		avgdP += (PB[i*nBasis + j] - PB1[i*nBasis + j])
			*(PB[i*nBasis + j] - PB1[i*nBasis + j]);
	}
	avgdP = sqrt(avgdP / nBasis / nBasis / 2.0);

	// compute dPA1 and dPB1
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		dPA1[i*nBasis + j] = PA[i*nBasis + j] - PA1[i*nBasis + j];
		dPB1[i*nBasis + j] = PB[i*nBasis + j] - PB1[i*nBasis + j];
	}

	switch (nIter){

	case 1: // update using simple damping for the first iteration
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA[i*nBasis + j] = alpha * PA[i*nBasis + j] +
				(1.0 - alpha) * PA1[i*nBasis + j];

			PB[i*nBasis + j] = alpha * PB[i*nBasis + j] +
				(1.0 - alpha) * PB1[i*nBasis + j];
		}
		break;

	default: // using 2-point DIIS for the rest

		// compute diis variables
		B00 = B11 = B01 = 0.0;

		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
#define DIIS_B(X,Y) B##X##Y += dPA##X[i*nBasis+j]*dPA##Y[i*nBasis+j]; \
	B##X##Y += dPB##X[i*nBasis + j] * dPB##Y[i*nBasis + j];

			DIIS_B(0, 0); DIIS_B(1, 1); DIIS_B(0, 1);
		}
		c1 = (B00 - B01) / (B11 - B01);
		c0 = 1.0 / (1.0 + c1);
		c1 = c1*c0;

		// use simple damping if DIIS fails
		if (B11>B00){
			c1 = alpha;
			c0 = 1.0 - alpha;
		}

		// apply DIIS updates
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA[i*nBasis + j] = c1*(PA[i*nBasis + j]) +
				c0*(PA0[i*nBasis + j] + dPA0[i*nBasis + j]);
			PB[i*nBasis + j] = c1*(PB[i*nBasis + j]) +
				c0*(PB0[i*nBasis + j] + dPB0[i*nBasis + j]);
		}
		break;

	}

	// store previous density matrix
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		dPA0[i*nBasis + j] = dPA1[i*nBasis + j];
		dPB0[i*nBasis + j] = dPB1[i*nBasis + j];

		PA0[i*nBasis + j] = PA1[i*nBasis + j];
		PB0[i*nBasis + j] = PB1[i*nBasis + j];

		PA1[i*nBasis + j] = PA[i*nBasis + j];
		PB1[i*nBasis + j] = PB[i*nBasis + j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}


//
// conv_diis3 : is a the 3-point version of the conv_diis2
//
// Oct, 2012 - Teepanis Chachiyo
//      Initial implementation and testing
//
float conv_diis3(int nBasis,      // number of basis
	float alpha,    // drag coefficient between 0 and 1
	float *PA,      // output alpha density matrix
	float *PB){     // output beta  density matrix

	static float *PA0 = NULL;        // input alpha density matrix
	static float *PB0 = NULL;        // input beta  density matrix
	static float *PA1 = NULL;        // input alpha density matrix
	static float *PB1 = NULL;        // input beta  density matrix
	static float *PA2 = NULL;        // input alpha density matrix
	static float *PB2 = NULL;        // input beta  density matrix
	static float *dPA0 = NULL;       // error alpha density matrix
	static float *dPB0 = NULL;       // error beta  density matrix
	static float *dPA1 = NULL;       // error alpha density matrix
	static float *dPB1 = NULL;       // error beta  density matrix
	static float *dPA2 = NULL;       // error alpha density matrix
	static float *dPB2 = NULL;       // error beta  density matrix
	static int nIter = 0;             // number of iteration
	int i, j;                        // generic loop indexes
	float avgdP;                   // average deviation
	float c2, c1, c0;                // diis coefficients
	float B00, B11, B22, B01, B02, B12; // diis variables

	// reset if requested
	if (nBasis == 0 && alpha == 0.0 && PA == NULL && PB == NULL){
		nIter = 0;
		if (PA0 != NULL){ free(PA0); PA0 = NULL; }
		if (PB0 != NULL){ free(PB0); PB0 = NULL; }
		if (PA1 != NULL){ free(PA1); PA1 = NULL; }
		if (PB1 != NULL){ free(PB1); PB1 = NULL; }
		if (PA2 != NULL){ free(PA2); PA2 = NULL; }
		if (PB2 != NULL){ free(PB2); PB2 = NULL; }

		if (dPA0 != NULL){ free(dPA0); dPA0 = NULL; }
		if (dPB0 != NULL){ free(dPB0); dPB0 = NULL; }
		if (dPA1 != NULL){ free(dPA1); dPA1 = NULL; }
		if (dPB1 != NULL){ free(dPB1); dPB1 = NULL; }
		if (dPA2 != NULL){ free(dPA2); dPA2 = NULL; }
		if (dPB2 != NULL){ free(dPB2); dPB2 = NULL; }
		return 0.0;
	}

	// check range of alpha
	if (alpha < 0.0 || alpha > 1.0){
		printf("conv_diis3 : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix
	////////////////////////////////////////////////////////
	if (nIter == 0){

		// check if array is NULL
		if (PA0 != NULL || PB0 != NULL ||
			PA1 != NULL || PB1 != NULL ||
			PA2 != NULL || PB2 != NULL ||
			dPA0 != NULL || dPB0 != NULL ||
			dPA1 != NULL || dPB1 != NULL ||
			dPA2 != NULL || dPB2 != NULL){
			printf("conv_diis3 : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		PA0 = calloc(nBasis*nBasis, sizeof(float));
		PB0 = calloc(nBasis*nBasis, sizeof(float));
		PA1 = calloc(nBasis*nBasis, sizeof(float));
		PB1 = calloc(nBasis*nBasis, sizeof(float));
		PA2 = calloc(nBasis*nBasis, sizeof(float));
		PB2 = calloc(nBasis*nBasis, sizeof(float));
		dPA0 = calloc(nBasis*nBasis, sizeof(float));
		dPB0 = calloc(nBasis*nBasis, sizeof(float));
		dPA1 = calloc(nBasis*nBasis, sizeof(float));
		dPB1 = calloc(nBasis*nBasis, sizeof(float));
		dPA2 = calloc(nBasis*nBasis, sizeof(float));
		dPB2 = calloc(nBasis*nBasis, sizeof(float));
		if (PA0 == NULL || PB0 == NULL ||
			PA1 == NULL || PB1 == NULL ||
			PA2 == NULL || PB2 == NULL ||
			dPA0 == NULL || dPB0 == NULL ||
			dPA1 == NULL || dPB1 == NULL ||
			dPA2 == NULL || dPB2 == NULL){
			printf("conv_diis3 : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA2[i*nBasis + j] = PA[i*nBasis + j];
			PB2[i*nBasis + j] = PB[i*nBasis + j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	////////////////////////////////////////////////////////
	/// subsequent iteration
	////////////////////////////////////////////////////////

	// compute rms deviation
	avgdP = 0.0;
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		avgdP += (PA[i*nBasis + j] - PA2[i*nBasis + j])
			*(PA[i*nBasis + j] - PA2[i*nBasis + j]);
		avgdP += (PB[i*nBasis + j] - PB2[i*nBasis + j])
			*(PB[i*nBasis + j] - PB2[i*nBasis + j]);
	}
	avgdP = sqrt(avgdP / nBasis / nBasis / 2.0);

	// compute dPA2 and dPB2
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		dPA2[i*nBasis + j] = PA[i*nBasis + j] - PA2[i*nBasis + j];
		dPB2[i*nBasis + j] = PB[i*nBasis + j] - PB2[i*nBasis + j];
	}

	switch (nIter){

	case 1: // update using simple damping for the first iteration
	case 2: // update using simple damping for the second iteration
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA[i*nBasis + j] = alpha * PA[i*nBasis + j] +
				(1.0 - alpha) * PA2[i*nBasis + j];

			PB[i*nBasis + j] = alpha * PB[i*nBasis + j] +
				(1.0 - alpha) * PB2[i*nBasis + j];
		}
		break;

	default: // using 3-point DIIS for the rest

		// compute diis variables
		B00 = B11 = B22 = B01 = B02 = B12 = 0.0;

		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			DIIS_B(0, 0); DIIS_B(1, 1); DIIS_B(2, 2);
			DIIS_B(0, 1); DIIS_B(0, 2); DIIS_B(1, 2);
		}
		c1 = (B02*B02 + B00*B12 - B01*B02 - B00*B22 + B01*B22 - B02*B12) /
			(B12*B12 - B01*B12 + B02*B11 + B01*B22 - B02*B12 - B11*B22);
		c2 = (B01*B01 - B00*B11 - B01*B02 + B02*B11 + B12*B00 - B12*B01) /
			(B12*B12 - B01*B12 + B02*B11 + B01*B22 - B02*B12 - B11*B22);
		c0 = 1.0 / (1.0 + c1 + c2);
		c1 = c1*c0;
		c2 = c2*c0;

		// use simple damping if DIIS fails
		if (B22>B11){
			c2 = alpha;
			c1 = 1.0 - alpha;
			c0 = 0.0;
		}

		// apply DIIS updates
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA[i*nBasis + j] = c2*(PA[i*nBasis + j]) +
				c1*(PA1[i*nBasis + j] + dPA1[i*nBasis + j]) +
				c0*(PA0[i*nBasis + j] + dPA0[i*nBasis + j]);
			PB[i*nBasis + j] = c2*(PB[i*nBasis + j]) +
				c1*(PB1[i*nBasis + j] + dPB1[i*nBasis + j]) +
				c0*(PB0[i*nBasis + j] + dPB0[i*nBasis + j]);
		}

		break;

	}

	// store previous density matrix
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		dPA0[i*nBasis + j] = dPA1[i*nBasis + j];
		dPB0[i*nBasis + j] = dPB1[i*nBasis + j];

		dPA1[i*nBasis + j] = dPA2[i*nBasis + j];
		dPB1[i*nBasis + j] = dPB2[i*nBasis + j];

		PA0[i*nBasis + j] = PA1[i*nBasis + j];
		PB0[i*nBasis + j] = PB1[i*nBasis + j];

		PA1[i*nBasis + j] = PA2[i*nBasis + j];
		PB1[i*nBasis + j] = PB2[i*nBasis + j];

		PA2[i*nBasis + j] = PA[i*nBasis + j];
		PB2[i*nBasis + j] = PB[i*nBasis + j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}

//
// conv_diis4 : is a the 4-point version of the conv_diis2
//
// Jan, 2013 - Teepanis Chachiyo
//      Initial implementation and testing
//
float conv_diis4(int nBasis,      // number of basis
	float alpha,    // drag coefficient between 0 and 1
	float *PA,      // output alpha density matrix
	float *PB){     // output beta  density matrix

	static float *PA0 = NULL;        // input alpha density matrix
	static float *PB0 = NULL;        // input beta  density matrix
	static float *PA1 = NULL;        // input alpha density matrix
	static float *PB1 = NULL;        // input beta  density matrix
	static float *PA2 = NULL;        // input alpha density matrix
	static float *PB2 = NULL;        // input beta  density matrix
	static float *PA3 = NULL;        // input alpha density matrix
	static float *PB3 = NULL;        // input beta  density matrix
	static float *dPA0 = NULL;       // error alpha density matrix
	static float *dPB0 = NULL;       // error beta  density matrix
	static float *dPA1 = NULL;       // error alpha density matrix
	static float *dPB1 = NULL;       // error beta  density matrix
	static float *dPA2 = NULL;       // error alpha density matrix
	static float *dPB2 = NULL;       // error beta  density matrix
	static float *dPA3 = NULL;       // error alpha density matrix
	static float *dPB3 = NULL;       // error beta  density matrix
	static int nIter = 0;             // number of iteration
	int i, j;                        // generic loop indexes
	float avgdP;                   // average deviation
	float c3, c2, c1, c0;             // diis coefficients
	float B00, B01, B02, B03, B11,
		B12, B13, B22, B23, B33;     // diis variables

	// reset if requested
	if (nBasis == 0 && alpha == 0.0 && PA == NULL && PB == NULL){
		nIter = 0;
		if (PA0 != NULL){ free(PA0); PA0 = NULL; }
		if (PB0 != NULL){ free(PB0); PB0 = NULL; }
		if (PA1 != NULL){ free(PA1); PA1 = NULL; }
		if (PB1 != NULL){ free(PB1); PB1 = NULL; }
		if (PA2 != NULL){ free(PA2); PA2 = NULL; }
		if (PB2 != NULL){ free(PB2); PB2 = NULL; }
		if (PA3 != NULL){ free(PA3); PA3 = NULL; }
		if (PB3 != NULL){ free(PB3); PB3 = NULL; }

		if (dPA0 != NULL){ free(dPA0); dPA0 = NULL; }
		if (dPB0 != NULL){ free(dPB0); dPB0 = NULL; }
		if (dPA1 != NULL){ free(dPA1); dPA1 = NULL; }
		if (dPB1 != NULL){ free(dPB1); dPB1 = NULL; }
		if (dPA2 != NULL){ free(dPA2); dPA2 = NULL; }
		if (dPB2 != NULL){ free(dPB2); dPB2 = NULL; }
		if (dPA3 != NULL){ free(dPA3); dPA3 = NULL; }
		if (dPB3 != NULL){ free(dPB3); dPB3 = NULL; }
		return 0.0;
	}

	// check range of alpha
	if (alpha < 0.0 || alpha > 1.0){
		printf("conv_diis3 : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix
	////////////////////////////////////////////////////////
	if (nIter == 0){

		// check if array is NULL
		if (PA0 != NULL || PB0 != NULL ||
			PA1 != NULL || PB1 != NULL ||
			PA2 != NULL || PB2 != NULL ||
			PA3 != NULL || PB3 != NULL ||
			dPA0 != NULL || dPB0 != NULL ||
			dPA1 != NULL || dPB1 != NULL ||
			dPA2 != NULL || dPB2 != NULL ||
			dPA3 != NULL || dPB3 != NULL){
			printf("conv_diis4 : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		PA0 = calloc(nBasis*nBasis, sizeof(float));
		PB0 = calloc(nBasis*nBasis, sizeof(float));
		PA1 = calloc(nBasis*nBasis, sizeof(float));
		PB1 = calloc(nBasis*nBasis, sizeof(float));
		PA2 = calloc(nBasis*nBasis, sizeof(float));
		PB2 = calloc(nBasis*nBasis, sizeof(float));
		PA3 = calloc(nBasis*nBasis, sizeof(float));
		PB3 = calloc(nBasis*nBasis, sizeof(float));
		dPA0 = calloc(nBasis*nBasis, sizeof(float));
		dPB0 = calloc(nBasis*nBasis, sizeof(float));
		dPA1 = calloc(nBasis*nBasis, sizeof(float));
		dPB1 = calloc(nBasis*nBasis, sizeof(float));
		dPA2 = calloc(nBasis*nBasis, sizeof(float));
		dPB2 = calloc(nBasis*nBasis, sizeof(float));
		dPA3 = calloc(nBasis*nBasis, sizeof(float));
		dPB3 = calloc(nBasis*nBasis, sizeof(float));
		if (PA0 == NULL || PB0 == NULL ||
			PA1 == NULL || PB1 == NULL ||
			PA2 == NULL || PB2 == NULL ||
			PA3 == NULL || PB3 == NULL ||
			dPA0 == NULL || dPB0 == NULL ||
			dPA1 == NULL || dPB1 == NULL ||
			dPA2 == NULL || dPB2 == NULL ||
			dPA3 == NULL || dPB3 == NULL){
			printf("conv_diis4 : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA3[i*nBasis + j] = PA[i*nBasis + j];
			PB3[i*nBasis + j] = PB[i*nBasis + j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	////////////////////////////////////////////////////////
	/// subsequent iteration
	////////////////////////////////////////////////////////

	// compute rms deviation
	avgdP = 0.0;
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		avgdP += (PA[i*nBasis + j] - PA3[i*nBasis + j])
			*(PA[i*nBasis + j] - PA3[i*nBasis + j]);
		avgdP += (PB[i*nBasis + j] - PB3[i*nBasis + j])
			*(PB[i*nBasis + j] - PB3[i*nBasis + j]);
	}
	avgdP = sqrt(avgdP / nBasis / nBasis / 2.0);

	// compute dPA3 and dPB3
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		dPA3[i*nBasis + j] = PA[i*nBasis + j] - PA3[i*nBasis + j];
		dPB3[i*nBasis + j] = PB[i*nBasis + j] - PB3[i*nBasis + j];
	}

	switch (nIter){

	case 1: // update using simple damping for the first iteration
	case 2: // update using simple damping for the second iteration
	case 3: // update using simple damping for the second iteration
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA[i*nBasis + j] = alpha * PA[i*nBasis + j] +
				(1.0 - alpha) * PA3[i*nBasis + j];

			PB[i*nBasis + j] = alpha * PB[i*nBasis + j] +
				(1.0 - alpha) * PB3[i*nBasis + j];
		}
		break;

	default: // using 4-point DIIS for the rest

		// compute diis variables
		B00 = B01 = B02 = B03 = B11 = B12 = B13 = B22 = B23 = B33 = 0.0;

		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			DIIS_B(0, 0); DIIS_B(0, 1); DIIS_B(0, 2); DIIS_B(0, 3);
			DIIS_B(1, 1); DIIS_B(1, 2); DIIS_B(1, 3);
			DIIS_B(2, 2); DIIS_B(2, 3);
			DIIS_B(3, 3);
		}

		float A[25] = { 0.0, 0.0, 0.0, 0.0, -1.0,
			0.0, 0.0, 0.0, 0.0, -1.0,
			0.0, 0.0, 0.0, 0.0, -1.0,
			0.0, 0.0, 0.0, 0.0, -1.0,
			-1.0, -1.0, -1.0, -1.0, 0.0 };
		float B[5] = { 0.0, 0.0, 0.0, 0.0, -1.0 };
		A[0] = B00; A[1] = B01; A[2] = B02; A[3] = B03;
		A[5] = B01; A[6] = B11; A[7] = B12; A[8] = B13;
		A[10] = B02; A[11] = B12; A[12] = B22; A[13] = B23;
		A[15] = B03; A[16] = B13; A[17] = B23; A[18] = B33;
		linear_solver(5, 1, A, B);
		c0 = B[0]; c1 = B[1]; c2 = B[2]; c3 = B[3];


		// use simple damping if DIIS fails
		if (B33>B22){
			c3 = alpha;
			c2 = 1.0 - alpha;
			c1 = 0.0;
			c0 = 0.0;
		}

		// apply DIIS updates
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++){
			PA[i*nBasis + j] = c3*(PA[i*nBasis + j]) +
				c2*(PA2[i*nBasis + j] + dPA2[i*nBasis + j]) +
				c1*(PA1[i*nBasis + j] + dPA1[i*nBasis + j]) +
				c0*(PA0[i*nBasis + j] + dPA0[i*nBasis + j]);
			PB[i*nBasis + j] = c3*(PB[i*nBasis + j]) +
				c2*(PB2[i*nBasis + j] + dPB2[i*nBasis + j]) +
				c1*(PB1[i*nBasis + j] + dPB1[i*nBasis + j]) +
				c0*(PB0[i*nBasis + j] + dPB0[i*nBasis + j]);
		}

		break;

	}

	// store previous density matrix
	for (i = 0; i < nBasis; i++)
	for (j = 0; j < nBasis; j++){
		dPA0[i*nBasis + j] = dPA1[i*nBasis + j];
		dPB0[i*nBasis + j] = dPB1[i*nBasis + j];

		dPA1[i*nBasis + j] = dPA2[i*nBasis + j];
		dPB1[i*nBasis + j] = dPB2[i*nBasis + j];

		dPA2[i*nBasis + j] = dPA3[i*nBasis + j];
		dPB2[i*nBasis + j] = dPB3[i*nBasis + j];

		PA0[i*nBasis + j] = PA1[i*nBasis + j];
		PB0[i*nBasis + j] = PB1[i*nBasis + j];

		PA1[i*nBasis + j] = PA2[i*nBasis + j];
		PB1[i*nBasis + j] = PB2[i*nBasis + j];

		PA2[i*nBasis + j] = PA3[i*nBasis + j];
		PB2[i*nBasis + j] = PB3[i*nBasis + j];

		PA3[i*nBasis + j] = PA[i*nBasis + j];
		PB3[i*nBasis + j] = PB[i*nBasis + j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}
