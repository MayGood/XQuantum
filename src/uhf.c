/*
 * uhf.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "mol.h"
#include "basis.h"
#include "option.h"
#include "int.h"
#include "matrix.h"
#include "conv.h"
#include "lin.h"

// normalizeC : normalized eigen vector
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static void normalizeC(int nBasis, float *S, float *C){
	float sum;
	int i, j, k;

	for (k = 0; k < nBasis; k++){
		sum = 0.0;
		for (i = 0; i < nBasis; i++)
		for (j = 0; j < nBasis; j++)
			sum += C[k*nBasis + i]
			* C[k*nBasis + j]
			* S[i*nBasis + j];
		sum = 1.0 / sqrt(sum);
		for (i = 0; i < nBasis; i++)
			C[k*nBasis + i] = sum * C[k*nBasis + i];
	}
	return;
}

// uhf_getDMatrix : compute density matrix
//
// July 10, 2010 - Teepanis Chachiyo
//    Migrate to uhf scheme
//
// 2008 - Teepanis Chachiyo
//    Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//    Not using static anymore
//
void uhf_getDMatrix(int nBasis, int nOcc, float *C, float *P){
	int i,j,k;
	float sum;

	for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			sum = 0.0;
			for(k=0; k < nOcc; k++){
				sum += C[k*nBasis+i] * C[k*nBasis+j];
			}
			P[i*nBasis+j] = sum;
		}
	return;
}

// getEtotal : compute total energy this is equal to
// electronic part + nuclei part
//
// July 12, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculations
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static float uhf_getEtotal(
	int nBasis, struct Molecule_t *mol,
	float *PA, float *FA,
	float *PB, float *FB,
	float *H){

	float E=0.0;
	int i,j;

	// compute nuclei repulsion energy
//	E += nuclei_coulomb(mol);

	// include electron energy given in (Szabo and Ostlund)
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		E += 0.5*((H[i*nBasis+j]+FA[i*nBasis+j])*PA[i*nBasis+j] +
			      (H[i*nBasis+j]+FB[i*nBasis+j])*PB[i*nBasis+j]);
	}

	return E;
}

extern void GTO_JK_Matrix_CUDA(
		int nBasis, float *P,
		struct GTOBasis_t *gto,
		float *Schwarz, float cutoff, float *G);

extern void GTO_JK_Matrix_CUDA2(
		int nBasis, float *P,
		struct GTOBasis_t *gto,
		float *Schwarz, float cutoff, float *G);

// uhf_getGMatrix : compute G matrix element for both alpha and beta spin
//
// Mar 6, 2013 - Teepanis Chachiyo
//  Use parallel version
//
// Nov 19, 2012 - Teepanis Chachiyo
//  Passing cutoff value as an argument
//
// July 10, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculations
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static void uhf_getGMatrix(
	int nBasis,
	struct GTOBasis_t *gto,
	float *Schwarz,
	float cutoff,
	float *PA, float *PB,
	float *GA, float *GB,
	struct option_t *opt){

	int i, j;
	struct timeval startTime;
	struct timeval endTime;

	gettimeofday(&startTime, NULL);
	// reset to zero
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		GA[i*nBasis+j] = 0.0;
		GB[i*nBasis+j] = 0.0;
	}

//	GTO_JK_Matrix_CUDA(nBasis, PA, gto, Schwarz, cutoff, GA);
//	GTO_JK_Matrix_CUDA2(nBasis, PA, gto, Schwarz, cutoff, GA);
	GTO_JK_Matrix(nBasis, PA, gto, Schwarz, cutoff, GA);
//	GTO_JK_Matrix1(nBasis, PA, gto, Schwarz, cutoff, GA);
	gettimeofday(&endTime, NULL);
	printf("getG cost time: %ld us.\n", 1000000 * (endTime.tv_sec - startTime.tv_sec) + endTime.tv_usec - startTime.tv_usec);
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		GB[i*nBasis+j] = GA[i*nBasis+j];
	}
}

float uhf(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nEA,                 // total number of spin up electrons
	int nEB,                 // total number of spin down electrons
	float *CA,              // returned molecular alpha spin orbital
	float *CB,              // returned molecular beta spin orbital
	float *eA,              // returned eigen values
	float *eB,              // returned eigen values
	struct option_t *opt){   // global option

	float Etot = 0.0;      // total electronic energy
	float dE = 0.0;        // energy change
	float avgdP = 0.0;     // average change in density matrix
	float gamma = 1.0;     // update coefficient
	float sum = 0.0;       // generic summation variable
	float realCutoff;    // schwarz cutoff
	float thisCutoff;    // current cutoff
	float cutoffA = 0.0;   // first stage cutoff
	float cutoffB = 0.0;   // intermediate cutoff
	float cutoffC = 0.0;   // final stage cutoff
	float sumA = 0.0;      // sum of alpha density matrix
	float sumB = 0.0;      // sum of beta density matrix

	// matrix elements
	float *PA = NULL;     // density matrix
	float *GA = NULL;     // EE matrix
	float *dPA = NULL;    // change in density matrix
	float *dGA = NULL;    // change in G matrix
	float *FA = NULL;     // fock matrix
	float *PB = NULL;     // density matrix
	float *GB = NULL;     // EE matrix
	float *dPB = NULL;    // change in density matrix
	float *dGB = NULL;    // change in G matrix
	float *FB = NULL;     // fock matrix
	float *H = NULL;      // h core matrix
	float *T = NULL;      // kinetic matrix
	float *V = NULL;      // nuclei potential matrix
	float *S = NULL;      // overlap matrix
	float *Schwarz = NULL;// Schwarz upper bound matrix

	struct timeval startTime;
	struct timeval endTime;

	int i,j,iter=0;
	int notConverged=1;

	FILE *fd;            // file descriptor for density matrix

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----            SOLVING HARTREE-FOCK EQUATION          -----\n"
	"-------------------------------------------------------------\n"
	);
	fflush(stdout);

	// report
	if(opt->RHF && opt->UHF){
		printf("uhf - error detect both RHF and UHF activated\n");
		exit(-1);
	}

	if(opt->UHF){
		printf("Requested Unrestricted Hartree-Fock calculations\n");
		printf("There are %d alpha spin and %d beta spin electrons\n", nEA, nEB);
	}

	if(opt->RHF){
		printf("Requested Restricted Hartree-Fock calculations\n");
		printf("There are %d electrons in the density matrix\n", nEA+nEB);
		if(nEA!=nEB){
			printf("uhf - error number of electron in each spin is not the same\n");
			exit(-1);
		}
	}
	printMolecule_XYZ(mol, stdout);

#define ALLOCATE(P)                                      \
P = calloc(nBasis*nBasis, sizeof(float));               \
if(P==NULL){                                             \
	printf("rhf: Error - Cannot allocate memory\n"); \
	exit(EXIT_FAILURE);                              \
}

	// memory allocation
	ALLOCATE(PA);
	ALLOCATE(GA);
	ALLOCATE(dPA);
	ALLOCATE(dGA);
	ALLOCATE(FA);
	ALLOCATE(PB);
	ALLOCATE(GB);
	ALLOCATE(dPB);
	ALLOCATE(dGB);
	ALLOCATE(FB);
	ALLOCATE(H);
	ALLOCATE(T);
	ALLOCATE(V);
	ALLOCATE(S);

#undef  ALLOCATE

	/////////////////////////////////////////////
	// Building necessary matrix elements     ///
	/////////////////////////////////////////////

	// report
	printf(
	"Computing 1-electron matrix elements ... \n"
	"    H - Core Hamiltonian                 \n"
	"    S - Overlap Matrix                   \n");

	gettimeofday(&startTime, NULL);
	// get kinetic matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		T[i*nBasis+j] = GTO_kinetic(i,j,gto);
		// symmetrize matrix
		T[j*nBasis+i] = T[i*nBasis+j];
	}

	// get nuclei matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		V[i*nBasis+j] = GTO_nuclei(i,j,gto,mol);
		// symmetrize matrix
		V[j*nBasis+i] = V[i*nBasis+j];
	}

	// get overlap matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		S[i*nBasis+j] = GTO_overlap(i,j,gto);
		// symmetrize matrix
		S[j*nBasis+i] = S[i*nBasis+j];
	}

	// build Hcore
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		H[i*nBasis+j] = T[i*nBasis+j] + V[i*nBasis+j];
	}
	gettimeofday(&endTime, NULL);
	printf("H and S time: %ld us.\n", 1000000 * (endTime.tv_sec - startTime.tv_sec) + endTime.tv_usec - startTime.tv_usec);

	//
	// Diagonalize core hamiltonian to guess density matrix
	//
	if(opt->SCFGuess == SCFGUESS_CORE){
		fflush(stdout);
		printf("Diagonalizing H for initial density matrix ...\n");
		gen_sym_eigen(nBasis, H, S, eA, CA);
		gen_sym_eigen(nBasis, H, S, eB, CB);
		normalizeC(nBasis, S, CA);
		normalizeC(nBasis, S, CB);
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);
	}

	//
	// set density matrix to diagonal as a guess
	//
	if(opt->SCFGuess == SCFGUESS_DIAG){
		fflush(stdout);
		printf("Use identity matrix as initial density matrix ...\n");
		sum=0.0;
		for(i=0; i<nBasis; i++) sum+=S[i*nBasis+i];
		for(i=0; i<nBasis; i++)
		for(j=0; j<nBasis; j++){
			if(i==j){
				PA[i*nBasis+j] = nEA/sum;
				PB[i*nBasis+j] = nEB/sum;
			}else{
				PA[i*nBasis+j] = 0.0;
				PB[i*nBasis+j] = 0.0;
			}
		}
	}


	// compute cutoff
#define GETCUTOFF(accuracy) ((accuracy)/(sumA*sumA+sumA*sumB+sumB*sumB)/50.0)
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		sumA += PA[i*nBasis+j];
		sumB += PB[i*nBasis+j];
	}
	realCutoff = GETCUTOFF(opt->SCFConv);
	opt->SCFCutoff = realCutoff;

	// 2-electron integral information
	printf("Processing 2E integrals ...\n");
	printf("Schwarz inequality screening cut-off %.8E\n", opt->SCFCutoff);
	printf("Primitive prefactor cut-off %0.8E\n", PRIMITIVE_CUTOFF);
	Schwarz = create_Schwarz(nBasis, gto);

	// scf loop
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   printf("Use 4-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DIIS3:   printf("Use 3-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DIIS2:   printf("Use 2-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DAMPING: printf("Use simple weighting convergence method\n"); break;
	default:
		printf("uhf - error no specific request for convergence method\n");
		exit(-1);
	break;
	}

	printf("Drag coefficient %f\n", opt->SCFDrag);
	printf("SCFConv %.2E\n", opt->SCFConv);
	printf("SCFMax %d iterations\n", opt->SCFMax);
	printf("Enter SCF loop ... \n");
	printf(
	"Iteration  Total Energy [Hartrees]  RMSD Density\n"
	"------------------------------------------------\n");
	fflush(stdout);

	// oscillation drag coefficient
	gamma = opt->SCFDrag;

	// call convergence function for the first time to initialize it
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   conv_diis4(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DIIS3:   conv_diis3(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DIIS2:   conv_diis2(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DAMPING: conv_damping(nBasis, gamma, PA, PB); break;
	default:
		printf("uhf - error unknown opt->convMethod\n");
		exit(-1);
	break;
	}

	// preparations
	for(i=0; i<nBasis; i++)
	for(j=0; j<nBasis; j++){

		// set G matrix to zero
		GA[i*nBasis+j] = 0.0;
		GB[i*nBasis+j] = 0.0;

		// set delta density to the initial density
		dPA[i*nBasis+j] = PA[i*nBasis+j];
		dPB[i*nBasis+j] = PB[i*nBasis+j];
	}

	// manage integral accuracy
	switch(opt->SCFAccuracy){
	case SCFACCURACY_1STEP:
		thisCutoff = realCutoff;
	break;
	case SCFACCURACY_3STEP:
		cutoffA    = GETCUTOFF(SCFACCURACY_3STEP_A);
		cutoffB    = GETCUTOFF(SCFACCURACY_3STEP_B);
		cutoffC    = GETCUTOFF(SCFACCURACY_3STEP_C);
		thisCutoff = cutoffA;
	break;
	default:
		printf("uhf - error unknown SCFAccuracy\n");
		exit(-1);
	break;
	}

	do{
		iter++;

		// compute delta G matrix
		if(opt->UHF)
			uhf_getGMatrix(nBasis, gto, Schwarz, thisCutoff, dPA, dPB, dGA, dGB, opt);
		if(opt->RHF)
			uhf_getGMatrix(nBasis, gto, Schwarz, thisCutoff, dPA, dPA, dGA, dGB, opt);

		// updates
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			// update G matrix
			GA[i*nBasis+j] = GA[i*nBasis+j] + dGA[i*nBasis+j];
			GB[i*nBasis+j] = GB[i*nBasis+j] + dGB[i*nBasis+j];

			// update fock matrix
			FA[i*nBasis+j] = H[i*nBasis+j] + GA[i*nBasis+j];
			FB[i*nBasis+j] = H[i*nBasis+j] + GB[i*nBasis+j];

			// saving current density matrix to dPA and dPB
			dPA[i*nBasis+j] = PA[i*nBasis+j];
			dPB[i*nBasis+j] = PB[i*nBasis+j];
		}

		// solve generalized eigen value problem and normalize orbital
		gen_sym_eigen(nBasis, FA, S, eA, CA);
		gen_sym_eigen(nBasis, FB, S, eB, CB);
		normalizeC(nBasis, S, CA);
		normalizeC(nBasis, S, CB);

		// get new P matrix
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);

		// compute energy and energ difference
		dE     = Etot;
		Etot   = uhf_getEtotal(nBasis, mol, PA, FA, PB, FB, H);
		dE     = Etot - dE;

		// update P matrix using convergence method
		switch(opt->convMethod){
		case CONVMETHOD_DIIS4:   avgdP = conv_diis4(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DIIS3:   avgdP = conv_diis3(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DIIS2:   avgdP = conv_diis2(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DAMPING: avgdP = conv_damping(nBasis, gamma, PA, PB); break;
		}

		// check convergence
		notConverged = fabs(dE) > opt->SCFConv || avgdP    > opt->SCFConv;

		// compute delta density matrix for the next step
		for(i=0; i<nBasis; i++)
		for(j=0; j<nBasis; j++){
			dPA[i*nBasis+j] = PA[i*nBasis+j] - dPA[i*nBasis+j];
			dPB[i*nBasis+j] = PB[i*nBasis+j] - dPB[i*nBasis+j];
		}

		printf(" %5d %20.4f %20.4E\n", iter, Etot, avgdP);

		// check if we have reached scfmax limit
		if(iter >= opt->SCFMax) break;

		// flush output
		fflush(stdout);

		// manage integral accuracy
		if(opt->SCFAccuracy == SCFACCURACY_3STEP){

			// check convergence
			if(!notConverged) break;

			// check if we need to switch accuracy
			if((thisCutoff == cutoffA && fabs(dE) <= SCFACCURACY_3STEP_A && avgdP <= SCFACCURACY_3STEP_A) ||
			   (thisCutoff == cutoffB && fabs(dE) <= SCFACCURACY_3STEP_B && avgdP <= SCFACCURACY_3STEP_B) ||
			   (thisCutoff == cutoffC && fabs(dE) <= SCFACCURACY_3STEP_C && avgdP <= SCFACCURACY_3STEP_C)){

				// switch accuracy
				     if(thisCutoff == cutoffA) thisCutoff = cutoffB;
				else if(thisCutoff == cutoffB) thisCutoff = cutoffC;
				else if(thisCutoff == cutoffC) thisCutoff = realCutoff;

				// rebuilding fock matrix
				//for(i=0; i < nBasis; i++)
				//for(j=0; j < nBasis; j++){
				//	 GA[i*nBasis+j] = 0.0;
				//	 GB[i*nBasis+j] = 0.0;
				//	dPA[i*nBasis+j] = PA[i*nBasis+j];
				//	dPB[i*nBasis+j] = PB[i*nBasis+j];
				//}
				printf("................. switch accuracy ..............\n");
				//printf("........ switch accuracy and rebuild fock matrix\n");
				fflush(stdout);

			}

		}

	}while(notConverged);

	// report
	if(notConverged){
		printf("SCF have not converged because iter >= SCFMax\n");
	}else{
		printf("Done SCF Loop Total Energy is %20.8f Hartrees\n", Etot);
	}
	fflush(stdout);

	// clear convergence routine
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   conv_diis4(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DIIS3:   conv_diis3(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DIIS2:   conv_diis2(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DAMPING: conv_damping(0, 0.0, NULL, NULL); break;
	}

	// clean memory
	free(Schwarz);
	free(PA);
	free(GA);
	free(dPA);
	free(dGA);
	free(FA);
	free(PB);
	free(GB);
	free(dPB);
	free(dGB);
	free(FB);
	free(H);
	free(T);
	free(V);
	free(S);

	if(notConverged) return 0.0;
	else             return Etot;
}
