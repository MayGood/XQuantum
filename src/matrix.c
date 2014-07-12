/*
 * matrix.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "matrix.h"
#include "int.h"
#include "quartet.h"
#include "multipole.h"

//
// GTO_moment : computes moment integral between
// two Gaussian basis functions. It has the following form:
//
//                mx     my     mz
// < GTO_i | (x-xc) (y-yc) (z-zc)  | GTO_j >
//
// 2011 - Aniwat Kesorn
//     Adapted from GTO_overlap
//
// Oct 13, 2012 - Teepanis Chachiyo Aniwat Kesorn
//     Add to Siam Quantum
//
float GTO_moment(int i,                            // ith basis
	int j,                            // jth basis
	const struct GTOBasis_t *gto,     // basis database
	int mx, int my, int mz,           // moment order
	float xc, float yc, float zc){ // moment center

	int iCnt, jCnt;   // contracted function
	float sum = 0.0;   // integral sum

	// looping over contracted functions
	for (iCnt = 0; iCnt < gto[i].nContract; iCnt++){
		for (jCnt = 0; jCnt < gto[j].nContract; jCnt++){
			sum = sum + gto[i].coef[iCnt] * gto[j].coef[jCnt] *
				gto[i].norm[iCnt] * gto[j].norm[jCnt] *
				moment(gto[i].exp[iCnt], gto[i].l, gto[i].m, gto[i].n,
				gto[i].x0, gto[i].y0, gto[i].z0,
				gto[j].exp[jCnt], gto[j].l, gto[j].m, gto[j].n,
				gto[j].x0, gto[j].y0, gto[j].z0,
				mx, my, mz, xc, yc, zc);
		}
	}
	return sum;
}

//
// GTO_overlap : computes overlap integral between
// two Gaussian basis functions. It has the following form:
// < GTO_i | GTO_j >
//
// Jan 20, 2008 - Teepanis Chachiyo
//     Iinitial implementation
//
float GTO_overlap(int i,                          // ith basis
	int j,                          // jth basis
	const struct GTOBasis_t *gto){  // basis database

	int iCnt, jCnt;   // contracted function
	float sum = 0.0;   // integral sum

	// looping over contracted functions
	for (iCnt = 0; iCnt < gto[i].nContract; iCnt++){
		for (jCnt = 0; jCnt < gto[j].nContract; jCnt++){
			sum = sum + gto[i].coef[iCnt] * gto[j].coef[jCnt] *
				gto[i].norm[iCnt] * gto[j].norm[jCnt] *
				overlap(gto[i].exp[iCnt], gto[i].l, gto[i].m, gto[i].n,
				gto[i].x0, gto[i].y0, gto[i].z0,
				gto[j].exp[jCnt], gto[j].l, gto[j].m, gto[j].n,
				gto[j].x0, gto[j].y0, gto[j].z0);
		}
	}
	return sum;
}

// GTO_kinetic : computes kinetic operator integral between
// two Gaussian basis functions. It has the following form:
// < GTO_i | -1/2 Laplacian | GTO_j >
//
// Feb 22, 2008 - Teepanis Chachiyo
//     Initial implementation
//
float GTO_kinetic(int i,                          // ith basis
	int j,                          // jth basis
	const struct GTOBasis_t *gto){  // basis database

	int iCnt, jCnt;  // contracted function
	float sum = 0.0;  // integral sum

	// looping over contracted functions
	for (iCnt = 0; iCnt < gto[i].nContract; iCnt++){
		for (jCnt = 0; jCnt < gto[j].nContract; jCnt++){
			sum = sum + gto[i].coef[iCnt] * gto[j].coef[jCnt] *
				gto[i].norm[iCnt] * gto[j].norm[jCnt] *
				kinetic(gto[i].exp[iCnt], gto[i].l, gto[i].m, gto[i].n,
				gto[i].x0, gto[i].y0, gto[i].z0,
				gto[j].exp[jCnt], gto[j].l, gto[j].m, gto[j].n,
				gto[j].x0, gto[j].y0, gto[j].z0);
		}
	}
	return sum;
}

// GTO_nuclei : computes nuclei attraction integral between
// two Gaussian basis functions. It has the following form:
// < GTO_i |   SUM Z_k/(r-R_k) |  GTO_j >
// where index k is a sum over all nuclei
//
// Jan 29, 2008 - Teepanis Chachiyo
//     Initial implementation
//
float GTO_nuclei(int i,                          // ith basis
	int j,                          // jth basis
	const struct GTOBasis_t *gto,   // basis database
	const struct Molecule_t *mol){  // molecule database

	int k;           // nucleus index
	int iCnt, jCnt;  // contracted functions
	float sum = 0.0;  // integral sum

	// looping over contracted functions
	for (iCnt = 0; iCnt < gto[i].nContract; iCnt++){
		for (jCnt = 0; jCnt < gto[j].nContract; jCnt++){

			// looping over nuclei
			for (k = 0; k < mol->nAtom; k++){
				sum = sum + gto[i].coef[iCnt] * gto[j].coef[jCnt] *
					gto[i].norm[iCnt] * gto[j].norm[jCnt] *
					(float)mol->Z[k] *
					nai(gto[i].x0, gto[i].y0, gto[i].z0,
					1.0,
					gto[i].l, gto[i].m, gto[i].n,
					gto[i].exp[iCnt],
					gto[j].x0, gto[j].y0, gto[j].z0,
					1.0,
					gto[j].l, gto[j].m, gto[j].n,
					gto[j].exp[jCnt],
					mol->x[k], mol->y[k], mol->z[k]);
			}
		}
	}
	return sum;

}

// contr_coulomb : compute contracted integral using THO.
// The coding structure is from PyQuante.
//
// Feb 19, 2008 - Teepanis Chachiyo
//     Original implementation.
//
float contr_eri(
	int lena, float *aexps, float *acoefs, float *anorms,
	float xa, float ya, float za, int la, int ma, int na,
	int lenb, float *bexps, float *bcoefs, float *bnorms,
	float xb, float yb, float zb, int lb, int mb, int nb,
	int lenc, float *cexps, float *ccoefs, float *cnorms,
	float xc, float yc, float zc, int lc, int mc, int nc,
	int lend, float *dexps, float *dcoefs, float *dnorms,
	float xd, float yd, float zd, int ld, int md, int nd){

	float val = 0.;
	int i, j, k, l;
	float EE;

	// proceed from highest exponent value
	for (i = 0; i<lena; i++)
	for (j = 0; j<lenb; j++)
	for (k = 0; k<lenc; k++)
	for (l = 0; l<lend; l++){
		// compute element
		EE = acoefs[i] * bcoefs[j] * ccoefs[k] * dcoefs[l]
			* eri(xa, ya, za, anorms[i], la, ma, na, aexps[i],
			xb, yb, zb, bnorms[j], lb, mb, nb, bexps[j],
			xc, yc, zc, cnorms[k], lc, mc, nc, cexps[k],
			xd, yd, zd, dnorms[l], ld, md, nd, dexps[l]);
		val += EE;
	}

	return val;
}

// create_Schwarz : computes the Schwarz matrix and return pointer to the
// matrix. This technique essentially reduces that size from N**4 problem
// to a managable N**2 problem.
//
// Douglas L. Strout and Gustavo E. Scuseria. "A quantitative study of the
// scaling properties of the Hartree-Fock method." J. Chem. Phys. (1995)
// Vol 102 page 8448
//
// Feb 18, 2008 - Teepanis Chachiyo
//     Initial implementation
//
float * create_Schwarz(
	int nBasis,                     // number of basis functions
	const struct GTOBasis_t *gto){  // basis set info

	int p, q;        // loop index
	float *sch;    // Schwarz matrix pointer
	float upBound; // upper bound from Schwarz inequality

	struct timeval startTime;
	struct timeval endTime;

	gettimeofday(&startTime, NULL);

	// allocate memory
	sch = calloc(nBasis*nBasis, sizeof(float));
	if (sch == NULL){
		printf("Cannot allocate Schwarz matrix\n");
		exit(-1);
	}

	// compute the matrix
	// There is a sch[p,q] and sch[q,p] symmetry, so we
	// need to take advantage of this property to gain
	// speed in the future.
	//
	// Feb 22, 2008 - Teepanis Chachiyo
	//
	for (p = 0; p < nBasis; p++){
		for (q = 0; q < nBasis; q++){

			// evaluate integral
			upBound =
				contr_eri(
				gto[p].nContract,
				gto[p].exp, gto[p].coef, gto[p].norm,
				gto[p].x0, gto[p].y0, gto[p].z0,
				gto[p].l, gto[p].m, gto[p].n,
				gto[q].nContract,
				gto[q].exp, gto[q].coef, gto[q].norm,
				gto[q].x0, gto[q].y0, gto[q].z0,
				gto[q].l, gto[q].m, gto[q].n,
				gto[p].nContract,
				gto[p].exp, gto[p].coef, gto[p].norm,
				gto[p].x0, gto[p].y0, gto[p].z0,
				gto[p].l, gto[p].m, gto[p].n,
				gto[q].nContract,
				gto[q].exp, gto[q].coef, gto[q].norm,
				gto[q].x0, gto[q].y0, gto[q].z0,
				gto[q].l, gto[q].m, gto[q].n);

			// make sure we have positive value
			sch[p*nBasis + q] = sqrt(fabs(upBound));
		}
	}
	gettimeofday(&endTime, NULL);
	printf("Schwarz time: %ld us.\n", 1000000 * (endTime.tv_sec - startTime.tv_sec) + endTime.tv_usec - startTime.tv_usec);
	return sch;
}

#define LNBASEI 2.4663034624
#define ERROR_BASE 1.5
#define ERROR_MIN       -120
#define ERROR_MAX       +120
#define ERROR_UNKNOWN   +121
#define ERROR_NOINFO    +122

#define MXZERO 0
#define MYZERO 1
#define MZZERO 2
#define MXMAX  3
#define MYMAX  4
#define MZMAX  5

#define PQIJ_ALLSAME     5
#define PQIJ_2PAIRS      4
#define PQIJ_PQPAIR      3
#define PQIJ_IJPAIR      2
#define PQIJ_PIQJPAIR    1
#define PQIJ_ALLDISTINCT 0

// isSameShell determine that two basis function i and j belong
// to the same shell. A definition of shell is a set of basis
// function with the same center and same contracted exponent,
// for example, sp shell, d shell.
//
// Mar 08, 2010 - Teepanis Chachiyo
//     Original Implementation
//
int isSameShell(
	int nBasis,                    // number of basis functions
	int i, int j,                  // which basis to test
	const struct GTOBasis_t *gto){ // basis set info

	int k;

	// test invalid range
	if (i >= nBasis || j >= nBasis){
		printf("isSameShell - invalid i,j\n");
		exit(-1);
	}

	// test center
	if (gto[i].x0 != gto[j].x0 ||
		gto[i].y0 != gto[j].y0 ||
		gto[i].z0 != gto[j].z0)
		return 0;

	// test contracted exponent
	if (gto[i].nContract != gto[j].nContract) return 0;

	for (k = 0; k < gto[i].nContract; k++)
	if (gto[i].exp[k] != gto[j].exp[k]) return 0;

	return 1;
}

// genGTOShell : generate the shell structure from basis set information.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
// March 8, 2013 - Teepanis Chachiyo
//    Check if the number of contraction does not exceed MAXCONTRACT
//
struct GTOShell_t * genGTOShell(
	int nBasis,                         // number of basis function
	const struct GTOBasis_t *gto,       // basis function data structure
	int *RnShell){                      // returned number of shell created

	int nShell;               // number of shell
	int p, q, i;                // generic counter
	struct GTOShell_t *shell; // pointer to shell data structure


	// check the number of contraction
	for (p = 0; p < nBasis; p++)
	if (gto[p].nContract > MAXCONTRACT){
		printf("genGTOShell: too many contractions, increase MAXCONTRACT\n");
		exit(-1);
	}

	// allocate memory
	shell = calloc(nBasis, sizeof(struct GTOShell_t));
	if (shell == NULL){
		printf("genGTOShell: cannot allocate memory\n");
		exit(-1);
	}

	// filling shell data
	nShell = 0;
	for (p = 0; p < nBasis; p++){

		// copy values
		shell[nShell].min = p;
		shell[nShell].x = gto[p].x0;
		shell[nShell].y = gto[p].y0;
		shell[nShell].z = gto[p].z0;
		shell[nShell].nContract = gto[p].nContract;
		for (i = 0; i < gto[p].nContract; i++)
			shell[nShell].exps[i] = gto[p].exp[i];

		// expand basis to cover this shell
		for (q = p; q < nBasis; q++) {
			if (isSameShell(nBasis, p, q, gto)){

				// determine maxL
				if (shell[nShell].maxL < (gto[q].l + gto[q].m + gto[q].n))
					shell[nShell].maxL = (gto[q].l + gto[q].m + gto[q].n);

				// copy values
				shell[nShell].max = q;
				shell[nShell].l[shell[nShell].nBasis] = gto[q].l;
				shell[nShell].m[shell[nShell].nBasis] = gto[q].m;
				shell[nShell].n[shell[nShell].nBasis] = gto[q].n;
				for (i = 0; i < gto[q].nContract; i++)
					shell[nShell].coef[i][shell[nShell].nBasis] =
					gto[q].norm[i] * gto[q].coef[i];

				shell[nShell].nBasis++;
			}
			else
				break;
		}

		p += (shell[nShell].nBasis - 1);
		nShell++;
	}
	shell = realloc(shell, nShell*sizeof(struct GTOShell_t));

	// determine shell type
	for (p = 0; p < nShell; p++)
		switch (shell[p].nBasis){
		case  1: shell[p].type = TYPES; break;
		case  3: shell[p].type = TYPEP; break;
		case  4: shell[p].type = TYPEL; break;
		case  6: shell[p].type = TYPED; break;
		case 10: shell[p].type = TYPEF; break;
		default: shell[p].type = TYPEX; break;
	}

	*RnShell = nShell;
	return shell;
}

// shellCpy : copies data between shell quartet
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void shellCpy(struct GTOShell_t *dest, struct GTOShell_t *src){
	int i, j;
	dest->x = src->x;
	dest->y = src->y;
	dest->z = src->z;
	dest->nBasis = src->nBasis;
	dest->nContract = src->nContract;
	dest->maxL = src->maxL;
	dest->min = src->min;
	dest->max = src->max;
	dest->type = src->type;

	for (i = 0; i < MAXBASIS; i++){
		dest->l[i] = src->l[i];
		dest->m[i] = src->m[i];
		dest->n[i] = src->n[i];
		dest->exps[i] = src->exps[i];
		for (j = 0; j < MAXBASIS; j++)
			dest->coef[i][j] = src->coef[i][j];
	}
}

// sortGTOShell: sort the GTOShell according their shell types
// in order to enhance "instruction-cache" localization
// when calling the computeQuartet and storeQuartet
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void sortGTOShell(int nShell, struct GTOShell_t *shell){
	struct GTOShell_t buffer;
	int P, Q;

	for (P = 0; P < nShell; P++)
	for (Q = 0; Q < (nShell - 1); Q++){
		if (shell[Q].type > shell[Q + 1].type){
			shellCpy(&buffer, shell + Q);
			shellCpy(shell + Q, shell + Q + 1);
			shellCpy(shell + Q + 1, &buffer);
		}
	}
}

// getErrExpQuartetMPole : compute the log(x) base 1.5 of maximum error within
// this shell quartet when using multipole approximation compared to the
// exact value of the integrals.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
signed char getErrExpQuartetMPole(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	const struct Multipole_t *mp,
	int nBasis,
	float *EEStore){

	int p, q, i, j;                            // basis function index
	const struct Multipole_t *m1, *m2;       // pointer to multipole structure
	const struct Multipole_t *m1ptr, *m2ptr; // starting points
	float K;                               // inverse of distance
	float maxErr = 0.0;

	m1ptr = mp + P->min * nBasis + Q->min;
	m2ptr = mp + I->min * nBasis + J->min;

	// loop thru all basis quartet in the shell
	for (p = 0; p < P->nBasis; p++)
	for (q = 0; q < Q->nBasis; q++)
	for (i = 0; i < I->nBasis; i++)
	for (j = 0; j < J->nBasis; j++){

		m1 = m1ptr + p*nBasis + q;
		m2 = m2ptr + i*nBasis + j;

		K = (m1->x - m2->x)*(m1->x - m2->x) +
			(m1->y - m2->y)*(m1->y - m2->y) +
			(m1->z - m2->z)*(m1->z - m2->z);

		// compute multipole value and its error
		if (K < MULTIPOLE_RADII2_CUTOFF)
			K = *EEStore - 0.0;
		else
			K = *EEStore - m1->q * m2->q * sqrt(1 / K);

		// get maximum
		if (fabs(K) > maxErr) maxErr = fabs(K);

		EEStore++;
	}

	// get log(x) base 1.5 of the maximum
	if (maxErr == 0.0) maxErr = ERROR_MIN; else maxErr = log(maxErr)*LNBASEI;
	if (maxErr < ERROR_MIN) maxErr = ERROR_MIN;
	if (maxErr > ERROR_MAX) maxErr = ERROR_MAX;

	return (signed char)maxErr;
}

// computeQuartetMPole : compute the shell quartet integral using
// multipole expansion
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void computeQuartetMPole(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	const struct Multipole_t *mp,
	int nBasis,
	float *EEStore){

	int p, q, i, j;                            // basis function index
	const struct Multipole_t *m1, *m2;       // pointer to multipole structure
	const struct Multipole_t *m1ptr, *m2ptr; // starting points
	float K;                               // inverse of distance

	m1ptr = mp + P->min * nBasis + Q->min;
	m2ptr = mp + I->min * nBasis + J->min;

	// loop thru all integral in the quartet
	for (p = 0; p < P->nBasis; p++)
	for (q = 0; q < Q->nBasis; q++)
	for (i = 0; i < I->nBasis; i++)
	for (j = 0; j < J->nBasis; j++){

		m1 = m1ptr + p*nBasis + q;
		m2 = m2ptr + i*nBasis + j;

		K = (m1->x - m2->x)*(m1->x - m2->x) +
			(m1->y - m2->y)*(m1->y - m2->y) +
			(m1->z - m2->z)*(m1->z - m2->z);

		// compute only the monopole-monopole term
		if (K < MULTIPOLE_RADII2_CUTOFF)
			*EEStore = 0.0;
		else
			*EEStore = m1->q * m2->q * sqrt(1 / K);

		EEStore++;
	}
}

extern float  genSetBxyzSxyzF_gpu(
	float xa, float ya, float za, int maxa, float alphaa,
	float xb, float yb, float zb, int maxb, float alphab,
	float xc, float yc, float zc, int maxc, float alphac,
	float xd, float yd, float zd, int maxd, float alphad,
	float *Bx, float *By, float *Bz,
	float *Sx, float *Sy, float *Sz,
	float *F, int maxL);

extern void computeEE_gpu(
	float r, int pC, int qC, int iC, int jC,
	int nEE, int pBasis, int qBasis, int iBasis, int jBasis,
	struct GTOShell_t *P, struct GTOShell_t *Q, struct GTOShell_t *I, struct GTOShell_t *J,
	unsigned int *mT, unsigned int *mX, unsigned int *mY, unsigned int *mZ,
	float *Bx, float *By, float *Bz,
	float *Sx, float *Sy, float *Sz,
	unsigned int *iBx, unsigned int *iBy, unsigned int *iBz,
	float *EEStore);

void printEEStore(float *EEStore, int nEE) {
	int i;
	for(i = 0; i < nEE; i++) {
		printf("%10.6f", EEStore[i]);
	}
	printf("\n");
}

//
// the global array needed for computeQuartetEE
//
unsigned int iBx[MAXSHELLINT];   // index for Bx storage
unsigned int iBy[MAXSHELLINT];   // index for By storage
unsigned int iBz[MAXSHELLINT];   // index for Bz storage
unsigned int  mX[MAXSHELLINT];   // maximum in x direction
unsigned int  mY[MAXSHELLINT];   // maximum in y direction
unsigned int  mZ[MAXSHELLINT];   // maximum in z direction
unsigned int  mT[MAXSHELLINT];   // type of mX,mY,mZ for loop optimization

float Bx[4 * (MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)];
float By[4 * (MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)];
float Bz[4 * (MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)];
float Sx[4 * (MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)];
float Sy[4 * (MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)];
float Sz[4 * (MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)*(MAXL + 1)];
float  F[4 * (MAXL + 1)];

// computeQuartetEE : compute the integral exactly using THO method.
// For some cases, the "unrolled" version has been created and will
// be much faster.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void computeQuartetEE(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	float *EEStore){

	float fEEStore[MAXSHELLINT];
	float * __restrict tBx, *__restrict tBy, *__restrict tBz;  // pointer to Bx, By, Bz
	float * __restrict tSx, *__restrict tSy, *__restrict tSz;  // pointer to Sx, Sy, Sz
	int kX, kY, kZ;                // summation loop
	float r;                      // generic float precision
	float rpqi;                   // prefactor*p*q*i coefficients
	int pC, qC, iC, jC;               // contraction index
	int p, q, i, j;                   // basis function index
	int nEE;                       // integral index
	register float EE;            // 2e integral
	register float tSum;          // intermediate variables for kX,kY,kZ loop

	struct timeval startTime;
	struct timeval endTime;
/*
#define SS 0
#define LL 1
#define DD 3
	// call pre "unrolled" case if available
	switch (256 * P->type + 64 * Q->type + 8 * I->type + J->type){

	case 256 * SS + 64 * SS + 8 * SS + SS: computeSSSSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * LL + LL: computeLLLLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * DD + DD: computeDDDDQuartetEE(P, Q, I, J, EEStore); return; break;

		// SL pair
	case 256 * SS + 64 * LL + 8 * LL + LL: computeSLLLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * LL + LL: computeLSLLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * SS + LL: computeLLSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * LL + SS: computeLLLSQuartetEE(P, Q, I, J, EEStore); return; break;

	case 256 * LL + 64 * SS + 8 * SS + SS: computeLSSSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * SS + SS: computeSLSSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * SS + 8 * LL + SS: computeSSLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * SS + 8 * SS + LL: computeSSSLQuartetEE(P, Q, I, J, EEStore); return; break;

	case 256 * SS + 64 * SS + 8 * LL + LL: computeSSLLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * SS + LL: computeSLSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * LL + SS: computeSLLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * SS + LL: computeLSSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * LL + SS: computeLSLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * SS + SS: computeLLSSQuartetEE(P, Q, I, J, EEStore); return; break;

		// DS pair
	case 256 * DD + 64 * SS + 8 * SS + SS: computeDSSSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * DD + 8 * SS + SS: computeSDSSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * SS + 8 * DD + SS: computeSSDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * SS + 8 * SS + DD: computeSSSDQuartetEE(P, Q, I, J, EEStore); return; break;

	case 256 * SS + 64 * DD + 8 * DD + SS: computeSDDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * DD + 8 * SS + DD: computeSDSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * SS + 8 * DD + DD: computeSSDDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * SS + DD: computeDSSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * DD + SS: computeDSDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * SS + SS: computeDDSSQuartetEE(P, Q, I, J, EEStore); return; break;

	case 256 * SS + 64 * DD + 8 * DD + DD: computeSDDDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * DD + DD: computeDSDDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * SS + DD: computeDDSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * DD + SS: computeDDDSQuartetEE(P, Q, I, J, EEStore); return; break;

		// DL pair
	case 256 * DD + 64 * LL + 8 * LL + LL: computeDLLLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * DD + 8 * LL + LL: computeLDLLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * DD + LL: computeLLDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * LL + DD: computeLLLDQuartetEE(P, Q, I, J, EEStore); return; break;

	case 256 * LL + 64 * DD + 8 * DD + DD: computeLDDDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * DD + DD: computeDLDDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * LL + DD: computeDDLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * DD + LL: computeDDDLQuartetEE(P, Q, I, J, EEStore); return; break;

	case 256 * LL + 64 * DD + 8 * DD + LL: computeLDDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * DD + 8 * LL + DD: computeLDLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * DD + DD: computeLLDDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * LL + DD: computeDLLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * DD + LL: computeDLDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * LL + LL: computeDDLLQuartetEE(P, Q, I, J, EEStore); return; break;

		// 2S and DL pair
	case 256 * SS + 64 * SS + 8 * DD + LL: computeSSDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * SS + 8 * LL + DD: computeSSLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * DD + 8 * SS + LL: computeSDSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * SS + DD: computeSLSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * DD + 8 * LL + SS: computeSDLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * DD + SS: computeSLDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * SS + LL: computeDSSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * SS + DD: computeLSSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * LL + SS: computeDSLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * DD + SS: computeLSDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * SS + SS: computeDLSSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * DD + 8 * SS + SS: computeLDSSQuartetEE(P, Q, I, J, EEStore); return; break;

		// 2L and SD pair
	case 256 * LL + 64 * LL + 8 * SS + DD: computeLLSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * LL + 8 * DD + SS: computeLLDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * LL + DD: computeLSLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * DD + 8 * LL + SS: computeLDLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * DD + LL: computeLSDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * DD + 8 * SS + LL: computeLDSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * LL + DD: computeSLLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * LL + SS: computeDLLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * DD + LL: computeSLDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * SS + LL: computeDLSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * DD + 8 * LL + LL: computeSDLLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * LL + LL: computeDSLLQuartetEE(P, Q, I, J, EEStore); return; break;

		// 2D and SL pair
	case 256 * DD + 64 * DD + 8 * SS + LL: computeDDSLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * DD + 8 * LL + SS: computeDDLSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * DD + LL: computeDSDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * DD + SS: computeDLDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * SS + 8 * LL + DD: computeDSLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * DD + 64 * LL + 8 * SS + DD: computeDLSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * DD + 8 * DD + LL: computeSDDLQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * DD + 8 * DD + SS: computeLDDSQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * DD + 8 * LL + DD: computeSDLDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * DD + 8 * SS + DD: computeLDSDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * SS + 64 * LL + 8 * DD + DD: computeSLDDQuartetEE(P, Q, I, J, EEStore); return; break;
	case 256 * LL + 64 * SS + 8 * DD + DD: computeLSDDQuartetEE(P, Q, I, J, EEStore); return; break;
	}
#undef SS
#undef LL
#undef DD
*/
	// index preparation
	nEE = 0;
	for (p = 0; p < P->nBasis; p++)
	for (q = 0; q < Q->nBasis; q++)
	for (i = 0; i < I->nBasis; i++)
	for (j = 0; j < J->nBasis; j++){

		// precompute index for efficientcy
		kZ = (J->maxL + 1);
		kY = kZ*(I->maxL + 1);
		kX = kY*(Q->maxL + 1);
		iBx[nEE] = 4 * (MAXL + 1)*(P->l[p] * kX + Q->l[q] * kY + I->l[i] * kZ + J->l[j]);
		iBy[nEE] = 4 * (MAXL + 1)*(P->m[p] * kX + Q->m[q] * kY + I->m[i] * kZ + J->m[j]);
		iBz[nEE] = 4 * (MAXL + 1)*(P->n[p] * kX + Q->n[q] * kY + I->n[i] * kZ + J->n[j]);
		mX[nEE] = P->l[p] + Q->l[q] + I->l[i] + J->l[j];
		mY[nEE] = P->m[p] + Q->m[q] + I->m[i] + J->m[j];
		mZ[nEE] = P->n[p] + Q->n[q] + I->n[i] + J->n[j];
		if (mX[nEE] == 0)                                  mT[nEE] = MXZERO;
		else if (mY[nEE] == 0)                             mT[nEE] = MYZERO;
		else if (mZ[nEE] == 0)                             mT[nEE] = MZZERO;
		else if (mX[nEE] > mY[nEE] && mX[nEE] > mZ[nEE]) mT[nEE] = MXMAX;
		else if (mY[nEE] > mZ[nEE])                      mT[nEE] = MYMAX;
		else                                            mT[nEE] = MZMAX;

		// reset EEStore to zero
		EEStore[nEE] = 0.0;
		fEEStore[nEE] = 0.0;

		// increase number of EE index
		nEE++;
	}
	//printf("================================================\n");

	for (pC = 0; pC < P->nContract; pC++){
		for (qC = 0; qC < Q->nContract; qC++){
			for (iC = 0; iC < I->nContract; iC++){
				for (jC = 0; jC < J->nContract; jC++){
			/*
					r = genSetBxyzSxyzF_gpu(P->x, P->y, P->z, P->maxL, P->exps[pC],
							Q->x, Q->y, Q->z, Q->maxL, Q->exps[qC],
							I->x, I->y, I->z, I->maxL, I->exps[iC],
							J->x, J->y, J->z, J->maxL, J->exps[jC],
							Bx, By, Bz, Sx, Sy, Sz, F, MAXL);
					/*/
					r = genSetBxyzSxyzFf(P->x, P->y, P->z, P->maxL, P->exps[pC],
						Q->x, Q->y, Q->z, Q->maxL, Q->exps[qC],
						I->x, I->y, I->z, I->maxL, I->exps[iC],
						J->x, J->y, J->z, J->maxL, J->exps[jC],
						Bx, By, Bz, Sx, Sy, Sz, F, MAXL);//*/
					//if(r-rg > 0.000001 || rg-r > 0.000001)
						//printf("%f\n", r-rg);
					if (r < PRIMITIVE_CUTOFF) continue;
/*
					computeEE_gpu(
						r, pC, qC, iC, jC, P->nBasis * Q->nBasis * I->nBasis * J->nBasis, P->nBasis, Q->nBasis, I->nBasis, J->nBasis,
						P, Q, I, J,	mT, mX, mY, mZ,
						Bx, By, Bz, Sx, Sy, Sz, iBx, iBy, iBz, fEEStore);
					for(p = 0; p < MAXSHELLINT; p++) {
						EEStore[p] = fEEStore[p];
					}
					//printEEStore(EEStore, P->nBasis * Q->nBasis * I->nBasis * J->nBasis);
/*/
					// loop all basis within shell
					gettimeofday(&startTime, NULL);
					nEE = 0;
					for (p = 0; p < P->nBasis; p++){
					for (q = 0; q < Q->nBasis; q++){
					for (i = 0; i < I->nBasis; i++){
						rpqi = r * P->coef[pC][p] * Q->coef[qC][q] * I->coef[iC][i];
						for (j = 0; j < J->nBasis; j++){


							// compute two-electron integral using summation
							EE = 0.0;

							switch (mT[nEE]){

							case MXZERO:
								tBy = By + iBy[nEE];
								tSz = Sz + iBz[nEE];
								for (kY = mY[nEE]; kY >= 0; kY--) EE += tBy[kY] * tSz[kY];
								break;

							case MYZERO:
								tBz = Bz + iBz[nEE];
								tSx = Sx + iBx[nEE];
								for (kZ = mZ[nEE]; kZ >= 0; kZ--) EE += tBz[kZ] * tSx[kZ];
								break;

							case MZZERO:
								tBx = Bx + iBx[nEE];
								tSy = Sy + iBy[nEE];
								for (kX = mX[nEE]; kX >= 0; kX--) EE += tBx[kX] * tSy[kX];
								break;

							case MXMAX:
								tBz = Bz + iBz[nEE];
								tSx = Sx + iBx[nEE];
								tBy = By + iBy[nEE];
								for (kY = mY[nEE]; kY >= 0; kY--){
									for (tSum = 0.0, kZ = mZ[nEE]; kZ >= 0; kZ--) tSum += tBz[kZ] * tSx[kY + kZ];
									EE += tBy[kY] * tSum;
								}
								break;

							case MYMAX:
								tBx = Bx + iBx[nEE];
								tSy = Sy + iBy[nEE];
								tBz = Bz + iBz[nEE];
								for (kZ = mZ[nEE]; kZ >= 0; kZ--){
									for (tSum = 0.0, kX = mX[nEE]; kX >= 0; kX--) tSum += tBx[kX] * tSy[kX + kZ];
									EE += tBz[kZ] * tSum;
								}
								break;

							case MZMAX:
								tBy = By + iBy[nEE];
								tSz = Sz + iBz[nEE];
								tBx = Bx + iBx[nEE];
								for (kX = mX[nEE]; kX >= 0; kX--){
									for (tSum = 0.0, kY = mY[nEE]; kY >= 0; kY--) tSum += tBy[kY] * tSz[kX + kY];
									EE += tBx[kX] * tSum;
								}
								break;
							}

							EEStore[nEE] += EE * rpqi * J->coef[jC][j];
							nEE++;
						}
					}
					}
					}
					gettimeofday(&endTime, NULL);
					//printf("CPU cost time: %ld us.\n", 1000000 * (endTime.tv_sec - startTime.tv_sec) + endTime.tv_usec - startTime.tv_usec);
					//printEEStore(EEStore, P->nBasis * Q->nBasis * I->nBasis * J->nBasis);//*/
				}
			}
		}
	}
}

// storeQuartetEE : put the integral in the fock matrix by multiplying with
// the electron densities. Pre "unrolled" cases are also available. This
// subroutine is only for the case "PQIJ_ALLDISTINCT".
//
// Feb 2013 - Teepanis Chachiyo
//  - Initial implementation and testing
//
// 24 Feb, 2013 - Teepanis Chachiyo
//  - Put the code in C preprocessor form and now supports
//    all type of shells
//  - Use different code for restricted and unrestricted case
//
void storeQuartetEE(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	char pqijT,
	int restricted,
	int nBasis,
	float *GA, float *GB,
	const float *PT, const float *PA, const float *PB,
	const float *EEStore){

	int p, q, i, j;      // basis index
	int pp, qq, ii, jj;  // multiplier index
	const float *PTptr, *PAptr, *PBptr;
	register float Ga, Gb, EE;
/*
#define SS 0
#define LL 1
#define DD 3
	if (pqijT == PQIJ_ALLDISTINCT){

		if (restricted)
			switch (256 * P->type + 64 * Q->type + 8 * I->type + J->type){

			case 256 * SS + 64 * SS + 8 * SS + SS: RstoreSSSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * LL + LL: RstoreLLLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * DD + DD: RstoreDDDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

				// SL pair
			case 256 * SS + 64 * LL + 8 * LL + LL: RstoreSLLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * LL + LL: RstoreLSLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * SS + LL: RstoreLLSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * LL + SS: RstoreLLLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

			case 256 * LL + 64 * SS + 8 * SS + SS: RstoreLSSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * SS + SS: RstoreSLSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * LL + SS: RstoreSSLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * SS + LL: RstoreSSSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

			case 256 * SS + 64 * SS + 8 * LL + LL: RstoreSSLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * SS + LL: RstoreSLSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * LL + SS: RstoreSLLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * SS + LL: RstoreLSSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * LL + SS: RstoreLSLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * SS + SS: RstoreLLSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

				// DS pair
			case 256 * DD + 64 * SS + 8 * SS + SS: RstoreDSSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * SS + SS: RstoreSDSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * DD + SS: RstoreSSDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * SS + DD: RstoreSSSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

			case 256 * SS + 64 * DD + 8 * DD + SS: RstoreSDDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * SS + DD: RstoreSDSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * DD + DD: RstoreSSDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * SS + DD: RstoreDSSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * DD + SS: RstoreDSDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * SS + SS: RstoreDDSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

			case 256 * SS + 64 * DD + 8 * DD + DD: RstoreSDDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * DD + DD: RstoreDSDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * SS + DD: RstoreDDSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * DD + SS: RstoreDDDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

				// DL pair
			case 256 * DD + 64 * LL + 8 * LL + LL: RstoreDLLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * LL + LL: RstoreLDLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * DD + LL: RstoreLLDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * LL + DD: RstoreLLLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

			case 256 * LL + 64 * DD + 8 * DD + DD: RstoreLDDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * DD + DD: RstoreDLDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * LL + DD: RstoreDDLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * DD + LL: RstoreDDDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

			case 256 * LL + 64 * DD + 8 * DD + LL: RstoreLDDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * LL + DD: RstoreLDLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * DD + DD: RstoreLLDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * LL + DD: RstoreDLLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * DD + LL: RstoreDLDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * LL + LL: RstoreDDLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

				// 2S and DL pair
			case 256 * SS + 64 * SS + 8 * DD + LL: RstoreSSDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * LL + DD: RstoreSSLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * SS + LL: RstoreSDSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * SS + DD: RstoreSLSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * LL + SS: RstoreSDLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * DD + SS: RstoreSLDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * SS + LL: RstoreDSSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * SS + DD: RstoreLSSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * LL + SS: RstoreDSLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * DD + SS: RstoreLSDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * SS + SS: RstoreDLSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * SS + SS: RstoreLDSSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

				// 2L and SD pair
			case 256 * LL + 64 * LL + 8 * SS + DD: RstoreLLSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * DD + SS: RstoreLLDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * LL + DD: RstoreLSLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * LL + SS: RstoreLDLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * DD + LL: RstoreLSDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * SS + LL: RstoreLDSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * LL + DD: RstoreSLLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * LL + SS: RstoreDLLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * DD + LL: RstoreSLDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * SS + LL: RstoreDLSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * LL + LL: RstoreSDLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * LL + LL: RstoreDSLLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;

				// 2D and SL pair
			case 256 * DD + 64 * DD + 8 * SS + LL: RstoreDDSLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * LL + SS: RstoreDDLSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * DD + LL: RstoreDSDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * DD + SS: RstoreDLDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * LL + DD: RstoreDSLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * SS + DD: RstoreDLSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * DD + LL: RstoreSDDLQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * DD + SS: RstoreLDDSQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * LL + DD: RstoreSDLDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * SS + DD: RstoreLDSDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * DD + DD: RstoreSLDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * DD + DD: RstoreLSDDQuartetEE(P, Q, I, J, nBasis, GA, PT, PA, EEStore); return; break;
		}

		else
			switch (256 * P->type + 64 * Q->type + 8 * I->type + J->type){

			case 256 * SS + 64 * SS + 8 * SS + SS: UstoreSSSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * LL + LL: UstoreLLLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * DD + DD: UstoreDDDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

				// SL pair
			case 256 * SS + 64 * LL + 8 * LL + LL: UstoreSLLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * LL + LL: UstoreLSLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * SS + LL: UstoreLLSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * LL + SS: UstoreLLLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

			case 256 * LL + 64 * SS + 8 * SS + SS: UstoreLSSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * SS + SS: UstoreSLSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * LL + SS: UstoreSSLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * SS + LL: UstoreSSSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

			case 256 * SS + 64 * SS + 8 * LL + LL: UstoreSSLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * SS + LL: UstoreSLSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * LL + SS: UstoreSLLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * SS + LL: UstoreLSSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * LL + SS: UstoreLSLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * SS + SS: UstoreLLSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

				// DS pair
			case 256 * DD + 64 * SS + 8 * SS + SS: UstoreDSSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * SS + SS: UstoreSDSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * DD + SS: UstoreSSDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * SS + DD: UstoreSSSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

			case 256 * SS + 64 * DD + 8 * DD + SS: UstoreSDDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * SS + DD: UstoreSDSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * DD + DD: UstoreSSDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * SS + DD: UstoreDSSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * DD + SS: UstoreDSDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * SS + SS: UstoreDDSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

			case 256 * SS + 64 * DD + 8 * DD + DD: UstoreSDDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * DD + DD: UstoreDSDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * SS + DD: UstoreDDSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * DD + SS: UstoreDDDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

				// DL pair
			case 256 * DD + 64 * LL + 8 * LL + LL: UstoreDLLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * LL + LL: UstoreLDLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * DD + LL: UstoreLLDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * LL + DD: UstoreLLLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

			case 256 * LL + 64 * DD + 8 * DD + DD: UstoreLDDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * DD + DD: UstoreDLDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * LL + DD: UstoreDDLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * DD + LL: UstoreDDDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

			case 256 * LL + 64 * DD + 8 * DD + LL: UstoreLDDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * LL + DD: UstoreLDLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * DD + DD: UstoreLLDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * LL + DD: UstoreDLLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * DD + LL: UstoreDLDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * LL + LL: UstoreDDLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

				// 2S and DL pair
			case 256 * SS + 64 * SS + 8 * DD + LL: UstoreSSDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * SS + 8 * LL + DD: UstoreSSLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * SS + LL: UstoreSDSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * SS + DD: UstoreSLSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * LL + SS: UstoreSDLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * DD + SS: UstoreSLDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * SS + LL: UstoreDSSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * SS + DD: UstoreLSSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * LL + SS: UstoreDSLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * DD + SS: UstoreLSDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * SS + SS: UstoreDLSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * SS + SS: UstoreLDSSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

				// 2L and SD pair
			case 256 * LL + 64 * LL + 8 * SS + DD: UstoreLLSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * LL + 8 * DD + SS: UstoreLLDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * LL + DD: UstoreLSLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * LL + SS: UstoreLDLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * DD + LL: UstoreLSDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * SS + LL: UstoreLDSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * LL + DD: UstoreSLLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * LL + SS: UstoreDLLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * DD + LL: UstoreSLDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * SS + LL: UstoreDLSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * LL + LL: UstoreSDLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * LL + LL: UstoreDSLLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;

				// 2D and SL pair
			case 256 * DD + 64 * DD + 8 * SS + LL: UstoreDDSLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * DD + 8 * LL + SS: UstoreDDLSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * DD + LL: UstoreDSDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * DD + SS: UstoreDLDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * SS + 8 * LL + DD: UstoreDSLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * DD + 64 * LL + 8 * SS + DD: UstoreDLSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * DD + LL: UstoreSDDLQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * DD + SS: UstoreLDDSQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * DD + 8 * LL + DD: UstoreSDLDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * DD + 8 * SS + DD: UstoreLDSDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * SS + 64 * LL + 8 * DD + DD: UstoreSLDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
			case 256 * LL + 64 * SS + 8 * DD + DD: UstoreLSDDQuartetEE(P, Q, I, J, nBasis, GA, GB, PT, PA, PB, EEStore); return; break;
		}


	}
#undef SS
#undef LL
#undef DD
*/

	jj = 1;
	ii = J->nBasis;
	qq = I->nBasis * J->nBasis;
	pp = Q->nBasis * I->nBasis * J->nBasis;

	const float *EEptr;

	//
	// the restricted version deals with only GA
	//
#define RSTOREJ_DIDI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
	for (p = 0; p < P->nBasis; p++) 									\
	for (q = 0; q < Q->nBasis; q++){ 									\
	PTptr = PT + I->min * nBasis + J->min; 								\
	EEptr = EEStore + p*pp + q*qq; 										\
	Ga = 0.0; 															\
	for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 		\
	for (j = 0; j < J->nBasis; j++){ 									\
		EE = EEptr[j*jj]; EE += EE; 									\
		Ga += PTptr[j] * EE; 											\
	} 																	\
	GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 						\
	GA[(q + Q->min)*nBasis + (p + P->min)] += Ga; 						\
	}

#define RSTOREJ_DISA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
		for (p = 0; p < P->nBasis; p++) 								\
		for (q = 0; q < Q->nBasis; q++){ 								\
	PTptr = PT + I->min * nBasis + J->min; 								\
	EEptr = EEStore + p*pp + q*qq; 										\
	Ga = 0.0; 															\
		for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 	\
		for (j = 0; j < J->nBasis; j++){ 								\
			EE = EEptr[j*jj]; 											\
			Ga += PTptr[j] * EE; 										\
		} 																\
		GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 					\
		GA[(q + Q->min)*nBasis + (p + P->min)] += Ga; 					\
		}

#define RSTOREJ_SADI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (q = 0; q < Q->nBasis; q++){ 							\
			PTptr = PT + I->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + q*qq; 								\
			Ga = 0.0; 													\
			for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; EE += EE; 							\
				Ga += PTptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 				\
			}

#define RSTOREJ_SASA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (q = 0; q < Q->nBasis; q++){ 							\
			PTptr = PT + I->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + q*qq; 								\
			Ga = 0.0; 													\
			for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; 										\
				Ga += PTptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 				\
			}

#define RSTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (i = 0; i < I->nBasis; i++){ 							\
			PAptr = PA + Q->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + i*ii; 								\
			Ga = 0.0; 													\
			for (q = 0; q < Q->nBasis; q++, EEptr += qq, PAptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; 										\
				Ga += PAptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (i + I->min)] -= Ga; 				\
			GA[(i + I->min)*nBasis + (p + P->min)] -= Ga; 				\
			}

#define RSTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (i = 0; i < I->nBasis; i++){ 							\
			PAptr = PA + Q->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + i*ii; 								\
			Ga = 0.0; 													\
			for (q = 0; q < Q->nBasis; q++, EEptr += qq, PAptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; 										\
				Ga += PAptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (i + I->min)] -= Ga; 				\
			}

			//
			// the unrestricted version deals with both GA and GB
			//
#define USTOREJ_DIDI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (q = 0; q < Q->nBasis; q++){ 							\
			PTptr = PT + I->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + q*qq; 								\
			Ga = 0.0; Gb = 0.0; 										\
			for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; EE += EE; 							\
				Ga += PTptr[j] * EE; 									\
				Gb += PTptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 				\
			GB[(p + P->min)*nBasis + (q + Q->min)] += Gb; 				\
			GA[(q + Q->min)*nBasis + (p + P->min)] += Ga; 				\
			GB[(q + Q->min)*nBasis + (p + P->min)] += Gb; 				\
			}

#define USTOREJ_DISA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (q = 0; q < Q->nBasis; q++){ 							\
			PTptr = PT + I->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + q*qq; 								\
			Ga = 0.0; Gb = 0.0; 										\
			for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; 										\
				Ga += PTptr[j] * EE; 									\
				Gb += PTptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 				\
			GB[(p + P->min)*nBasis + (q + Q->min)] += Gb; 				\
			GA[(q + Q->min)*nBasis + (p + P->min)] += Ga; 				\
			GB[(q + Q->min)*nBasis + (p + P->min)] += Gb; 				\
			}

#define USTOREJ_SADI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (q = 0; q < Q->nBasis; q++){ 							\
			PTptr = PT + I->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + q*qq; 								\
			Ga = 0.0; Gb = 0.0; 										\
			for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; EE += EE; 							\
				Ga += PTptr[j] * EE; 									\
				Gb += PTptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 				\
			GB[(p + P->min)*nBasis + (q + Q->min)] += Gb; 				\
			}

#define USTOREJ_SASA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (q = 0; q < Q->nBasis; q++){ 							\
			PTptr = PT + I->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + q*qq; 								\
			Ga = 0.0; Gb = 0.0; 										\
			for (i = 0; i < I->nBasis; i++, EEptr += ii, PTptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; 										\
				Ga += PTptr[j] * EE; 									\
				Gb += PTptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (q + Q->min)] += Ga; 				\
			GB[(p + P->min)*nBasis + (q + Q->min)] += Gb; 				\
			}

#define USTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (i = 0; i < I->nBasis; i++){ 							\
			PAptr = PA + Q->min * nBasis + J->min; 						\
			PBptr = PB + Q->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + i*ii; 								\
			Ga = 0.0; Gb = 0.0; 										\
			for (q = 0; q < Q->nBasis; q++, EEptr += qq, PAptr += nBasis, PBptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; 										\
				Ga += PAptr[j] * EE; 									\
				Gb += PBptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (i + I->min)] -= Ga; 				\
			GB[(p + P->min)*nBasis + (i + I->min)] -= Gb; 				\
			GA[(i + I->min)*nBasis + (p + P->min)] -= Ga; 				\
			GB[(i + I->min)*nBasis + (p + P->min)] -= Gb; 				\
			}

#define USTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj) 						\
			for (p = 0; p < P->nBasis; p++) 							\
			for (i = 0; i < I->nBasis; i++){ 							\
			PAptr = PA + Q->min * nBasis + J->min; 						\
			PBptr = PB + Q->min * nBasis + J->min; 						\
			EEptr = EEStore + p*pp + i*ii; 								\
			Ga = 0.0; Gb = 0.0; 										\
			for (q = 0; q < Q->nBasis; q++, EEptr += qq, PAptr += nBasis, PBptr += nBasis) 	\
			for (j = 0; j < J->nBasis; j++){ 							\
				EE = EEptr[j*jj]; 										\
				Ga += PAptr[j] * EE; 									\
				Gb += PBptr[j] * EE; 									\
			} 															\
			GA[(p + P->min)*nBasis + (i + I->min)] -= Ga; 				\
			GB[(p + P->min)*nBasis + (i + I->min)] -= Gb; 				\
			}

			if (restricted)
				switch (pqijT){
				case PQIJ_ALLSAME:
					RSTOREJ_SASA(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREK_SAME(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					break;

				case PQIJ_2PAIRS:
					RSTOREJ_SASA(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREJ_SASA(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					RSTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					break;

				case PQIJ_PQPAIR:
					RSTOREJ_SADI(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREJ_DISA(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					RSTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREK_DIFF(P, Q, J, I, p, q, j, i, pp, qq, jj, ii);
					break;

				case PQIJ_IJPAIR:
					RSTOREJ_DISA(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREJ_SADI(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					RSTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREK_DIFF(Q, P, I, J, q, p, i, j, qq, pp, ii, jj);
					break;

				case PQIJ_PIQJPAIR:
					RSTOREJ_DIDI(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREK_SAME(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREK_SAME(P, Q, J, I, p, q, j, i, pp, qq, jj, ii);
					RSTOREK_SAME(Q, P, I, J, q, p, i, j, qq, pp, ii, jj);
					RSTOREK_SAME(Q, P, J, I, q, p, j, i, qq, pp, jj, ii);
					break;

				case PQIJ_ALLDISTINCT:
					RSTOREJ_DIDI(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREJ_DIDI(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					RSTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					RSTOREK_DIFF(P, Q, J, I, p, q, j, i, pp, qq, jj, ii);
					RSTOREK_DIFF(Q, P, I, J, q, p, i, j, qq, pp, ii, jj);
					RSTOREK_DIFF(Q, P, J, I, q, p, j, i, qq, pp, jj, ii);
					break;
			}

			else
				switch (pqijT){
				case PQIJ_ALLSAME:
					USTOREJ_SASA(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREK_SAME(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					break;

				case PQIJ_2PAIRS:
					USTOREJ_SASA(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREJ_SASA(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					USTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					break;

				case PQIJ_PQPAIR:
					USTOREJ_SADI(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREJ_DISA(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					USTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREK_DIFF(P, Q, J, I, p, q, j, i, pp, qq, jj, ii);
					break;

				case PQIJ_IJPAIR:
					USTOREJ_DISA(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREJ_SADI(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					USTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREK_DIFF(Q, P, I, J, q, p, i, j, qq, pp, ii, jj);
					break;

				case PQIJ_PIQJPAIR:
					USTOREJ_DIDI(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREK_SAME(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREK_SAME(P, Q, J, I, p, q, j, i, pp, qq, jj, ii);
					USTOREK_SAME(Q, P, I, J, q, p, i, j, qq, pp, ii, jj);
					USTOREK_SAME(Q, P, J, I, q, p, j, i, qq, pp, jj, ii);
					break;

				case PQIJ_ALLDISTINCT:
					USTOREJ_DIDI(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREJ_DIDI(I, J, P, Q, i, j, p, q, ii, jj, pp, qq);
					USTOREK_DIFF(P, Q, I, J, p, q, i, j, pp, qq, ii, jj);
					USTOREK_DIFF(P, Q, J, I, p, q, j, i, pp, qq, jj, ii);
					USTOREK_DIFF(Q, P, I, J, q, p, i, j, qq, pp, ii, jj);
					USTOREK_DIFF(Q, P, J, I, q, p, j, i, qq, pp, jj, ii);
					break;
			}
}

// getMemCutoff : get cutoff for memeory of the schwarz matrix so that
// the amount of memory need for storage does not exceed user's defined
// maximum value.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
float getMemCutoff(
	int nShell,
	const struct GTOShell_t *shell,
	float fixedCutoff,
	const float *schwarz_shell,
	const struct option_t *opt){

	int P,Q,I,J,N;                    // shell index
	float cutmem;                    // cutoff value for available memory
	unsigned long maxEE,nEEBank;      // maximum possible number of integrals
	float r;                         // generic float precision

	// increase cutmem until it is enough for available memory
	cutmem = fixedCutoff;
	maxEE  = (unsigned long)(opt->maxMem)*1000000/sizeof(float);

	// set corse search flag
	r = 1;
	do{
CutMemSearch:
		nEEBank = 0;

		// loop all distinct shell sequence
		for(P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I;
		for(J=0; J <= N; J++){


			// screen integral using Schwarz inequality at shell level
			if(schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J] < cutmem)
				continue;

			// accumulate the number of shell permutation
			nEEBank += shell[P].nBasis * shell[Q].nBasis *
			           shell[I].nBasis * shell[J].nBasis;
			if(nEEBank > maxEE){
				if(r==1) // corse search
					cutmem = cutmem*2.00;
				else     // fine search
					cutmem = cutmem*1.05;
				goto CutMemSearch;
			}

		}
		}

		// done corse search and begin fine search
		if(r==1 && nEEBank < maxEE){ r = 0; cutmem/=2.0; goto CutMemSearch; }

	}while(nEEBank > maxEE);

	// cutmem should not be smaller than fixedCutoff
	if(cutmem < fixedCutoff) cutmem = fixedCutoff;

	return cutmem;
}

// GTO_JK_Matrix_Quartet : calculate G matrix for unrestricted hartree-fock
// calculations. The equations are given in (Szabo and Ostlund, 1989)
//
// for Alpha Spin
// [GB]pq = SUM [PT]ij*(pq|ij) - [PA]ij*(pi|qj)
//
// for Beta Spin
// [GA]pg = SUM [PT]ij*(pq|ij) - [PB]ij*(pi|qj)
//
// where PT is the toal spin density and is defined by
// PT = PA+PB
//
// As of Feb 2013, this is the faster version of the code the compute
// fock matrix.
//
// Feb 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
// Feb 24, 2013 - Teepanis Chachiyo
//     No longer call storeJK and use storeQuartetEE instead
//
void GTO_JK_Matrix_Quartet(
	int nBasis,                    // number of basis functions
	const float *PA,              // density matrix for spin up
	const float *PB,              // density matrix for spin down
	const struct GTOBasis_t *gto,  // basis set info
	const float *schwarz_basis,   // pointer to schwarz matrix
	float fixedCutoff,            // cutoff to ignore
	float *GA,                    // return G for spin up
	float *GB,                    // return G for spin down
	struct option_t *opt){         // global option

	static int firstVisit=1;              // first visit flag
	static struct GTOShell_t *shell=NULL; // shell data structure
	static int nShell;                    // the number of shell
	int P,Q,I,J,N;                        // shell loop index
	int p,q;                              // basis loop index
	static float *schwarz_shell=NULL;    // schwarz matrix at shell level
	static float *PM=NULL;               // maximum density at shell level
	static float *PT=NULL;               // total density matrix
	float r;                             // generic float precision
	float maxPShell;                     // maximum density in the shell
	float EEStore[MAXSHELLINT];          // electron integral storage
	char pqijT;                           // type of shell quartet
	int nEE;                              // number of integral within shell

	static struct Multipole_t *mp=NULL;   // multipole data
	static signed char *erBank=NULL;      // multipole error bank at shell level
	signed char *erBankPtr=NULL;          // pointer to current value
	signed char erExp;                    // current value of error

	static float *eeBank=NULL;           // global 2e storage
	float *eeBankPtr=NULL;               // pointer to current 2e storage
	static float memCutoff=0.0;          // schwartz cutoff to store in memory

#define ALLOC(array,item,type)                                        \
array=calloc(item,sizeof(type));                                      \
if(array==NULL){                                                      \
	printf("GTO_JK_Matrix_Quartet - error cannot allocate memory\n"); \
	exit(-1);                                                         \
}

	// reset call
	if(nBasis==0){
		firstVisit = 1;
		if(shell != NULL){ free(shell); shell=NULL; }
		if(schwarz_shell != NULL){ free(schwarz_shell); schwarz_shell=NULL; }
		if(PM != NULL){ free(PM); PM=NULL; }
		if(PT != NULL){ free(PT); PT=NULL; }
		if(mp != NULL){ free(mp); mp=NULL; }
		if(erBank != NULL){ free(erBank); erBank=NULL; }
		if(eeBank != NULL){ free(eeBank); eeBank=NULL; }
		return;
	}

	// first visit preparations
	if(firstVisit){

		// generate shell
		shell = genGTOShell(nBasis, gto, &nShell);
		sortGTOShell(nShell, shell);

		// generate multipole data
		mp = genMultipole(nBasis, gto);

		// allocate memory
		ALLOC(schwarz_shell, nShell*nShell, float);
		ALLOC(PM, nShell*nShell, float);
		ALLOC(PT, nBasis*nBasis, float);

		// compute schwarz at shell level
		for(P=0; P < nShell; P++)
		for(Q=0; Q < nShell; Q++){
			r = 0.0;
			for(p=shell[P].min; p <= shell[P].max; p++)
			for(q=shell[Q].min; q <= shell[Q].max; q++)
				if(r <schwarz_basis[p*nBasis+q]) r = schwarz_basis[p*nBasis+q];
			schwarz_shell[P*nShell+Q] = r;
		}

		// count the number shells needed and allocate memory
		unsigned long nEE=0;
		for(P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I;
		for(J=0; J <= N; J++){
			r = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
			if(r < opt->SCFCutoff) continue;
			nEE++;
		}
		}
		ALLOC(erBank, nEE, signed char);
		memset(erBank, ERROR_UNKNOWN, nEE);

		// determine memory cutoff, number of integrals, and allocate memory
		memCutoff = getMemCutoff(nShell, shell, fixedCutoff, schwarz_shell, opt);
		nEE=0;
		for(P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I;
		for(J=0; J <= N; J++){
			if(schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J] < memCutoff)
				continue;

			nEE += shell[P].nBasis * shell[Q].nBasis *
			       shell[I].nBasis * shell[J].nBasis;
		}
		}
		ALLOC(eeBank, nEE, float);

//		setupGPU(nShell, nBasis, shell);
	}

	// compute total density matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		PT[p*nBasis+q] = PA[p*nBasis+q] + PB[p*nBasis+q];

	// compute schwarz and maximum density at shell level
	for(P=0; P < nShell; P++)
	for(Q=0; Q < nShell; Q++){
		r = 0.0;
		for(p=shell[P].min; p <= shell[P].max; p++)
		for(q=shell[Q].min; q <= shell[Q].max; q++){
			if(r<fabs(PT[p*nBasis+q])) r = fabs(PT[p*nBasis+q]);
			if(r<fabs(PA[p*nBasis+q])) r = fabs(PA[p*nBasis+q]);
			if(r<fabs(PB[p*nBasis+q])) r = fabs(PB[p*nBasis+q]);
		}
		PM[P*nShell+Q] = r;
	}

	// main loop
	erBankPtr = erBank;
	eeBankPtr = eeBank;
	/*
	//setupGPU(nShell, nBasis);
	computeQuartetEE_gpu(
		shell, schwarz_shell, nShell,
		opt->RHF, nBasis, GA, GB,
		PT, PA, PB);
	//resetGPU();//*/
	//*
	for(P=0; P < nShell; P++)
	for(Q=0; Q <= P; Q++)
	for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I;
	for(J=0; J <= N; J++){

		// schwarz screening
		r = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
		if(r < opt->SCFCutoff) continue;
		erExp = *erBankPtr; erBankPtr++;

		if(r < fixedCutoff) continue;

		// maximum density in this shell
		maxPShell = PM[P*nShell+Q] + PM[I*nShell+J] +
		            PM[P*nShell+I] + PM[P*nShell+J] +
		            PM[Q*nShell+I] + PM[Q*nShell+J];

		//computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);

		//*
		// load from memory if possible
		if(r >= memCutoff){
			nEE = shell[P].nBasis * shell[Q].nBasis *
			      shell[I].nBasis * shell[J].nBasis;
			if(firstVisit){
				computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);
				memcpy(eeBankPtr, EEStore, sizeof(float) * nEE);
				eeBankPtr += nEE;
			}else{
				// screening again with density weighted
				if(r*maxPShell < fixedCutoff){ eeBankPtr += nEE; continue; }

				memcpy(EEStore, eeBankPtr, sizeof(float) * nEE);
				eeBankPtr += nEE;
			}
		}else{
			// screening again with density weighted
			if(r*maxPShell < fixedCutoff) continue;

			// compute the entire shell quartet
			if(erExp!=ERROR_UNKNOWN && pow(ERROR_BASE,erExp)*maxPShell < fixedCutoff)
				computeQuartetMPole(shell+P,shell+Q,shell+I,shell+J,mp,nBasis,EEStore);
			else
				computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);
		}
///
		// store multipole error
		if(erExp==ERROR_UNKNOWN)
			erBankPtr[-1]=getErrExpQuartetMPole(shell+P,shell+Q,shell+I,shell+J,
			                                    mp,nBasis,EEStore);

		// determine P Q I J type for storage
		if((P==Q)&&(I==J)&&(P==I)){ pqijT = PQIJ_ALLSAME;
		}else if((P==Q)&&(I==J)){   pqijT = PQIJ_2PAIRS;
		}else if(P==Q){             pqijT = PQIJ_PQPAIR;
		}else if(I==J){             pqijT = PQIJ_IJPAIR;
		}else if((P==I)&&(Q==J)){   pqijT = PQIJ_PIQJPAIR;
		}else{                      pqijT = PQIJ_ALLDISTINCT;
		}

		// store in the G matrix
		storeQuartetEE(shell+P,shell+Q,shell+I,shell+J,
		               pqijT,opt->RHF,nBasis,GA,GB,PT,PA,PB,EEStore);

	}
	}//*/

	// handle restricted case
	if(opt->RHF) for(p=0; p < nBasis*nBasis; p++) GB[p] = GA[p];

	// symmetrize G matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < p; q++){
		GA[q*nBasis+p] = GA[p*nBasis+q];
		GB[q*nBasis+p] = GB[p*nBasis+q];
	}

	// reset flag
	if(firstVisit) firstVisit=0;
}

// getMemCutoff_Parallel: this is  the parallel version of the getMemCutoff
//
// Mar 3, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
float getMemCutoff_Parallel(
	int childID,
	int nCPU,
	int nShell,
	const struct GTOShell_t *shell,
	float fixedCutoff,
	const float *schwarz_shell,
	const struct option_t *opt){

	int PQ, P, Q, I, J, N;                 // shell index
	float cutmem;                    // cutoff value for available memory
	unsigned long maxEE, nEEBank;      // maximum possible number of integrals
	float r;                         // generic float precision

	// increase cutmem until it is enough for available memory
	cutmem = fixedCutoff;
	maxEE = (unsigned long)(opt->maxMem) * 1000000 / sizeof(float);

	// set corse search flag
	r = 1;
	do{
	CutMemSearch:
		nEEBank = 0;

		// loop all distinct shell sequence
		for (PQ = 0, P = 0; P < nShell; P++)
		for (Q = 0; Q <= P; Q++, PQ++) if (PQ%nCPU == childID)
		for (I = 0; I <= P; I++){
			if (I == P) N = Q; else N = I;
			for (J = 0; J <= N; J++){


				// screen integral using Schwarz inequality at shell level
				if (schwarz_shell[P*nShell + Q] * schwarz_shell[I*nShell + J] < cutmem)
					continue;

				// accumulate the number of shell permutation
				nEEBank += shell[P].nBasis * shell[Q].nBasis *
					shell[I].nBasis * shell[J].nBasis;
				if (nEEBank > maxEE){
					if (r == 1) // corse search
						cutmem = cutmem*2.00;
					else     // fine search
						cutmem = cutmem*1.05;
					goto CutMemSearch;
				}

			}
		}

		// done corse search and begin fine search
		if (r == 1 && nEEBank < maxEE){ r = 0; cutmem /= 2.0; goto CutMemSearch; }

	} while (nEEBank > maxEE);

	// cutmem should not be smaller than fixedCutoff
	if (cutmem < fixedCutoff) cutmem = fixedCutoff;

	return cutmem;
}

// GTO_JK_Matrix_Quartet_Parallel: this is a parallel version of the
// GTO_JK_Matrix_Quartet. See the original subroutine for information.
//
// Mar 3, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void GTO_JK_Matrix_Quartet_Parallel(
	int childID,                   // child id number
	int nBasis,                    // number of basis functions
	const float *PA,              // density matrix for spin up
	const float *PB,              // density matrix for spin down
	const struct GTOBasis_t *gto,  // basis set info
	const float *schwarz_basis,   // pointer to schwarz matrix
	float fixedCutoff,            // cutoff to ignore
	float *GA,                    // return G for spin up
	float *GB,                    // return G for spin down
	struct option_t *opt){         // global option

	static int firstVisit = 1;              // first visit flag
	static struct GTOShell_t *shell = NULL; // shell data structure
	static int nShell;                    // the number of shell
	int PQ, P, Q, I, J, N;                     // shell loop index
	int p, q;                              // basis loop index
	static float *schwarz_shell = NULL;    // schwarz matrix at shell level
	static float *PM = NULL;               // maximum density at shell level
	static float *PT = NULL;               // total density matrix
	float r;                             // generic float precision
	float maxPShell;                     // maximum density in the shell
	float EEStore[MAXSHELLINT];          // electron integral storage
	char pqijT;                           // type of shell quartet
	int nEE;                              // number of integral within shell

	static struct Multipole_t *mp = NULL;   // multipole data
	static signed char *erBank = NULL;      // multipole error bank at shell level
	signed char *erBankPtr = NULL;          // pointer to current value
	signed char erExp;                    // current value of error

	static float *eeBank = NULL;           // global 2e storage
	float *eeBankPtr = NULL;               // pointer to current 2e storage
	static float memCutoff = 0.0;          // schwartz cutoff to store in memory

	struct timeval startTime;
	struct timeval curTime;

#define ALLOC(array,item,type)                                        \
	array = calloc(item, sizeof(type));                               \
	if (array == NULL){	                                               \
	printf("GTO_JK_Matrix_Quartet - error cannot allocate memory\n");  \
	exit(-1);                                                          \
	}

		// reset call
	if (nBasis == 0){
		firstVisit = 1;
		if (shell != NULL) {
			free(shell);
			shell = NULL;
		}
		if (schwarz_shell != NULL) {
			free(schwarz_shell);
			schwarz_shell = NULL;
		}
		if (PM != NULL) {
			free(PM);
			PM = NULL;
		}
		if (PT != NULL) {
			free(PT);
			PT = NULL;
		}
		if (mp != NULL) {
			free(mp);
			mp = NULL;
		}
		if (erBank != NULL) {
			free(erBank);
			erBank = NULL;
		}
		if (eeBank != NULL) {
			free(eeBank);
			eeBank = NULL;
		}
		return;
	}
	gettimeofday(&startTime, NULL);

	// first visit preparations
	if (firstVisit){

		// generate shell
		shell = genGTOShell(nBasis, gto, &nShell);
		sortGTOShell(nShell, shell);

		// generate multipole data
		mp = genMultipole(nBasis, gto);

		// allocate memory
		ALLOC(schwarz_shell, nShell*nShell, float);
		ALLOC(PM, nShell*nShell, float);
		ALLOC(PT, nBasis*nBasis, float);

		// compute schwarz at shell level
		for (P = 0; P < nShell; P++) {
			for (Q = 0; Q < nShell; Q++) {
				r = 0.0;
				for (p = shell[P].min; p <= shell[P].max; p++) {
					for (q = shell[Q].min; q <= shell[Q].max; q++) {
						if (r <schwarz_basis[p*nBasis + q]) {
							r = schwarz_basis[p*nBasis + q];
						}
					}
				}
				schwarz_shell[P*nShell + Q] = r;
			}
		}

		// count the number shells needed and allocate memory
		unsigned long nEE = 0;
		for (PQ = 0, P = 0; P < nShell; P++) {
			for (Q = 0; Q <= P; Q++, PQ++) {
				if (PQ % (opt->nCPU) == childID) {
					for (I = 0; I <= P; I++){
						if (I == P)
							N = Q;
						else
							N = I;
						for (J = 0; J <= N; J++) {
							r = schwarz_shell[P*nShell + Q] * schwarz_shell[I*nShell + J];
							if (r < opt->SCFCutoff)
								continue;
							nEE++;
						}
					}
				}
			}
		}
		ALLOC(erBank, nEE, signed char);
		memset(erBank, ERROR_UNKNOWN, nEE);

		// determine memory cutoff, number of integrals, and allocate memory
		memCutoff = getMemCutoff_Parallel(childID, opt->nCPU, nShell, shell,
			fixedCutoff, schwarz_shell, opt);
		nEE = 0;
		for (PQ = 0, P = 0; P < nShell; P++) {
			for (Q = 0; Q <= P; Q++, PQ++) {
				if (PQ % (opt->nCPU) == childID) {
					for (I = 0; I <= P; I++) {
						if (I == P)
							N = Q;
						else
							N = I;
						for (J = 0; J <= N; J++) {
							if (schwarz_shell[P*nShell + Q] * schwarz_shell[I*nShell + J] < memCutoff)
								continue;

							nEE += shell[P].nBasis * shell[Q].nBasis *
								shell[I].nBasis * shell[J].nBasis;
						}
					}
				}
			}
		}
		ALLOC(eeBank, nEE, float);

	}

	gettimeofday(&curTime, NULL);
	printf("first visit preparations time: %ld us.\n", 1000000 * (curTime.tv_sec - startTime.tv_sec) + curTime.tv_usec - startTime.tv_usec);

	// compute total density matrix
	for (p = 0; p < nBasis; p++) {
		for (q = 0; q < nBasis; q++) {
			PT[p*nBasis + q] = PA[p*nBasis + q] + PB[p*nBasis + q];
		}
	}

	// compute schwarz and maximum density at shell level
	for (P = 0; P < nShell; P++) {
		for (Q = 0; Q < nShell; Q++) {
			r = 0.0;
			for (p = shell[P].min; p <= shell[P].max; p++) {
				for (q = shell[Q].min; q <= shell[Q].max; q++) {
					if (r<fabs(PT[p*nBasis + q]))
						r = fabs(PT[p*nBasis + q]);
					if (r<fabs(PA[p*nBasis + q]))
						r = fabs(PA[p*nBasis + q]);
					if (r<fabs(PB[p*nBasis + q]))
						r = fabs(PB[p*nBasis + q]);
				}
			}
			PM[P*nShell + Q] = r;
		}
	}

	gettimeofday(&startTime, NULL);
	// main loop
	erBankPtr = erBank;
	eeBankPtr = eeBank;
	for (PQ = 0, P = 0; P < nShell; P++) {
		for (Q = 0; Q <= P; Q++, PQ++) {
			if (PQ % (opt->nCPU) == childID) {
				for (I = 0; I <= P; I++) {
					if (I == P)
						N = Q;
					else
						N = I;
					for (J = 0; J <= N; J++) {

						// schwarz screening
						r = schwarz_shell[P*nShell + Q] * schwarz_shell[I*nShell + J];
						if (r < opt->SCFCutoff)
							continue;
						erExp = *erBankPtr;
						erBankPtr++;

						if (r < fixedCutoff)
							continue;

						// maximum density in this shell
						maxPShell = PM[P*nShell + Q] + PM[I*nShell + J] +
							PM[P*nShell + I] + PM[P*nShell + J] +
							PM[Q*nShell + I] + PM[Q*nShell + J];

						// load from memory if possible
						if (r >= memCutoff) {
							nEE = shell[P].nBasis * shell[Q].nBasis *
								shell[I].nBasis * shell[J].nBasis;
							if (firstVisit) {
								computeQuartetEE(shell + P, shell + Q, shell + I, shell + J, EEStore);
								memcpy(eeBankPtr, EEStore, sizeof(float)* nEE);
								eeBankPtr += nEE;
							}
							else {
								// screening again with density weighted
								if (r*maxPShell < fixedCutoff) {
									eeBankPtr += nEE;
									continue;
								}

								memcpy(EEStore, eeBankPtr, sizeof(float)* nEE);
								eeBankPtr += nEE;
							}
						}
						else{
							// screening again with density weighted
							if (r*maxPShell < fixedCutoff)
								continue;

							// compute the entire shell quartet
							if (erExp != ERROR_UNKNOWN && pow(ERROR_BASE, erExp)*maxPShell < fixedCutoff) {
								computeQuartetMPole(shell + P, shell + Q, shell + I, shell + J, mp, nBasis, EEStore);
							}
							else {
								computeQuartetEE(shell + P, shell + Q, shell + I, shell + J, EEStore);
							}
						}

						// store multipole error
						if (erExp == ERROR_UNKNOWN) {
							erBankPtr[-1] = getErrExpQuartetMPole(shell + P, shell + Q, shell + I, shell + J,
							mp, nBasis, EEStore);
						}

						// determine P Q I J type for storage
						if ((P == Q) && (I == J) && (P == I)){
							pqijT = PQIJ_ALLSAME;
						}
						else if ((P == Q) && (I == J)){
							pqijT = PQIJ_2PAIRS;
						}
						else if (P == Q){
							pqijT = PQIJ_PQPAIR;
						}
						else if (I == J){
							pqijT = PQIJ_IJPAIR;
						}
						else if ((P == I) && (Q == J)){
							pqijT = PQIJ_PIQJPAIR;
						}
						else{
							pqijT = PQIJ_ALLDISTINCT;
						}

						// store in the G matrix
						storeQuartetEE(shell + P, shell + Q, shell + I, shell + J,
							pqijT, opt->RHF, nBasis, GA, GB, PT, PA, PB, EEStore);

					}
				}
			}
		}
	}
	gettimeofday(&curTime, NULL);
	printf("main loop time: %ld us.\n", 1000000 * (curTime.tv_sec - startTime.tv_sec) + curTime.tv_usec - startTime.tv_usec);

	// handle restricted case
	if (opt->RHF) {
		for (p = 0; p < nBasis*nBasis; p++) {
			GB[p] = GA[p];
		}
	}

	// symmetrize G matrix
	for (p = 0; p < nBasis; p++) {
		for (q = 0; q < p; q++) {
			GA[q*nBasis + p] = GA[p*nBasis + q];
			GB[q*nBasis + p] = GB[p*nBasis + q];
		}
	}

	// reset flag
	if (firstVisit)
		firstVisit = 0;
}
#undef ALLOC

extern void GTO_JK_Matrix_CUDA1(
		unsigned int nEE, unsigned int *iter_p, unsigned int *iter_q, unsigned int *iter_i, unsigned int *iter_j,
		int nBasis, float *P,
		struct GTOBasis_t *gto, float *G);
// GTO_JK_Matrix : compute 2J-K matrix elements taking into account of
// the symmetry in the Coulomb integrals. This type of matrix is for
// closed-shell calculations.
//
// It computes repulsion integral between molecular orbital
// described by density matrix P and a pair of basis i,j. It has the following
// form: [RETURN]pq = SUM Pij*(2 EE(piqj) - EE(pijq))
//
// Note: For closed-shell system: Gpq = [RETURN]pq / 2.0
//
// Important Definition
// --------------------
//
// (piqj) = Int[  dr1 dr2 phi_p*(r1) phi_i*(r2) (1/r12) phi_q(r1) phi_j(r2) ]
//
// Jan 29, 2008 - Teepanis Chachiyo
//     Original Implementation
//
// Feb 18, 2008 - Teepanis Chachiyo
//     Optimization on G[p,q] multiplications
//
void GTO_JK_Matrix(
	int nBasis,                  // number of basis functions
	float *P,                   // density matrix
	struct GTOBasis_t *gto,      // basis set info
	float *Schwarz,             // pointer to schwarz matrix
	float cutoff,               // cutoff to ignore
	float *G){                  // return matrix

	int p,q,i,j,maxj;       // looping variables
	float EE;              // coulomb integral
	float upBound;         // Schwarz upper bound

	// algorithm based on my algorithm during graduate years
	// which rooted from Thijssen.
	// This is for looping all unique (piqj)
	// HF.tar.gz  -- Teepanis Chachiyo
	for(p=0; p < nBasis; p++){
		for(q=0; q < p+1; q++){
			for(i=0; i < p+1; i++){
				if(i == p){
					maxj = q+1;
				}else{
					maxj = i+1;
				}
				for(j=0; j < maxj; j++){
#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
					// screen integral using Schwarz inequality
					upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
					if(upBound < SCHWARZ_CUTOFF) continue;
					upBound = sqrt(upBound);
					if(fabs(upBound*P(p,q)) < cutoff &&
					   fabs(upBound*P(p,i)) < cutoff &&
					   fabs(upBound*P(p,j)) < cutoff &&
					   fabs(upBound*P(q,i)) < cutoff &&
					   fabs(upBound*P(q,j)) < cutoff &&
					   fabs(upBound*P(i,j)) < cutoff) continue;
					// compute two-electron integral
					EE =
					contr_eri(
					     gto[p].nContract,
					     gto[p].exp, gto[p].coef, gto[p].norm,
					     gto[p].x0,        gto[p].y0,  gto[p].z0,
					     gto[p].l,         gto[p].m,   gto[p].n,
					     gto[q].nContract,
					     gto[q].exp, gto[q].coef, gto[q].norm,
					     gto[q].x0,        gto[q].y0,  gto[q].z0,
					     gto[q].l,         gto[q].m,   gto[q].n,
					     gto[i].nContract,
					     gto[i].exp, gto[i].coef, gto[i].norm,
					     gto[i].x0,        gto[i].y0,  gto[i].z0,
					     gto[i].l,         gto[i].m,   gto[i].n,
					     gto[j].nContract,
					     gto[j].exp, gto[j].coef, gto[j].norm,
					     gto[j].x0,        gto[j].y0,  gto[j].z0,
					     gto[j].l,         gto[j].m,   gto[j].n);

					// perform matrix symmetry rotations
					if((p==q)&&(i==j)&&(p==i)){  // all same
						G(p,q) = G(p,q) +     EE*P(i,j);
					}else if((p==q)&&(i==j)){    // 2 pairs
						G(p,p) = G(p,p) + 2.0*EE*P(i,i);
						G(i,i) = G(i,i) + 2.0*EE*P(p,p);

						G(p,i) = G(p,i) -     EE*P(i,p);
						G(i,p) = G(i,p) -     EE*P(i,p);
					}else if(p==q){              // pq pair
						G(p,p) = G(p,p) + 4.0*EE*P(i,j);

						G(i,j) = G(i,j) + 2.0*EE*P(p,p);
						G(j,i) = G(j,i) + 2.0*EE*P(p,p);

						G(p,j) = G(p,j) -     EE*P(i,p);
						G(j,p) = G(j,p) -     EE*P(i,p);

						G(i,p) = G(i,p) -     EE*P(p,j);
						G(p,i) = G(p,i) -     EE*P(p,j);
					}else if(i==j){              // ij pair
						G(i,i) = G(i,i) + 4.0*EE*P(p,q);

						G(p,q) = G(p,q) + 2.0*EE*P(i,i);
						G(q,p) = G(q,p) + 2.0*EE*P(i,i);

						G(p,i) = G(p,i) -     EE*P(i,q);
						G(i,p) = G(i,p) -     EE*P(i,q);

						G(q,i) = G(q,i) -     EE*P(i,p);
						G(i,q) = G(i,q) -     EE*P(i,p);
					}else if((p==i)&&(q==j)){    // pi-qj pair
						G(p,q) = G(p,q) + 3.0*EE*P(p,q);
						G(q,p) = G(q,p) + 3.0*EE*P(p,q);

						G(p,p) = G(p,p) -     EE*P(q,q);
						G(q,q) = G(q,q) -     EE*P(p,p);
					}else{                       // all distinct
						G(p,q) = G(p,q) + 4.0*EE*P(i,j);
						G(q,p) = G(q,p) + 4.0*EE*P(i,j);

						G(i,j) = G(i,j) + 4.0*EE*P(p,q);
						G(j,i) = G(j,i) + 4.0*EE*P(p,q);

						G(p,j) = G(p,j) -     EE*P(i,q);
						G(j,p) = G(j,p) -     EE*P(i,q);

						G(p,i) = G(p,i) -     EE*P(j,q);
						G(i,p) = G(i,p) -     EE*P(j,q);

						G(q,i) = G(q,i) -     EE*P(j,p);
						G(i,q) = G(i,q) -     EE*P(j,p);

						G(q,j) = G(q,j) -     EE*P(i,p);
						G(j,q) = G(j,q) -     EE*P(i,p);
					}
				}
			}
		}
	}
#undef G
#undef P
}

#define MAXINT 2147483647
void GTO_JK_Matrix1(
	int nBasis,                  // number of basis functions
	float *P,                   // density matrix
	struct GTOBasis_t *gto,      // basis set info
	float *Schwarz,             // pointer to schwarz matrix
	float cutoff,               // cutoff to ignore
	float *G){                  // return matrix

	int p,q,i,j,maxj;       // looping variables
	float EE;              // coulomb integral
	float upBound;         // Schwarz upper bound

	int nEE = 0;
	unsigned int *iter_p, *iter_q, *iter_i, *iter_j;
#define ALLOC(array,item,type)                                        \
array=calloc(item,sizeof(type));                                      \
if(array==NULL){                                                      \
	printf("GTO_JK_Matrix_Quartet - error cannot allocate memory\n"); \
	exit(-1);                                                         \
}
	ALLOC(iter_p, MAXINT, unsigned int);
	ALLOC(iter_q, MAXINT, unsigned int);
	ALLOC(iter_i, MAXINT, unsigned int);
	ALLOC(iter_j, MAXINT, unsigned int);
	nEE = 0;

	// algorithm based on my algorithm during graduate years
	// which rooted from Thijssen.
	// This is for looping all unique (piqj)
	// HF.tar.gz  -- Teepanis Chachiyo
	for(p=0; p < nBasis; p++){
		for(q=0; q < p+1; q++){
			for(i=0; i < p+1; i++){
				if(i == p){
					maxj = q+1;
				}else{
					maxj = i+1;
				}
				for(j=0; j < maxj; j++){
#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
					// screen integral using Schwarz inequality
					upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
					if(upBound < SCHWARZ_CUTOFF) continue;
					upBound = sqrt(upBound);
					if(fabs(upBound*P(p,q)) < cutoff &&
					   fabs(upBound*P(p,i)) < cutoff &&
					   fabs(upBound*P(p,j)) < cutoff &&
					   fabs(upBound*P(q,i)) < cutoff &&
					   fabs(upBound*P(q,j)) < cutoff &&
					   fabs(upBound*P(i,j)) < cutoff) continue;
					iter_p[nEE] = p;
					iter_q[nEE] = q;
					iter_i[nEE] = i;
					iter_j[nEE] = j;
					nEE++;
				}
			}
		}
	}
	realloc(iter_p, sizeof(unsigned int)*nEE);
	realloc(iter_q, sizeof(unsigned int)*nEE);
	realloc(iter_i, sizeof(unsigned int)*nEE);
	realloc(iter_j, sizeof(unsigned int)*nEE);

	GTO_JK_Matrix_CUDA1(nEE, iter_p, iter_q, iter_i, iter_j, nBasis, P, gto, G);

#undef G
#undef P

#undef ALLOC
}

#undef MAXINT
