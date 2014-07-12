/*
 * quartet.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef QUARTET_H_
#define QUARTET_H_

#define MAXL          4
#define MAXSHELLINT  50625

#define MAXBASIS    16
#define MAXCONTRACT 16
struct GTOShell_t{
	float x, y, z;                             // shell center
	int nBasis;                               // number of basis in the shell
	int nContract;                            // number contract functions
	int maxL;                                 // maximum angular index
	int min, max;                              // basis index
	int l[MAXBASIS], m[MAXBASIS], n[MAXBASIS];  // angular index
	float exps[MAXCONTRACT];                 // exponents
	float coef[MAXCONTRACT][MAXBASIS];       // coefficients
	int type;                                 // predefined type
#define TYPES 0
#define TYPEL 1
#define TYPEP 2
#define TYPED 3
#define TYPEF 4
#define TYPEG 5
#define TYPEX 6
};

#define PROTOTYPE(a,b,c,d)										\
	void compute##a##b##c##d##QuartetEE(						\
	const struct GTOShell_t *P,									\
	const struct GTOShell_t *Q,									\
	const struct GTOShell_t *I,									\
	const struct GTOShell_t *J,									\
	float *EEStore);											\
																\
	void Rstore##a##b##c##d##QuartetEE(							\
	const struct GTOShell_t *P,									\
	const struct GTOShell_t *Q,									\
	const struct GTOShell_t *I,									\
	const struct GTOShell_t *J,									\
	int nBasis,													\
	float *GA,													\
	const float *PT, const float *PA,							\
	const float *EEStore);										\
																\
	void Ustore##a##b##c##d##QuartetEE(							\
	const struct GTOShell_t *P,									\
	const struct GTOShell_t *Q,									\
	const struct GTOShell_t *I,									\
	const struct GTOShell_t *J,									\
	int nBasis,													\
	float *GA, float *GB,										\
	const float *PT, const float *PA, const float *PB,		\
	const float *EEStore);


PROTOTYPE(S, S, S, S);
PROTOTYPE(L, L, L, L);
PROTOTYPE(D, D, D, D);

// LS pair
PROTOTYPE(L, S, S, S);
PROTOTYPE(S, L, S, S);
PROTOTYPE(S, S, L, S);
PROTOTYPE(S, S, S, L);

PROTOTYPE(S, L, L, L);
PROTOTYPE(L, S, L, L);
PROTOTYPE(L, L, S, L);
PROTOTYPE(L, L, L, S);

PROTOTYPE(S, S, L, L);
PROTOTYPE(S, L, S, L);
PROTOTYPE(S, L, L, S);
PROTOTYPE(L, S, S, L);
PROTOTYPE(L, S, L, S);
PROTOTYPE(L, L, S, S);

// SD pair
PROTOTYPE(D, S, S, S);
PROTOTYPE(S, D, S, S);
PROTOTYPE(S, S, D, S);
PROTOTYPE(S, S, S, D);

PROTOTYPE(S, D, D, D);
PROTOTYPE(D, S, D, D);
PROTOTYPE(D, D, S, D);
PROTOTYPE(D, D, D, S);

PROTOTYPE(S, S, D, D);
PROTOTYPE(S, D, S, D);
PROTOTYPE(S, D, D, S);
PROTOTYPE(D, S, S, D);
PROTOTYPE(D, S, D, S);
PROTOTYPE(D, D, S, S);

// LD pair
PROTOTYPE(D, L, L, L);
PROTOTYPE(L, D, L, L);
PROTOTYPE(L, L, D, L);
PROTOTYPE(L, L, L, D);

PROTOTYPE(L, D, D, D);
PROTOTYPE(D, L, D, D);
PROTOTYPE(D, D, L, D);
PROTOTYPE(D, D, D, L);

PROTOTYPE(L, L, D, D);
PROTOTYPE(L, D, L, D);
PROTOTYPE(L, D, D, L);
PROTOTYPE(D, L, L, D);
PROTOTYPE(D, L, D, L);
PROTOTYPE(D, D, L, L);

// 2S and LD pair
PROTOTYPE(S, S, L, D);
PROTOTYPE(S, S, D, L);
PROTOTYPE(S, L, S, D);
PROTOTYPE(S, D, S, L);
PROTOTYPE(S, L, D, S);
PROTOTYPE(S, D, L, S);
PROTOTYPE(L, S, S, D);
PROTOTYPE(D, S, S, L);
PROTOTYPE(L, S, D, S);
PROTOTYPE(D, S, L, S);
PROTOTYPE(L, D, S, S);
PROTOTYPE(D, L, S, S);

// 2L and SD pair
PROTOTYPE(L, L, S, D);
PROTOTYPE(L, L, D, S);
PROTOTYPE(L, S, L, D);
PROTOTYPE(L, D, L, S);
PROTOTYPE(L, S, D, L);
PROTOTYPE(L, D, S, L);
PROTOTYPE(S, L, L, D);
PROTOTYPE(D, L, L, S);
PROTOTYPE(S, L, D, L);
PROTOTYPE(D, L, S, L);
PROTOTYPE(S, D, L, L);
PROTOTYPE(D, S, L, L);

// 2D and SL pair
PROTOTYPE(D, D, S, L);
PROTOTYPE(D, D, L, S);
PROTOTYPE(D, S, D, L);
PROTOTYPE(D, L, D, S);
PROTOTYPE(D, S, L, D);
PROTOTYPE(D, L, S, D);
PROTOTYPE(S, D, D, L);
PROTOTYPE(L, D, D, S);
PROTOTYPE(S, D, L, D);
PROTOTYPE(L, D, S, D);
PROTOTYPE(S, L, D, D);
PROTOTYPE(L, S, D, D);
#undef  PROTOTYPE

#endif /* QUARTET_H_ */
