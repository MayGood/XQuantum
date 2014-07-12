/*
 * fgamma.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <math.h>
#include "fgamma.h"

// fgamma_s : compute Fm(x) using series expansion
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
//		Initial implementation and testing
//
#define CONV 1.0E-15
#define ITMAX 100
float fgamma_s(int m, float x){
	float a;
	float sum = 0.0;
	int i;

	m = m + m + 1;
	a = 1.0 / m;
	x = 2.0*x;
	for (i = 1; i < ITMAX; i++){
		sum += a;
		a = a*x / (m + i + i);
		if (a<CONV) break;
	}
	return exp(-x / 2.0)*sum;
}
#undef CONV
#undef ITMAX

// fgamma_steed : compute Fm(x) using continued fraction
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation and testing
//
#define ITMAX 100
#define CONV  1.0e-15
float fgamma_steed(int m, float x){
	int i;
	float a, bi, ai, D, df, f;

	// compute continued fraction
	// using Steed's method
	a = m + 0.5;
	bi = x + 1.0 - a;
	ai = 0.0;
	D = f = df = 1.0 / bi;
	for (i = 1; i <= ITMAX; i++){
		ai += a - i - i + 1;
		bi += 2.0;
		D = 1.0 / (bi + ai*D);
		df *= (bi*D - 1.0);
		f += df;
		if (fabs(df / f) < CONV) break;
	}
	// compute 1*3*5...(2m-1)/(2x)^m
	D = 1.0;
	a = 1.0 / (x + x);
	for (i = 1; i <= (m + m - 1); i += 2)
		D *= i*a;

#define SQRT_PI_OVER_2 (float)0.8862269254528
	D *= SQRT_PI_OVER_2 / sqrt(x);
#undef  SQRT_PI_OVER_2
	return D - 0.5*exp(-x)*f;
}
#undef ITMAX
#undef CONV

// fgamma_0 :  compute F0(x) using erf and sqrt relation.
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation & testing
//
// May 19, 2011 - Teepanis Chachiyo
//      handle case x approaches zero
float fgamma_0(float x){

	float t, e;

	// case x approaches zero
	if (x < 0.0005)
		return    1.0 - x / 3.0 + x*x / 10.0;

	t = sqrt(x);
	return (float)0.8862269254528*erf(t) / t;
}

// fgamma :  compute Fm(x) using power serie or continued
// fraction depending on the range of x.
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation & testing
//
float fgamma(int m, float x){
	if (x < (m + 0.5))
		return fgamma_s(m, x);
	else
		return fgamma_steed(m, x);
}

// fgamma_set : compute a range of Fm(x) from m=0 upto
// (including) m=max using recurrence relation.
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation  & testing
//
void fgamma_set(int max, float x, float *fset){
	float a, b, c;
	int i, j;

	// case x approaches zero
	if (x < 0.0005){
		for (i = 0; i <= max; i++)
			fset[i] = 1.0 / (i + i + 1.0)
			- x / (i + i + 3.0)
			+ x*x / (i + i + i + i + 10.0);
		return;
	}

	// compute F0(x) first
	fset[0] = fgamma_0(x);

	if (max < 1) return;

	// forward
	a = 1.0 / (x + x);
	c = a*exp(-x);
	for (i = 1; i <= max; i++){
		if ((b = (i + i - 1)*a)>1.0) break;
		fset[i] = b*fset[i - 1] - c;
	}

	if (max < i) return;

	// backward
	fset[max] = fgamma(max, x);
	a = 1.0 / a;
	c = c*a;
	for (j = max - 1; j >= i; j--)
		fset[j] = (a*fset[j + 1] + c) / (j + j + 1);
}
