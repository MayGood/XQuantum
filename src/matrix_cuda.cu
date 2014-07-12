/*
 * matrix_cuda.cu
 *
 *  Created on: Jun 7, 2014
 *      Author: xxhx
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "basis.h"
#include "matrix.h"

/**
 * This macro checks return value of the CUDA runtime call and exits
 * the application if the call failed.
 */
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }

__device__ float atomicFloatAdd(float* address, float val)
{
    int* address_as_i = (int*)address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_i, assumed,
        		__float_as_int(val +
        				__int_as_float(assumed)));
    } while (assumed != old);
    return __int_as_float(old);
}

__device__ float atomicFloatSub(float* address, float val)
{
    int* address_as_i = (int*)address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_i, assumed,
        		__float_as_int(__int_as_float(assumed) - val));
    } while (assumed != old);
    return __int_as_float(old);
}

/*
 * fgamma part
 */
#define CONV 1.0E-15
#define ITMAX 100
__device__ float fgamma_s_device(int m, float x){
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

#define ITMAX 100
#define CONV  1.0e-15
__device__ float fgamma_steed_device(int m, float x){
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

__device__ float fgamma_0_device(float x){

	float t, e;

	// case x approaches zero
	if (x < 0.0005)
		return    1.0 - x / 3.0 + x*x / 10.0;

	t = sqrt(x);
	return (float)0.8862269254528*erf(t) / t;
}

__device__ float fgamma_device(int m, float x){
	if (x < (m + 0.5))
		return fgamma_s_device(m, x);
	else
		return fgamma_steed_device(m, x);
}

__device__ void fgamma_set_device(int max, float x, float *fset){
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
	fset[0] = fgamma_0_device(x);

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
	fset[max] = fgamma_device(max, x);
	a = 1.0 / a;
	c = c*a;
	for (j = max - 1; j >= i; j--)
		fset[j] = (a*fset[j + 1] + c) / (j + j + 1);
}
/* end of fgamma */
#define MAXL 8
__device__ void eriB_1d_device(float *B, int l1, int l2, int l3, int l4,
	float Px, float Ax, float Bx,
	float Qx, float Cx, float Dx,
	float eta1, float eta2,
	float zeta){

	float PA, PB, QC, QD, QP;   // auxilary

	// special cases
#define ERI1111(Ax,Bx,Cx,Dx,Px,Qx,eta1,eta2,zeta) 							\
	PA = Px - Ax; 															\
	PB = Px - Bx; 															\
	QC = Qx - Cx; 															\
	QD = Qx - Dx; 															\
	QP = Qx - Px; 															\
	eta1 *= 0.5; 															\
	eta2 *= 0.5; 															\
	zeta *= 0.5; 															\
	B[0] = (eta1 + PA*PB)*(eta2 + QC*QD); 									\
	B[1] = eta1*(PA + PB)*(QP*(eta2 + QC*QD) + eta2*(QC + QD)); 			\
	B[1] -= eta1*eta1*(eta2 + QC*QD); 										\
	B[1] -= eta2*(eta1 + PA*PB)*(eta2 + QP*(QC + QD)); 						\
	B[1] *= zeta; 															\
	B[2] = eta2*eta2*(eta1 + PA*PB) + eta1*eta1*(eta2 + QC*QD); 			\
	B[2] -= eta1*eta2*(PA + PB)*(QC + QD); 									\
	B[2] *= QP*QP; 															\
	B[2] += 3.0*QP*eta1*eta2*(eta1*(QC + QD) - eta2*(PA + PB)); 			\
	B[2] += 3.0*eta1*eta1*eta2*eta2; 										\
	B[2] *= zeta*zeta; 														\
	B[3] = QP*eta2*(PA + PB) - QP*eta1*(QC + QD) - 6.0*eta1*eta2; 			\
	B[3] *= QP*QP*zeta*zeta*zeta*eta1*eta2; 								\
	B[4] = QP*QP*eta1*eta2*zeta*zeta; 										\
	B[4] = B[4] * B[4];

#define ERI2010(Ax,Cx,Px,Qx,eta1,eta2,zeta) 								\
	PA = Px - Ax; 															\
	QC = Qx - Cx; 															\
	QP = Qx - Px; 															\
	eta1 *= 0.25; 															\
	eta2 *= 0.25; 															\
	B[0] = QC*(PA*PA + eta1 + eta1); 										\
	B[1] = 2.0*eta1*eta2*(PA + PA - QP) + 2.0*eta1*QC*(PA*QP - eta1); 		\
	B[1] -= PA*PA*QP*eta2; 													\
	B[1] *= zeta; 															\
	B[2] = eta1*eta2*6.0 - 2.0*PA*QP*eta2 + QC*QP*eta1; 					\
	B[2] *= QP*eta1*zeta*zeta; 												\
	QP *= QP*QP; 															\
	zeta *= zeta*zeta; 														\
	B[3] = -QP*zeta*eta1*eta1*eta2;

#define ERI2100(Ax,Bx,Px,Qx,eta1,zeta) 										\
	PA = Px - Ax; 															\
	PB = Px - Bx; 															\
	QP = Qx - Px; 															\
	eta1 *= 0.25; 															\
	B[0] = 2.0*eta1*(PA + PA + PB) + PA*PA*PB; 								\
	B[1] = eta1*(6.0*QP - 4.0*PA - PB - PB) + PA*QP*(PA + PB + PB); 		\
	B[1] *= zeta*eta1; 														\
	B[2] = QP*(PA + PA + PB) - 6.0*eta1; 									\
	B[2] *= QP*zeta*zeta*eta1*eta1; 										\
	B[3] = QP*eta1*zeta; 													\
	B[3] *= B[3] * B[3];

#define ERI0111(Bx,Cx,Dx,Px,Qx,eta1,eta2,zeta) 								\
	PB = Px - Bx; 															\
	QC = Qx - Cx; 															\
	QD = Qx - Dx; 															\
	QP = Qx - Px; 															\
	zeta *= 0.25; 															\
	eta1 *= 0.50; 															\
	eta2 *= 0.50; 															\
	B[0] = PB*(eta2 + QC*QD); 												\
	B[1] = eta1*eta2*(QC + QD + QP); 										\
	B[1] += eta1*QC*QD*QP; 													\
	B[1] -= eta2*PB*QP*(QC + QD); 											\
	B[1] -= eta2*eta2*PB; 													\
	B[1] *= 2.0*zeta; 														\
	B[2] = eta1*QP*(QC + QD) + 3.0*eta1*eta2 - PB*QP*eta2; 					\
	B[2] *= -4.0*QP*eta2*zeta*zeta; 										\
	QP *= QP*QP; 															\
	zeta *= zeta*zeta; 														\
	B[3] = 8.0*QP*eta1*eta2*eta2*zeta;

#define ERI2000(Ax,Px,Qx,eta1,zeta) 										\
	PA = Px - Ax; 															\
	QP = Qx - Px; 															\
	B[0] = 0.5*eta1 + PA*PA; 												\
	B[1] = -0.125*zeta*eta1*(eta1 - 4 * QP*PA); 							\
	B[2] = 0.25*QP*zeta*eta1; 												\
	B[2] = B[2] * B[2];

#define ERI1100(Ax,Bx,Px,Qx,eta1,zeta) 										\
	PA = Px - Ax; 															\
	PB = Px - Bx; 															\
	QP = Qx - Px; 															\
	B[0] = 0.5*eta1 + PA*PB; 												\
	B[1] = -0.125*zeta*eta1*(eta1 - 2 * QP*(PA + PB)); 						\
	B[2] = 0.25*QP*eta1*zeta; 												\
	B[2] = B[2] * B[2];

#define ERI1010(Ax,Cx,Px,Qx,eta1,eta2,zeta) 								\
	PA = Px - Ax; 															\
	QC = Qx - Cx; 															\
	QP = Qx - Px; 															\
	B[0] = PA*QC; 															\
	B[1] = 0.125*zeta*(eta1*eta2 + 2 * QP*(QC*eta1 - PA*eta2)); 			\
	B[2] = 0.25*QP*zeta; 													\
	B[2] = -eta1*eta2*B[2] * B[2];

#define ERI1000(Ax,Px,Qx,eta1,zeta) 										\
	B[0] = Px - Ax; 														\
	B[1] = 0.25*(Qx - Px)*eta1*zeta;

	/*
	I=l1+l2+l3+l4;
	switch(I){
	case 2: if(l1==2){         ERI2000(Ax,Px,Qx,eta1,zeta);    return; //(20|00)
	}else if(l2==2){   ERI2000(Bx,Px,Qx,eta1,zeta);    return; //(02|00)
	}else if(l3==2){   ERI2000(Cx,Qx,Px,eta2,zeta);    return; //(00|20)
	}else if(l4==2){   ERI2000(Dx,Qx,Px,eta2,zeta);    return; //(00|02)
	}else if(!(l3+l4)){ERI1100(Ax,Bx,Px,Qx,eta1,zeta); return; //(11|00)
	}else if(!(l1+l2)){ERI1100(Cx,Dx,Qx,Px,eta2,zeta); return; //(00|11)
	}else if(!(l2+l4)){ERI1010(Ax,Cx,Px,Qx,eta1,eta2,zeta); return; //(10|10)
	}else if(!(l1+l4)){ERI1010(Bx,Cx,Px,Qx,eta1,eta2,zeta); return; //(01|10)
	}else if(!(l2+l3)){ERI1010(Ax,Dx,Px,Qx,eta1,eta2,zeta); return; //(10|01)
	}else if(!(l1+l3)){ERI1010(Bx,Dx,Px,Qx,eta1,eta2,zeta); return; //(01|01)
	}
	break;

	case 1: if(l1==1){          ERI1000(Ax,Px,Qx,eta1,zeta); return; //(10|00)
	}else if(l2==1){    ERI1000(Bx,Px,Qx,eta1,zeta); return; //(01|00)
	}else if(l3==1){    ERI1000(Cx,Qx,Px,eta2,zeta); return; //(00|10)
	}else{              ERI1000(Dx,Qx,Px,eta2,zeta); return; //(00|01)
	}
	break;

	case 0: B[0] = 1.0; return;
	break;
	}
	*/

	switch (64 * l1 + 16 * l2 + 4 * l3 + l4){
	case 64 * 0 + 16 * 0 + 4 * 0 + 0: B[0] = 1.0; return; break;

	case 64 * 1 + 16 * 0 + 4 * 0 + 0: ERI1000(Ax, Px, Qx, eta1, zeta); return; break;
	case 64 * 0 + 16 * 1 + 4 * 0 + 0: ERI1000(Bx, Px, Qx, eta1, zeta); return; break;
	case 64 * 0 + 16 * 0 + 4 * 1 + 0: ERI1000(Cx, Qx, Px, eta2, zeta); return; break;
	case 64 * 0 + 16 * 0 + 4 * 0 + 1: ERI1000(Dx, Qx, Px, eta2, zeta); return; break;

	case 64 * 2 + 16 * 0 + 4 * 0 + 0: ERI2000(Ax, Px, Qx, eta1, zeta); return; break;
	case 64 * 0 + 16 * 2 + 4 * 0 + 0: ERI2000(Bx, Px, Qx, eta1, zeta); return; break;
	case 64 * 0 + 16 * 0 + 4 * 2 + 0: ERI2000(Cx, Qx, Px, eta2, zeta); return; break;
	case 64 * 0 + 16 * 0 + 4 * 0 + 2: ERI2000(Dx, Qx, Px, eta2, zeta); return; break;

	case 64 * 1 + 16 * 1 + 4 * 0 + 0: ERI1100(Ax, Bx, Px, Qx, eta1, zeta);      return; break;
	case 64 * 0 + 16 * 0 + 4 * 1 + 1: ERI1100(Cx, Dx, Qx, Px, eta2, zeta);      return; break;
	case 64 * 1 + 16 * 0 + 4 * 1 + 0: ERI1010(Ax, Cx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 0 + 16 * 1 + 4 * 1 + 0: ERI1010(Bx, Cx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 1 + 16 * 0 + 4 * 0 + 1: ERI1010(Ax, Dx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 0 + 16 * 1 + 4 * 0 + 1: ERI1010(Bx, Dx, Px, Qx, eta1, eta2, zeta); return; break;

	case 64 * 0 + 16 * 1 + 4 * 1 + 1: ERI0111(Bx, Cx, Dx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 1 + 16 * 0 + 4 * 1 + 1: ERI0111(Ax, Cx, Dx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 1 + 16 * 1 + 4 * 0 + 1: ERI0111(Dx, Ax, Bx, Qx, Px, eta2, eta1, zeta); return; break;
	case 64 * 1 + 16 * 1 + 4 * 1 + 0: ERI0111(Cx, Ax, Bx, Qx, Px, eta2, eta1, zeta); return; break;

	case 64 * 2 + 16 * 1 + 4 * 0 + 0: ERI2100(Ax, Bx, Px, Qx, eta1, zeta); return; break;
	case 64 * 1 + 16 * 2 + 4 * 0 + 0: ERI2100(Bx, Ax, Px, Qx, eta1, zeta); return; break;
	case 64 * 0 + 16 * 0 + 4 * 2 + 1: ERI2100(Cx, Dx, Qx, Px, eta2, zeta); return; break;
	case 64 * 0 + 16 * 0 + 4 * 1 + 2: ERI2100(Dx, Cx, Qx, Px, eta2, zeta); return; break;

	case 64 * 2 + 16 * 0 + 4 * 1 + 0: ERI2010(Ax, Cx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 2 + 16 * 0 + 4 * 0 + 1: ERI2010(Ax, Dx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 0 + 16 * 2 + 4 * 1 + 0: ERI2010(Bx, Cx, Px, Qx, eta1, eta2, zeta); return; break;
	case 64 * 0 + 16 * 2 + 4 * 0 + 1: ERI2010(Bx, Dx, Px, Qx, eta1, eta2, zeta); return; break;

	case 64 * 1 + 16 * 0 + 4 * 2 + 0: ERI2010(Cx, Ax, Qx, Px, eta2, eta1, zeta); return; break;
	case 64 * 1 + 16 * 0 + 4 * 0 + 2: ERI2010(Dx, Ax, Qx, Px, eta2, eta1, zeta); return; break;
	case 64 * 0 + 16 * 1 + 4 * 2 + 0: ERI2010(Cx, Bx, Qx, Px, eta2, eta1, zeta); return; break;
	case 64 * 0 + 16 * 1 + 4 * 0 + 2: ERI2010(Dx, Bx, Qx, Px, eta2, eta1, zeta); return; break;

	case 64 * 1 + 16 * 1 + 4 * 1 + 1: ERI1111(Ax, Bx, Cx, Dx, Px, Qx, eta1, eta2, zeta); return; break;
	}
}

#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+	\
	(y1 - y2)*(y1 - y2) +							\
	(z1 - z2)*(z1 - z2))
__device__ float eri_device(
		float xa, float ya, float za, float norma,
		int la, int ma, int na, float alphaa,
		float xb, float yb, float zb, float normb,
		int lb, int mb, int nb, float alphab,
		float xc, float yc, float zc, float normc,
		int lc, int mc, int nc, float alphac,
		float xd, float yd, float zd, float normd,
		int ld, int md, int nd, float alphad) {

	float Bx[4 * MAXL], By[4 * MAXL], Bz[4 * MAXL], F[4 * MAXL];
	float rab2, rcd2, rpq2, xp, yp, zp, xq, yq, zq, sum;
	float eta1, eta2, zeta;
	int kX, kY, kZ;
	float t;

	eta1 = 1.0 / (alphaa + alphab);
	eta2 = 1.0 / (alphac + alphad);
	zeta = 4.0*(alphaa + alphab)*(alphac + alphad)
		/ (alphaa + alphab + alphac + alphad);

	xp = (alphaa*xa + alphab*xb)*eta1;
	yp = (alphaa*ya + alphab*yb)*eta1;
	zp = (alphaa*za + alphab*zb)*eta1;
	xq = (alphac*xc + alphad*xd)*eta2;
	yq = (alphac*yc + alphad*yd)*eta2;
	zq = (alphac*zc + alphad*zd)*eta2;

	rab2 = DIST2(xa, ya, za, xb, yb, zb);
	rcd2 = DIST2(xc, yc, zc, xd, yd, zd);
	rpq2 = DIST2(xp, yp, zp, xq, yq, zq);

	// compute fgamma
	t = 0.25*rpq2*zeta;
	fgamma_set_device(la + lb + lc + ld + ma + mb + mc + md + na + nb + nc + nd, t, F);

	// integrate in 3 directions
	eriB_1d_device(Bx, la, lb, lc, ld, xp, xa, xb, xq, xc, xd, eta1, eta2, zeta);
	eriB_1d_device(By, ma, mb, mc, md, yp, ya, yb, yq, yc, yd, eta1, eta2, zeta);
	eriB_1d_device(Bz, na, nb, nc, nd, zp, za, zb, zq, zc, zd, eta1, eta2, zeta);

	sum = 0.0;
	for (kX = 0; kX < la + lb + lc + ld + 1; kX++)
	for (kY = 0; kY < ma + mb + mc + md + 1; kY++)
	for (kZ = 0; kZ < na + nb + nc + nd + 1; kZ++)
		sum += Bx[kX] * By[kY] * Bz[kZ] * F[kX + kY + kZ];

#define TWO_PI_POW2_5 34.986836655249725
	return  TWO_PI_POW2_5
		*eta1*eta2
		/ sqrt(alphaa + alphab + alphac + alphad)
		*exp(-alphaa*alphab*rab2*eta1
		- alphac*alphad*rcd2*eta2)
		*sum
		*norma*normb*normc*normd;
#undef  TWO_PI_POW2_5
}
#undef MAXL

__device__ float contr_eri_device(
		int lena, float *aexps, float *acoefs, float *anorms,
		float xa, float ya, float za, int la, int ma, int na,
		int lenb, float *bexps, float *bcoefs, float *bnorms,
		float xb, float yb, float zb, int lb, int mb, int nb,
		int lenc, float *cexps, float *ccoefs, float *cnorms,
		float xc, float yc, float zc, int lc, int mc, int nc,
		int lend, float *dexps, float *dcoefs, float *dnorms,
		float xd, float yd, float zd, int ld, int md, int nd) {

	float val = 0.0;
	int i, j, k, l;
	float EE;

	// proceed from highest exponent value
	for(i = 0; i<lena; i++)
	for(j = 0; j<lenb; j++)
	for(k = 0; k<lenc; k++)
	for(l = 0; l<lend; l++) {
		// compute element
		EE = acoefs[i] * bcoefs[j] * ccoefs[k] * dcoefs[l]
			* eri_device(xa, ya, za, anorms[i], la, ma, na, aexps[i],
			xb, yb, zb, bnorms[j], lb, mb, nb, bexps[j],
			xc, yc, zc, cnorms[k], lc, mc, nc, cexps[k],
			xd, yd, zd, dnorms[l], ld, md, nd, dexps[l]);
		val += EE;
	}

	return val;
}

__global__ void GTO_JK_Matrix_kernel(
		int nBasis, float *P,
		struct GTOBasis_t *gto,
		float *Schwarz, float cutoff, float *G) {

	int p, q, i, j, maxj;
	float upBound;
	float EE;

	int idx, idy, maxx;
	maxx = nBasis*(nBasis+1)/2;
	idx = blockIdx.x * blockDim.x + threadIdx.x;
	idy = blockIdx.y * blockDim.y + threadIdx.y;

	if(idx >= (int)(maxx/2+1) || idy >= maxx)
		return;

	if(idx > idy) {
		idx = maxx - idx;
		idy = idx + idy;
	}

	p = (int)(sqrt(1.0+8.0*idx) / 2.0 - 0.5);
	q = idx - (p+1)*p/2;
	i = (int)(sqrt(1.0+8.0*idy) / 2.0 - 0.5);
	j = idy - (i+1)*i/2;
//	if(i == p){
//		maxj = q;
//	}else{
//		maxj = i;
//	}
//
//	if(p >= nBasis || q > p || i > p ||j > maxj)
//		return;

#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
	upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
	if(upBound < SCHWARZ_CUTOFF) return;
	upBound = sqrt(upBound);
	if(fabs(upBound*P(p,q)) < cutoff &&
		fabs(upBound*P(p,i)) < cutoff &&
		fabs(upBound*P(p,j)) < cutoff &&
		fabs(upBound*P(q,i)) < cutoff &&
		fabs(upBound*P(q,j)) < cutoff &&
		fabs(upBound*P(i,j)) < cutoff) return;

	EE =
		contr_eri_device(
			gto[p].nContract,
			gto[p].exp, gto[p].coef, gto[p].norm,
			gto[p].x0, gto[p].y0, gto[p].z0,
			gto[p].l, gto[p].m, gto[p].n,
			gto[q].nContract,
			gto[q].exp, gto[q].coef, gto[q].norm,
			gto[q].x0, gto[q].y0, gto[q].z0,
			gto[q].l, gto[q].m, gto[q].n,
			gto[i].nContract,
			gto[i].exp, gto[i].coef, gto[i].norm,
			gto[i].x0, gto[i].y0, gto[i].z0,
			gto[i].l, gto[i].m, gto[i].n,
			gto[j].nContract,
			gto[j].exp, gto[j].coef, gto[j].norm,
			gto[j].x0, gto[j].y0, gto[j].z0,
			gto[j].l, gto[j].m, gto[j].n);

	// perform matrix symmetry rotations
	if((p==q)&&(i==j)&&(p==i)){  // all same
		atomicAdd(G+nBasis*p + q, EE*P(i,j));
	}else if((p==q)&&(i==j)){    // 2 pairs
		atomicAdd(G+nBasis*p + p, 2.0*EE*P(i,i));
		atomicAdd(G+nBasis*i + i, 2.0*EE*P(p,p));

		atomicAdd(G+nBasis*p + i, -EE*P(i,p));
		atomicAdd(G+nBasis*i + p, -EE*P(i,p));
	}else if(p==q){              // pq pair
		atomicAdd(G+nBasis*p + p, 4.0*EE*P(i,j));

		atomicAdd(G+nBasis*i + j, 2.0*EE*P(p,p));
		atomicAdd(G+nBasis*j + i, 2.0*EE*P(p,p));

		atomicAdd(G+nBasis*p + j, -EE*P(i,p));
		atomicAdd(G+nBasis*j + p, -EE*P(i,p));

		atomicAdd(G+nBasis*i + p, -EE*P(p,j));
		atomicAdd(G+nBasis*p + i, -EE*P(p,j));
	}else if(i==j){              // ij pair
		atomicAdd(G+nBasis*i + i, 4.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + q, 2.0*EE*P(i,i));
		atomicAdd(G+nBasis*q + p, 2.0*EE*P(i,i));

		atomicAdd(G+nBasis*p + i, -EE*P(i,q));
		atomicAdd(G+nBasis*i + p, -EE*P(i,q));

		atomicAdd(G+nBasis*q + i, -EE*P(i,p));
		atomicAdd(G+nBasis*i + q, -EE*P(i,p));
	}else if((p==i)&&(q==j)){    // pi-qj pair
		atomicAdd(G+nBasis*p + q, 3.0*EE*P(p,q));
		atomicAdd(G+nBasis*q + p, 3.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + p, -EE*P(q,q));
		atomicAdd(G+nBasis*q + q, -EE*P(p,p));
	}else{                       // all distinct
		atomicAdd(G+nBasis*p + q, 4.0*EE*P(i,j));
		atomicAdd(G+nBasis*q + p, 4.0*EE*P(i,j));

		atomicAdd(G+nBasis*i + j, 4.0*EE*P(p,q));
		atomicAdd(G+nBasis*j + i, 4.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + j, -EE*P(i,q));
		atomicAdd(G+nBasis*j + p, -EE*P(i,q));

		atomicAdd(G+nBasis*p + i, -EE*P(j,q));
		atomicAdd(G+nBasis*i + p, -EE*P(j,q));

		atomicAdd(G+nBasis*q + i, -EE*P(j,p));
		atomicAdd(G+nBasis*i + q, -EE*P(j,p));

		atomicAdd(G+nBasis*q + j, -EE*P(i,p));
		atomicAdd(G+nBasis*j + q, -EE*P(i,p));
	}
#undef G
#undef P
}

__global__ void eri_device_kernel(int *count,
		int a, int b, int c, int d,
		int nBasis, float *P,
		struct GTOBasis_t *gto, float *G) {

	int i, j, k, l;
	float EE;
	int id = threadIdx.x;
	l = id % 3;
	id = id / 3;
	k = id % 3;
	id = id / 3;
	j = id % 3;
	id = id / 3;
	i = id % 3;
//	i = blockIdx.x;
//	j = blockIdx.y;
//	k = threadIdx.x;
//	l = threadIdx.y;
//	atomicAdd(count, 1);
	EE = gto[a].coef[i] * gto[b].coef[j] * gto[c].coef[k] * gto[d].coef[l]
	    * eri_device(gto[a].x0, gto[a].y0, gto[a].z0, gto[a].norm[i], gto[a].l, gto[a].m, gto[a].n, gto[a].exp[i],
		gto[b].x0, gto[b].y0, gto[b].z0, gto[b].norm[j], gto[b].l, gto[b].m, gto[b].n, gto[b].exp[j],
		gto[c].x0, gto[c].y0, gto[c].z0, gto[c].norm[k], gto[c].l, gto[c].m, gto[c].n, gto[c].exp[k],
		gto[d].x0, gto[d].y0, gto[d].z0, gto[d].norm[l], gto[d].l, gto[d].m, gto[d].n, gto[d].exp[l]);

#define P(x,y) P[nBasis*x+y]
	if((a==b)&&(c==d)&&(a==c)){  // all same
		atomicAdd(G+nBasis*a + b, EE*P(c,d));
	}else if((a==b)&&(c==d)){    // 2 pairs
		atomicAdd(G+nBasis*a + a, 2.0*EE*P(c,c));
		atomicAdd(G+nBasis*c + c, 2.0*EE*P(a,a));

		atomicAdd(G+nBasis*a + c, -EE*P(c,a));
		atomicAdd(G+nBasis*c + a, -EE*P(c,a));
	}else if(a==b){              // ab pair
		atomicAdd(G+nBasis*a + a, 4.0*EE*P(c,d));

		atomicAdd(G+nBasis*c + d, 2.0*EE*P(a,a));
		atomicAdd(G+nBasis*d + c, 2.0*EE*P(a,a));

		atomicAdd(G+nBasis*a + d, -EE*P(c,a));
		atomicAdd(G+nBasis*d + a, -EE*P(c,a));

		atomicAdd(G+nBasis*c + a, -EE*P(a,d));
		atomicAdd(G+nBasis*a + c, -EE*P(a,d));
	}else if(c==d){              // cd pair
		atomicAdd(G+nBasis*c + c, 4.0*EE*P(a,b));

		atomicAdd(G+nBasis*a + b, 2.0*EE*P(c,c));
		atomicAdd(G+nBasis*b + a, 2.0*EE*P(c,c));

		atomicAdd(G+nBasis*a + c, -EE*P(c,b));
		atomicAdd(G+nBasis*c + a, -EE*P(c,b));

		atomicAdd(G+nBasis*b + c, -EE*P(c,a));
		atomicAdd(G+nBasis*c + b, -EE*P(c,a));
	}else if((a==c)&&(b==d)){    // ac-bd pair
		atomicAdd(G+nBasis*a + b, 3.0*EE*P(a,b));
		atomicAdd(G+nBasis*b + a, 3.0*EE*P(a,b));

		atomicAdd(G+nBasis*a + a, -EE*P(b,b));
		atomicAdd(G+nBasis*b + b, -EE*P(a,a));
	}else{                       // all distinct
		atomicAdd(G+nBasis*a + b, 4.0*EE*P(c,d));
		atomicAdd(G+nBasis*b + a, 4.0*EE*P(c,d));

		atomicAdd(G+nBasis*c + d, 4.0*EE*P(a,b));
		atomicAdd(G+nBasis*d + c, 4.0*EE*P(a,b));

		atomicAdd(G+nBasis*a + d, -EE*P(c,b));
		atomicAdd(G+nBasis*d + a, -EE*P(c,b));

		atomicAdd(G+nBasis*a + c, -EE*P(d,b));
		atomicAdd(G+nBasis*c + a, -EE*P(d,b));

		atomicAdd(G+nBasis*b + c, -EE*P(d,a));
		atomicAdd(G+nBasis*c + b, -EE*P(d,a));

		atomicAdd(G+nBasis*b + d, -EE*P(c,a));
		atomicAdd(G+nBasis*d + b, -EE*P(c,a));
	}
#undef P
}

__global__ void GTO_JK_Matrix_kernel1(
		unsigned int nEE, unsigned int *iter_p, unsigned int *iter_q, unsigned int *iter_i, unsigned int *iter_j,
		int nBasis, float *P,
		struct GTOBasis_t *gto, float *G) {

	int p, q, i, j;
	float EE;

	int id;
	id = blockDim.x * blockIdx.x + threadIdx.x;

	if(id >= nEE)
		return;
	p = iter_p[id];
	q = iter_q[id];
	i = iter_i[id];
	j = iter_j[id];

#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
	EE =
		contr_eri_device(
			gto[p].nContract,
			gto[p].exp, gto[p].coef, gto[p].norm,
			gto[p].x0, gto[p].y0, gto[p].z0,
			gto[p].l, gto[p].m, gto[p].n,
			gto[q].nContract,
			gto[q].exp, gto[q].coef, gto[q].norm,
			gto[q].x0, gto[q].y0, gto[q].z0,
			gto[q].l, gto[q].m, gto[q].n,
			gto[i].nContract,
			gto[i].exp, gto[i].coef, gto[i].norm,
			gto[i].x0, gto[i].y0, gto[i].z0,
			gto[i].l, gto[i].m, gto[i].n,
			gto[j].nContract,
			gto[j].exp, gto[j].coef, gto[j].norm,
			gto[j].x0, gto[j].y0, gto[j].z0,
			gto[j].l, gto[j].m, gto[j].n);

	// perform matrix symmetry rotations
	if((p==q)&&(i==j)&&(p==i)){  // all same
		atomicAdd(G+nBasis*p + q, EE*P(i,j));
	}else if((p==q)&&(i==j)){    // 2 pairs
		atomicAdd(G+nBasis*p + p, 2.0*EE*P(i,i));
		atomicAdd(G+nBasis*i + i, 2.0*EE*P(p,p));

		atomicAdd(G+nBasis*p + i, -EE*P(i,p));
		atomicAdd(G+nBasis*i + p, -EE*P(i,p));
	}else if(p==q){              // pq pair
		atomicAdd(G+nBasis*p + p, 4.0*EE*P(i,j));

		atomicAdd(G+nBasis*i + j, 2.0*EE*P(p,p));
		atomicAdd(G+nBasis*j + i, 2.0*EE*P(p,p));

		atomicAdd(G+nBasis*p + j, -EE*P(i,p));
		atomicAdd(G+nBasis*j + p, -EE*P(i,p));

		atomicAdd(G+nBasis*i + p, -EE*P(p,j));
		atomicAdd(G+nBasis*p + i, -EE*P(p,j));
	}else if(i==j){              // ij pair
		atomicAdd(G+nBasis*i + i, 4.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + q, 2.0*EE*P(i,i));
		atomicAdd(G+nBasis*q + p, 2.0*EE*P(i,i));

		atomicAdd(G+nBasis*p + i, -EE*P(i,q));
		atomicAdd(G+nBasis*i + p, -EE*P(i,q));

		atomicAdd(G+nBasis*q + i, -EE*P(i,p));
		atomicAdd(G+nBasis*i + q, -EE*P(i,p));
	}else if((p==i)&&(q==j)){    // pi-qj pair
		atomicAdd(G+nBasis*p + q, 3.0*EE*P(p,q));
		atomicAdd(G+nBasis*q + p, 3.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + p, -EE*P(q,q));
		atomicAdd(G+nBasis*q + q, -EE*P(p,p));
	}else{                       // all distinct
		atomicAdd(G+nBasis*p + q, 4.0*EE*P(i,j));
		atomicAdd(G+nBasis*q + p, 4.0*EE*P(i,j));

		atomicAdd(G+nBasis*i + j, 4.0*EE*P(p,q));
		atomicAdd(G+nBasis*j + i, 4.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + j, -EE*P(i,q));
		atomicAdd(G+nBasis*j + p, -EE*P(i,q));

		atomicAdd(G+nBasis*p + i, -EE*P(j,q));
		atomicAdd(G+nBasis*i + p, -EE*P(j,q));

		atomicAdd(G+nBasis*q + i, -EE*P(j,p));
		atomicAdd(G+nBasis*i + q, -EE*P(j,p));

		atomicAdd(G+nBasis*q + j, -EE*P(i,p));
		atomicAdd(G+nBasis*j + q, -EE*P(i,p));
	}
#undef G
#undef P
}

__device__ unsigned int dev_nEE = 0;

__global__ void integral_screening_kernel(
		unsigned int *iter_p, unsigned int *iter_q, unsigned int *iter_i, unsigned int *iter_j,
		int nBasis, float *P,
		float *Schwarz, float cutoff) {

	unsigned int p, q, i, j, maxj;
	float upBound;

	unsigned int idx, idy;
	idx = blockIdx.x * blockDim.x + threadIdx.x;
	idy = blockIdx.y * blockDim.y + threadIdx.y;

	p = idx / nBasis;
	q = idx % nBasis;
	i = idy / nBasis;
	j = idy % nBasis;

	if(i == p)
		maxj = q;
	else
		maxj = i;

	if((p >= nBasis) || (q > p) || (i > p) || (j > maxj))
		return;

#define P(x,y) P[nBasis*x+y]
	upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
	if(upBound < SCHWARZ_CUTOFF) return;
	upBound = sqrt(upBound);
	if(fabs(upBound*P(p,q)) < cutoff &&
		fabs(upBound*P(p,i)) < cutoff &&
		fabs(upBound*P(p,j)) < cutoff &&
		fabs(upBound*P(q,i)) < cutoff &&
		fabs(upBound*P(q,j)) < cutoff &&
		fabs(upBound*P(i,j)) < cutoff) return;

	idx = atomicAdd(&dev_nEE, 1);
	iter_p[idx] = p;
	iter_q[idx] = q;
	iter_i[idx] = i;
	iter_j[idx] = j;

#undef P
}

__global__ void GTO_JK_Matrix_kernel2(
		unsigned int *iter_p, unsigned int *iter_q, unsigned int *iter_i, unsigned int *iter_j,
		int nBasis, float *P,
		struct GTOBasis_t *gto, float *G) {

	unsigned int p, q, i, j;
	float EE;

	unsigned int id;
	id = blockDim.x * blockIdx.x + threadIdx.x;

	if(id >= dev_nEE)
		return;
	p = iter_p[id];
	q = iter_q[id];
	i = iter_i[id];
	j = iter_j[id];

#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]

	EE =
		contr_eri_device(
			gto[p].nContract,
			gto[p].exp, gto[p].coef, gto[p].norm,
			gto[p].x0, gto[p].y0, gto[p].z0,
			gto[p].l, gto[p].m, gto[p].n,
			gto[q].nContract,
			gto[q].exp, gto[q].coef, gto[q].norm,
			gto[q].x0, gto[q].y0, gto[q].z0,
			gto[q].l, gto[q].m, gto[q].n,
			gto[i].nContract,
			gto[i].exp, gto[i].coef, gto[i].norm,
			gto[i].x0, gto[i].y0, gto[i].z0,
			gto[i].l, gto[i].m, gto[i].n,
			gto[j].nContract,
			gto[j].exp, gto[j].coef, gto[j].norm,
			gto[j].x0, gto[j].y0, gto[j].z0,
			gto[j].l, gto[j].m, gto[j].n);

	// perform matrix symmetry rotations
	if((p==q)&&(i==j)&&(p==i)){  // all same
		atomicAdd(G+nBasis*p + q, EE*P(i,j));
	}else if((p==q)&&(i==j)){    // 2 pairs
		atomicAdd(G+nBasis*p + p, 2.0*EE*P(i,i));
		atomicAdd(G+nBasis*i + i, 2.0*EE*P(p,p));

		atomicAdd(G+nBasis*p + i, -EE*P(i,p));
		atomicAdd(G+nBasis*i + p, -EE*P(i,p));
	}else if(p==q){              // pq pair
		atomicAdd(G+nBasis*p + p, 4.0*EE*P(i,j));

		atomicAdd(G+nBasis*i + j, 2.0*EE*P(p,p));
		atomicAdd(G+nBasis*j + i, 2.0*EE*P(p,p));

		atomicAdd(G+nBasis*p + j, -EE*P(i,p));
		atomicAdd(G+nBasis*j + p, -EE*P(i,p));

		atomicAdd(G+nBasis*i + p, -EE*P(p,j));
		atomicAdd(G+nBasis*p + i, -EE*P(p,j));
	}else if(i==j){              // ij pair
		atomicAdd(G+nBasis*i + i, 4.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + q, 2.0*EE*P(i,i));
		atomicAdd(G+nBasis*q + p, 2.0*EE*P(i,i));

		atomicAdd(G+nBasis*p + i, -EE*P(i,q));
		atomicAdd(G+nBasis*i + p, -EE*P(i,q));

		atomicAdd(G+nBasis*q + i, -EE*P(i,p));
		atomicAdd(G+nBasis*i + q, -EE*P(i,p));
	}else if((p==i)&&(q==j)){    // pi-qj pair
		atomicAdd(G+nBasis*p + q, 3.0*EE*P(p,q));
		atomicAdd(G+nBasis*q + p, 3.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + p, -EE*P(q,q));
		atomicAdd(G+nBasis*q + q, -EE*P(p,p));
	}else{                       // all distinct
		atomicAdd(G+nBasis*p + q, 4.0*EE*P(i,j));
		atomicAdd(G+nBasis*q + p, 4.0*EE*P(i,j));

		atomicAdd(G+nBasis*i + j, 4.0*EE*P(p,q));
		atomicAdd(G+nBasis*j + i, 4.0*EE*P(p,q));

		atomicAdd(G+nBasis*p + j, -EE*P(i,q));
		atomicAdd(G+nBasis*j + p, -EE*P(i,q));

		atomicAdd(G+nBasis*p + i, -EE*P(j,q));
		atomicAdd(G+nBasis*i + p, -EE*P(j,q));

		atomicAdd(G+nBasis*q + i, -EE*P(j,p));
		atomicAdd(G+nBasis*i + q, -EE*P(j,p));

		atomicAdd(G+nBasis*q + j, -EE*P(i,p));
		atomicAdd(G+nBasis*j + q, -EE*P(i,p));
	}
#undef G
#undef P
}

#define BLOCK_W 16
extern "C" void GTO_JK_Matrix_CUDA(
		int nBasis, float *P,
		struct GTOBasis_t *gto,
		float *Schwarz, float cutoff, float *G) {

	dim3 grids((nBasis*(nBasis+1)/4 + 1 + BLOCK_W-1) / BLOCK_W, (nBasis*(nBasis+1)/2 + BLOCK_W-1) / BLOCK_W, 1);
	dim3 blocks(BLOCK_W, BLOCK_W);

	float *dev_P;
	struct GTOBasis_t *dev_gto;
	float *dev_Schwarz;
	float *dev_G;

	CUDA_CHECK_RETURN(cudaSetDevice(0));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_P, sizeof(float) * nBasis * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_gto, sizeof(struct GTOBasis_t) * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_Schwarz, sizeof(float) * nBasis * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_G, sizeof(float) * nBasis * nBasis));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_P, P, sizeof(float) * nBasis * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_gto, gto, sizeof(struct GTOBasis_t) * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_Schwarz, Schwarz, sizeof(float) * nBasis * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(cudaMemset(dev_G, 0, sizeof(float) * nBasis * nBasis));

	GTO_JK_Matrix_kernel<<<grids, blocks>>>(nBasis, dev_P, dev_gto, dev_Schwarz, cutoff, dev_G);

//	CUDA_CHECK_RETURN(cudaThreadSynchronize());	// Wait for the GPU launched work to complete
	CUDA_CHECK_RETURN(cudaGetLastError());
	CUDA_CHECK_RETURN(cudaMemcpy(G, dev_G, sizeof(float) * nBasis * nBasis, cudaMemcpyDeviceToHost));

	CUDA_CHECK_RETURN(cudaFree((void*) dev_P));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_gto));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_Schwarz));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_G));
	CUDA_CHECK_RETURN(cudaDeviceReset());
}

extern "C" void GTO_JK_Matrix_CUDA1(
		unsigned int nEE, unsigned int *iter_p, unsigned int *iter_q, unsigned int *iter_i, unsigned int *iter_j,
		int nBasis, float *P,
		struct GTOBasis_t *gto, float *G) {

	float *dev_P;
	struct GTOBasis_t *dev_gto;
	float *dev_G;
	unsigned int *dev_iter_p, *dev_iter_q, *dev_iter_i, *dev_iter_j;

	struct timeval startTime;
	struct timeval endTime;

	gettimeofday(&startTime, NULL);

	CUDA_CHECK_RETURN(cudaSetDevice(0));

	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_p, sizeof(unsigned int) * nEE));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_q, sizeof(unsigned int) * nEE));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_i, sizeof(unsigned int) * nEE));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_j, sizeof(unsigned int) * nEE));

	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_P, sizeof(float) * nBasis * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_gto, sizeof(struct GTOBasis_t) * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_G, sizeof(float) * nBasis * nBasis));

	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_iter_p, iter_p, sizeof(unsigned int) * nEE, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_iter_q, iter_q, sizeof(unsigned int) * nEE, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_iter_i, iter_i, sizeof(unsigned int) * nEE, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_iter_j, iter_j, sizeof(unsigned int) * nEE, cudaMemcpyHostToDevice));

	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_P, P, sizeof(float) * nBasis * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_gto, gto, sizeof(struct GTOBasis_t) * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(cudaMemset(dev_G, 0, sizeof(float) * nBasis * nBasis));

	GTO_JK_Matrix_kernel1<<<(nEE + 1023) / 1024, 1024>>>(nEE, dev_iter_p, dev_iter_q, dev_iter_i, dev_iter_j, nBasis, dev_P, dev_gto, dev_G);

//	CUDA_CHECK_RETURN(cudaThreadSynchronize());	// Wait for the GPU launched work to complete
	CUDA_CHECK_RETURN(cudaGetLastError());
	CUDA_CHECK_RETURN(cudaMemcpy(G, dev_G, sizeof(float) * nBasis * nBasis, cudaMemcpyDeviceToHost));

	gettimeofday(&endTime, NULL);
	printf(" GPU cost time: %ld us.\n", 1000000 * (endTime.tv_sec - startTime.tv_sec) + endTime.tv_usec - startTime.tv_usec);

	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_p));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_q));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_i));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_j));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_P));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_gto));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_G));
	CUDA_CHECK_RETURN(cudaDeviceReset());
}

#define MAX_N_EE 102400000
extern "C" void GTO_JK_Matrix_CUDA2(
		int nBasis, float *P,
		struct GTOBasis_t *gto,
		float *Schwarz, float cutoff, float *G) {

	dim3 grids((nBasis*nBasis + BLOCK_W-1) / BLOCK_W, (nBasis*nBasis + BLOCK_W-1) / BLOCK_W, 1);
	dim3 blocks(BLOCK_W, BLOCK_W);

	float *dev_P;
	struct GTOBasis_t *dev_gto;
	float *dev_Schwarz;
	unsigned int *dev_iter_p, *dev_iter_q, *dev_iter_i, *dev_iter_j;
	float *dev_G;

	CUDA_CHECK_RETURN(cudaSetDevice(0));

	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_p, sizeof(unsigned int) * MAX_N_EE));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_q, sizeof(unsigned int) * MAX_N_EE));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_i, sizeof(unsigned int) * MAX_N_EE));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_iter_j, sizeof(unsigned int) * MAX_N_EE));

	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_P, sizeof(float) * nBasis * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_gto, sizeof(struct GTOBasis_t) * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_Schwarz, sizeof(float) * nBasis * nBasis));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_G, sizeof(float) * nBasis * nBasis));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_P, P, sizeof(float) * nBasis * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_gto, gto, sizeof(struct GTOBasis_t) * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(dev_Schwarz, Schwarz, sizeof(float) * nBasis * nBasis, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(cudaMemset(dev_G, 0, sizeof(float) * nBasis * nBasis));

	CUDA_CHECK_RETURN(cudaMemset(dev_iter_p, 0, sizeof(unsigned int) * MAX_N_EE));
	CUDA_CHECK_RETURN(cudaMemset(dev_iter_q, 0, sizeof(unsigned int) * MAX_N_EE));
	CUDA_CHECK_RETURN(cudaMemset(dev_iter_i, 0, sizeof(unsigned int) * MAX_N_EE));
	CUDA_CHECK_RETURN(cudaMemset(dev_iter_j, 0, sizeof(unsigned int) * MAX_N_EE));

	integral_screening_kernel<<<grids, blocks>>>(dev_iter_p, dev_iter_q, dev_iter_i, dev_iter_j, nBasis,
			dev_P, dev_Schwarz, cutoff);
	cudaDeviceSynchronize();

	GTO_JK_Matrix_kernel2<<<65535, 1024>>>(dev_iter_p, dev_iter_q, dev_iter_i, dev_iter_j, nBasis,
			dev_P, dev_gto, dev_G);

//	CUDA_CHECK_RETURN(cudaThreadSynchronize());	// Wait for the GPU launched work to complete
	CUDA_CHECK_RETURN(cudaGetLastError());
	CUDA_CHECK_RETURN(cudaMemcpy(G, dev_G, sizeof(float) * nBasis * nBasis, cudaMemcpyDeviceToHost));

	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_p));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_q));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_i));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_iter_j));

	CUDA_CHECK_RETURN(cudaFree((void*) dev_P));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_gto));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_Schwarz));
	CUDA_CHECK_RETURN(cudaFree((void*) dev_G));
	CUDA_CHECK_RETURN(cudaDeviceReset());
}
#undef MAX_N_EE
