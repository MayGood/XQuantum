/*
 * int.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef INT_H_
#define INT_H_

#define PRIMITIVE_CUTOFF 1.0E-15

float moment(
	float alpha1, int l1, int m1, int n1, // gto(1) exponent and angular index
	float xa, float ya, float za,       // gto(1) center
	float alpha2, int l2, int m2, int n2, // gto(2) exponent and angular index
	float xb, float yb, float zb,       // gto(2) center
	int mx, int my, int mz,                // moment angular index
	float xc, float yc, float zc);      // moment center

float overlap(float alpha1, int l1, int m1, int n1,
	float xa, float ya, float za,
	float alpha2, int l2, int m2, int n2,
	float xb, float yb, float zb);

float kinetic(float alpha1, int l1, int m1, int n1,
	float xa, float ya, float za,
	float alpha2, int l2, int m2, int n2,
	float xb, float yb, float zb);

float nai(float x1, float y1, float z1, float norm1,
	int l1, int m1, int n1, float alpha1,
	float x2, float y2, float z2, float norm2,
	int l2, int m2, int n2, float alpha2,
	float x3, float y3, float z3);

void eriB_1d(float *B, int l1, int l2, int l3, int l4,
	float Px, float Ax, float Bx,
	float Qx, float Cx, float Dx,
	float eta1, float eta2,
	float zeta);

float eri(float xa, float ya, float za, float norma,
	int la, int ma, int na, float alphaa,
	float xb, float yb, float zb, float normb,
	int lb, int mb, int nb, float alphab,
	float xc, float yc, float zc, float normc,
	int lc, int mc, int nc, float alphac,
	float xd, float yd, float zd, float normd,
	int ld, int md, int nd, float alphad);

float  genSetBxyzSxyzFf(
	float xa,  float ya,  float za, int maxa, float alphaa,
	float xb,  float yb,  float zb, int maxb, float alphab,
	float xc,  float yc,  float zc, int maxc, float alphac,
	float xd,  float yd,  float zd, int maxd, float alphad,
	float *Bx, float *By, float *Bz,
	float *Sx, float *Sy, float *Sz,
	float *F,int maxL);

float  genSetBxyzSxyzF(
	float xa, float ya, float za, int maxa, float alphaa,
	float xb, float yb, float zb, int maxb, float alphab,
	float xc, float yc, float zc, int maxc, float alphac,
	float xd, float yd, float zd, int maxd, float alphad,
	float *Bx, float *By, float *Bz,
	float *Sx, float *Sy, float *Sz,
	float *F, int maxL);

#endif /* INT_H_ */
