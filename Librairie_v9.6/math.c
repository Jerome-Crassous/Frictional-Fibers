#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

//---------------------------------------------------------------------------------------
double MyRand(){						//retourne un nombre entre 0 et 1
	long			k;
	static	long	idum=4325;		//une valeur au pif pour seed
	double			ans;
	idum ^= 123459876;
	k=idum/127773;
	idum=16807*(idum-k*127773)-2836*k;
	if (idum < 0) idum += 2147483647;
	ans=((double) idum)/2147483647.;
	idum ^= 123459876;
	return ans;
}

//---------------------------------------------------------------------------------------
double MyRandWithSeed(long seed) {						//retourne un nombre entre 0 et 1
	long			k;
	double			ans;

	seed ^= 123459876;
	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - 2836 * k;
	if (seed < 0) seed += 2147483647;

	seed ^= 123459876;
	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - 2836 * k;
	if (seed < 0) seed += 2147483647;

	ans = ((double)seed) / 2147483647.;
	return ans;
}

//---------------------------------------------------------------------------------------
double DotProduct(double* v1, double* v2) {
	return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

//---------------------------------------------------------------------------------------
double CosTheta(double* v1, double* v2) {
	return DotProduct(v1,v2)/Norm(v1)/Norm(v2);
}

//---------------------------------------------------------------------------------------
double Distance(double* v1, double* v2) {
	
	int		k;
	double  v[3];

	for (k = 0; k < 3; k++) 
		v[k]=v2[k]-v1[k];

	return Norm(v);
}

//---------------------------------------------------------------------------------------
double Norm(double* v) {
	return sqrt(DMax(0., DotProduct(v, v)));
}

//---------------------------------------------------------------------------------------
double Normalize(double* v) {

	int			k;
	double		norm;

	norm = Norm(v);
	for(k=0;k<3;k++) v[k]/=norm;
	return norm;
}

//---------------------------------------------------------------------------------------
void CrossProduct(double* v1, double* v2, double* cross) {

	cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
	cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
	cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return;
}

//---------------------------------------------------------------------------------------
void VectorMultiply(double* v, double scale) {
	v[0] *= scale, v[1] *= scale, v[2] *= scale;
	return;
}

//---------------------------------------------------------------------------------------
int IMin(int a, int b) {
	if (a < b)
		return a;
	else
		return b;
}

//---------------------------------------------------------------------------------------
int IMax(int a, int b) {
	if (a < b)
		return b;
	else
		return a;
}

//---------------------------------------------------------------------------------------
double DMin(double a, double b) {
	if (a < b)
		return a;
	else
		return b;
}

//---------------------------------------------------------------------------------------
double DMax(double a, double b) {
	if (a < b)
		return b;
	else
		return a;
}

//---------------------------------------------------------------------------------------
double DAbs(double a) {
	if (a < 0)
		return -a;
	else
		return a;
}

//---------------------------------------------------------------------------------------
double DSign(double a) {
	if (a > 0)
		return 1.;
	else
		return -1.;
}

//-------------------------------------------------------------------------
double Angle(double* v1, double* v2) {
	double aux;

	aux = fabs(DotProduct(v1, v2)) / sqrt(DotProduct(v1, v1) * DotProduct(v2, v2));
	return(acos(aux));
}

//----------------------------------------------------------------
void SetRandomOriention(double* e, double* m1) { //--------- genere 1 vecteurs unitaire perpendiculaire à e

	int				k;
	double			aux, eLoc[3]= { 0, 0, 0 };
		
	if (Norm(e) < EPS) {
		printf("Norm nulle dans SetRandomOrientation\n");
		exit(777);
	}

	//--------- copie de e normalisée
	for (k = 0; k < 3; k++) eLoc[k] = e[k];
	Normalize(eLoc);

	//--------- genere un m1 non parallele à eLoc
	for (;;) {
		for (k = 0; k < 3; k++) m1[k] = MyRand();
		aux = DotProduct(eLoc, m1);
		if (fabs(aux) > EPS) break;
	}

	//--------- m1 perpendiculaire à eLoc
	for (k = 0; k < 3; k++) 
		m1[k] -= aux * e[k];
	Normalize(m1);
}

//-----------------------------------------------------------------
void Rotate_Point_Around_Axis(double *x, double *pointOnAxis, double* direction, double angle){

	int		k;
	double	p[3], aux[3];
	double cosA = cos(angle);
	double sinA = sin(angle);

	for(k=0;k<3;k++) 
		p[k] = x[k] - pointOnAxis[k];
	CrossProduct(direction, p, aux);
	for (k = 0;k < 3;k++)
		x[k] = pointOnAxis[k] + p[k] * cosA + aux[k] * sinA + direction[k] * DotProduct(direction, p) * (1.0f - cosA);
}

//-----------------------------------------------------------------
void Shift_Point(double* x, double* shift) {

	int		k;

	for (k = 0;k < 3;k++)
		x[k] += shift[k];
}