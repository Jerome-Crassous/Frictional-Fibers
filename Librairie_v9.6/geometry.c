#include <stdio.h>
#include <math.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

//---------------------------------------------------------------------
void Calculate_Relaxed_M1_Bar(fiber* theFiberPtr) {

	unsigned int			i, k, n;
	double					bis, T[3], dT[3], Darboux[3], aux[3];

	//---------------------------
	CalculateLeLib(theFiberPtr);

	//----------------- guess m1_bar[0]
	if (Norm(theFiberPtr->m1_bar[0]) < 1.e-3){
		for (k = 0; k < 3; k++)
			theFiberPtr->m1_bar[0][k] = MyRand();
	}
	bis = DotProduct(theFiberPtr->m1_bar[0], theFiberPtr->e[0]);
	for (k = 0; k < 3; k++)
		theFiberPtr->m1_bar[0][k] -= bis * theFiberPtr->e[0][k];
	Normalize(theFiberPtr->m1_bar[0]);

	//----------------- m1_bar[i]
	for (i = 1; i < theFiberPtr->n - 1; i++) {
		for (k = 0; k < 3; k++)
			T[k] = (theFiberPtr->xt[i + 1][k] - theFiberPtr->xt[i - 1][k]) / theFiberPtr->l0;
		Normalize(T);
		for (k = 0; k < 3; k++)
			dT[k] = (theFiberPtr->xt[i + 1][k] - 2. * theFiberPtr->xt[i][k] + theFiberPtr->xt[i - 1][k])
			/ theFiberPtr->l0 / theFiberPtr->l0;
		CrossProduct(T, dT, Darboux);

		CrossProduct(Darboux, theFiberPtr->m1_bar[i - 1], aux);
		for (k = 0; k < 3; k++)
			theFiberPtr->m1_bar[i][k] = theFiberPtr->m1_bar[i - 1][k] + aux[k] * theFiberPtr->l0;
		bis = DotProduct(theFiberPtr->m1_bar[i], theFiberPtr->e[i]);
		for (k = 0; k < 3; k++)
			theFiberPtr->m1_bar[i][k] -= bis * theFiberPtr->e[i][k];
		Normalize(theFiberPtr->m1_bar[i]);

		for (n = 0; n < 50; n++) {
			for (k = 0; k < 3; k++)
				theFiberPtr->m1_bar[i][k] = (theFiberPtr->m1_bar[i - 1][k] + theFiberPtr->m1_bar[i][k]) / 2.;
			CrossProduct(Darboux, theFiberPtr->m1_bar[i], aux);
			for (k = 0; k < 3; k++)
				theFiberPtr->m1_bar[i][k] = theFiberPtr->m1_bar[i - 1][k] + aux[k] * theFiberPtr->l0;
			bis = DotProduct(theFiberPtr->m1_bar[i], theFiberPtr->e[i]);
			for (k = 0; k < 3; k++)
				theFiberPtr->m1_bar[i][k] -= bis * theFiberPtr->e[i][k];
			Normalize(theFiberPtr->m1_bar[i]);
		}
	}
}

//---------------------------------------------------------------------
void Calculate_M1_From_M1Bar_And_Theta(fiber* theFiberPtr) {

	unsigned int			i, k;

	if (theFiberPtr->n == 1) {		//----------------- cas bille
		;
	}
	else {							//----------------- cas fibres
		for (i = 0; i < theFiberPtr->n; i++) {
			for (k = 0; k < 3; k++)
				theFiberPtr->m1[i][k] = theFiberPtr->m1_bar[i][k];

			Rotate_Point_Around_Axis(theFiberPtr->m1[i],
				theFiberPtr->xt[i], theFiberPtr->e[i], theFiberPtr->thetat[i][2]);
		}
	}
}

//---------------------------------------------------------------------
void CalculateLeLib(fiber* theFiberPtr) {

	unsigned int				i, k;

	//-----------------------------------------
	for (i = 0; i + 1 < theFiberPtr->n; i++) {
		theFiberPtr->ltm[i] = 0;
		for (k = 0; k < 3; k++)
			theFiberPtr->e[i][k] = theFiberPtr->xtm[i + 1][k] - theFiberPtr->xtm[i][k];
		theFiberPtr->ltm[i] = Norm(theFiberPtr->e[i]);
	}

	for (i = 0; i + 1 < theFiberPtr->n ; i++) {
		theFiberPtr->lt[i] = 0;
		for (k = 0; k < 3; k++)
			theFiberPtr->e[i][k] = theFiberPtr->xt[i + 1][k] - theFiberPtr->xt[i][k];
		theFiberPtr->lt[i] = Normalize(theFiberPtr->e[i]);
	}
}

//----------------------------------------------------------------------
double KappaLib(struct fiber theFiber, unsigned int i, int direction) {

	unsigned int	k;
	double			dT[3] = { 0,0,0 }, m[3] = { 0,0,0 }, vect[3] = { 0,0,0 }, vect2[3] = { 0,0,0 };

	//--------------------------------------------
	if ((direction != 1) && (direction != 2)) {
		printf("mauvaise direction dans Omega\n");
		return 0.;
	}

	if ((i >= 1) && (i <= theFiber.n - 2)) {
		for (k = 0; k < 3; k++)
			dT[k] = (theFiber.xt[i + 1][k] - 2. * theFiber.xt[i][k] + theFiber.xt[i - 1][k]) / (theFiber.l0* theFiber.l0);			
		if (direction == 1) {
			for (k = 0; k < 3; k++)
				m[k] = (theFiber.m1[i - 1][k] + theFiber.m1[i][k]) / 2.;
		}
		else {
			CrossProduct(theFiber.e[i-1], theFiber.m1[i-1], vect);
			CrossProduct(theFiber.e[i], theFiber.m1[i], vect2);
			for (k = 0; k < 3; k++)
				m[k] = (vect[k] + vect2[k]) / 2.;
		}
		return  DotProduct(dT, m);
	}
	else
		return 0.;
}

//-----------------------------------------------------------
double ZMaxLib(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int			iFiber, i;
	float					zMax = -1.e10;

	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			for (i = 0; i < theFiberPtr[iFiber].n; i++) {
				zMax = DMax(zMax, theFiberPtr[iFiber].xt[i][2]);
			}
		}
	}
	return zMax;
}

//-----------------------------------------------------------
double ZMinLib(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int			iFiber, i;
	float					zMmin = 1.e10;

	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			for (i = 0; i < theFiberPtr[iFiber].n; i++) {
				zMmin = DMin(zMmin, theFiberPtr[iFiber].xt[i][2]);
			}
		}
	}
	return zMmin;
}