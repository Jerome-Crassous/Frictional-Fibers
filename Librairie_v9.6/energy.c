#include <stdlib.h>
#include <stdio.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

void ForceBendingOneNodeLib(fiber theFiber, int iNode, double* forceLocal);

//---------------------------------------------------------------------
double KineticEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int		iFiber, i, k;
	double				energy = 0;

	//-------------------------------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			for (i = 0; i < theFiberPtr[iFiber].n; i++) {
				for (k = 0; k < 3; k++)
					energy += theFiberPtr[iFiber].masse
					* (theFiberPtr[iFiber].xt[i][k] - theFiberPtr[iFiber].xtm[i][k])
					* (theFiberPtr[iFiber].xt[i][k] - theFiberPtr[iFiber].xtm[i][k]);
			}
			for (i = 0; i < theFiberPtr[iFiber].n - 1; i++)
				energy += theFiberPtr[iFiber].j
				* (theFiberPtr[iFiber].thetat[i][2] - theFiberPtr[iFiber].thetatm[i][2])
				* (theFiberPtr[iFiber].thetat[i][2] - theFiberPtr[iFiber].thetatm[i][2]);/**/
		}
	}
	return energy / (2. * theParameterPtr->dt * theParameterPtr->dt);
}

//---------------------------------------------------------------------
double BendingEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int		i, iFiber;
	double				aux = 0, kappa, energy = 0;

	//-------------------------------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			aux = 0;
			for (i = 1; i < theFiberPtr[iFiber].n - 1; i++) {
				kappa = KappaLib(theFiberPtr[iFiber], i, 1) - theFiberPtr[iFiber].kappa1_bar[i];
				aux += kappa * kappa;
				kappa = KappaLib(theFiberPtr[iFiber], i, 2) - theFiberPtr[iFiber].kappa2_bar[i];
				aux += kappa * kappa;
			}
			aux *= theFiberPtr[iFiber].bending * theFiberPtr[iFiber].l0 / 2.;
			energy += aux;
		}
	}
	return energy;
}

//---------------------------------------------------------------------
double TwistEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int		iFiber, i;
	double				aux, beta, energy = 0, vect[3];

	//-------------------------------------------------------------------
	for (iFiber = 0; iFiber + 1 < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			aux = 0;
			for (i = 1; i +2 < theFiberPtr[iFiber].n; i++) {
				CrossProduct(theFiberPtr[iFiber].m1[i], theFiberPtr[iFiber].m1[i - 1], vect);
				beta = DotProduct(theFiberPtr[iFiber].e[i - 1], vect) + DotProduct(theFiberPtr[iFiber].e[i], vect);
				aux += beta * beta;
			}
			aux *= theFiberPtr[iFiber].c / 8. / theFiberPtr[iFiber].l0;
			energy += aux;
		}
	}
	return energy;
}

//---------------------------------------------------------------------
double TractionEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int		i, iFiber;
	double				energy = 0;

	//-------------------------------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			for (i = 0; i < theFiberPtr[iFiber].n - 1; i++) {
				energy += theFiberPtr[iFiber].k0
					* (theFiberPtr[iFiber].lt[i] - theFiberPtr[iFiber].l0)
					* (theFiberPtr[iFiber].lt[i] - theFiberPtr[iFiber].l0)
					/ 2. / theFiberPtr[iFiber].l0;
			}
		}
	}
	return energy;
}

//-----------------------------------------------------------
double Elastic_Normal_Contact_Energy(contact* theContactPtr, parameter* theParameterPtr) {

	unsigned int  		iContact;
	double				energy = 0;

	//--------------------------------------------------------------------	
	for (iContact = 0; iContact < theParameterPtr->nContact; iContact++)
		energy += theContactPtr[iContact].deltatm * theContactPtr[iContact].deltatm;
	return energy * theParameterPtr->kn / 2.;
}

//-----------------------------------------------------------
double Elastic_Tangential_Contact_Energy(contact* theContactPtr, parameter* theParameterPtr) {

	unsigned int  		iContact;
	double				energy = 0;

	//--------------------------------------------------------------------	
	for (iContact = 0; iContact < theParameterPtr->nContact; iContact++)
		energy += DotProduct(theContactPtr[iContact].u_t, theContactPtr[iContact].u_t);
	return energy * theParameterPtr->kt / 2.;
}

//---------------------------------------------------------------------
void InternalForcesLib(fiber theFiber, unsigned int i, int direction, double* force) {

	unsigned int	k;
	double			force1[3], force2[3];

	InternalTractionForcesLib(theFiber, i, direction, force1);
	InternalBendingForcesLib(theFiber, i, direction, force2);

	for (k = 0; k < 3; k++)
		force[k] = force1[k] + force2[k];

	return;
}

//---------------------------------------------------------------------
void InternalTractionForcesLib(fiber theFiber, unsigned int i, int direction, double* force) {

	unsigned int				k;

	for (k = 0; k < 3; k++) force[k] = 0;

	//---------------------------------------
	if ((direction == 1) && (i < (theFiber.n - 1))) {
		for (k = 0; k < 3; k++)
			force[k] += theFiber.e[i][k] * theFiber.k0 * (theFiber.lt[i] - theFiber.l0) / theFiber.l0;
	}
	else if ((direction == -1) && (i > 0)) {
		for (k = 0; k < 3; k++)
			force[k] -= theFiber.e[i - 1][k] * theFiber.k0 * (theFiber.lt[i - 1] - theFiber.l0) / theFiber.l0;
	}
	return;
}

//---------------------------------------------------------------------
void ForceBendingOneNodeLib(fiber theFiber, int iNode, double* forceLocal) {

	unsigned int				k;
	double						kappa_red, m1[3], m2[3], t[3];

	//-----------------------------
	for (k = 0; k < 3; k++) forceLocal[k] = 0.;

	//-----------------------------
	if ((iNode < 1) || (iNode + 1 > theFiber.n)) return;

	for (k = 0; k < 3; k++) {
		m1[k] = (theFiber.m1[iNode - 1][k] + theFiber.m1[iNode][k]) / 2.;
		t[k] = (theFiber.e[iNode - 1][k] + theFiber.e[iNode][k]) / 2.;
	}
	CrossProduct(t, m1, m2);

	kappa_red = KappaLib(theFiber, iNode, 1) - theFiber.kappa1_bar[iNode];
	for (k = 0; k < 3; k++) forceLocal[k] += kappa_red * m1[k];
	kappa_red = KappaLib(theFiber, iNode, 2) - theFiber.kappa2_bar[iNode];
	for (k = 0; k < 3; k++) forceLocal[k] += kappa_red * m2[k];
}

//---------------------------------------------------------------------
void InternalBendingForcesLib(fiber theFiber, unsigned int i, int direction, double* force) {

	unsigned int				k;
	double						forceLocal[3];

	for (k = 0; k < 3; k++) force[k] = 0.;

	if (direction == 1) {
		ForceBendingOneNodeLib(theFiber, i + 1, forceLocal);
		for (k = 0; k < 3; k++)
			force[k] -= forceLocal[k];
		ForceBendingOneNodeLib(theFiber, i, forceLocal);
		for (k = 0; k < 3; k++)
			force[k] += forceLocal[k];
	}
	else if (direction == -1) {
		ForceBendingOneNodeLib(theFiber, i, forceLocal);
		for (k = 0; k < 3; k++)
			force[k] += forceLocal[k];
		ForceBendingOneNodeLib(theFiber, i - 1, forceLocal);
		for (k = 0; k < 3; k++)
			force[k] -= forceLocal[k];
	}

	for (k = 0; k < 3; k++)
		force[k] *= theFiber.bending / theFiber.l0;
}

//---------------------------------------------------------------------
void InternalBendingMomentLib(fiber theFiber, unsigned int i, int direction, double* moment) {

	unsigned int				k;
	double						kappa1_red, kappa2_red, m1[3], m2[3], t[3];

	//-----------------------------
	for (k = 0; k < 3; k++) moment[k] = 0.;

	if ((i == 0) || ((i + 1) == theFiber.n) || (theFiber.status != STATUS_FREE))
		return;

	//-------------------------------------------------------------------
	kappa1_red = KappaLib(theFiber, i, 1) - theFiber.kappa1_bar[i];
	kappa2_red = KappaLib(theFiber, i, 2) - theFiber.kappa2_bar[i];

	printf("kappa1_red %e kappa2_red %e\n", kappa1_red, kappa2_red);
	for (k = 0;k < 3;k++) {
		m1[k] = (theFiber.m1[i - 1][k] + theFiber.m1[i][k]) / 2.;
		t[k] = (theFiber.e[i - 1][k] + theFiber.e[i][k]) / 2.;
	}
	CrossProduct(t, m1, m2);

	for (k = 0;k < 3;k++) {
		moment[k] = (kappa1_red * m2[k] - kappa2_red * m1[k]);
		if (direction == 1)
			moment[k] *= theFiber.bending;
		else if (direction == -1)
			moment[k] *= -theFiber.bending;
	}
}