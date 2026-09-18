#include	"main.h"

extern FILE* g_FilePtr;
extern parameter* gParameterPtr;
extern fiber* gFiber;
extern contact* gContact;

//---------------------------------------------------------------------
void SetFibersParametersElastica(void) {

	unsigned int				iFiber, missingSegment, n;

	printf("SetFibersParametersElastica started\n");

	//--------------- la fibre
	iFiber = 0;
	n = 60;
	AllocateOneFiberLib(gFiber, iFiber, n);
	gFiber[iFiber].status = STATUS_FREE;
	gFiber[iFiber].radius = 1.;
	gFiber[iFiber].l0 = 1.;
	gFiber[iFiber].k0 = 1.;
	gFiber[iFiber].masse = 1.;
	gFiber[iFiber].bending = gParameterPtr->bending;
	gFiber[iFiber].bendingDamping = gParameterPtr->bendingDamping;
	gFiber[iFiber].c = 1. * gParameterPtr->bending;
	gFiber[iFiber].j = gFiber[iFiber].masse;

	iFiber = 1;
	n = 2;
	AllocateOneFiberLib(gFiber, iFiber, n);
	gFiber[iFiber].status = STATUS_FIXED;
	gFiber[iFiber].radius = 1.;
	gFiber[iFiber].l0 = 10.;
	gFiber[iFiber].k0 = 1.;
	gFiber[iFiber].masse = 1.;
	gFiber[iFiber].bending = gParameterPtr->bending;
	gFiber[iFiber].bendingDamping = gParameterPtr->bendingDamping;
	gFiber[iFiber].c = 1. * gParameterPtr->bending;
	gFiber[iFiber].j = gFiber[iFiber].masse;

	iFiber = 2;
	n = 2;
	AllocateOneFiberLib(gFiber, iFiber, n);
	gFiber[iFiber].status = STATUS_FIXED;
	gFiber[iFiber].radius = 1.;
	gFiber[iFiber].l0 = 10.;
	gFiber[iFiber].k0 = 1.;
	gFiber[iFiber].masse = 1.;
	gFiber[iFiber].bending = gParameterPtr->bending;
	gFiber[iFiber].bendingDamping = gParameterPtr->bendingDamping;
	gFiber[iFiber].c = 1. * gParameterPtr->bending;
	gFiber[iFiber].j = gFiber[iFiber].masse;

	//---------- compte les segments
	gParameterPtr->nSegment = 0;
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++)
		gParameterPtr->nSegment += gFiber[iFiber].n;
	printf("nSegment avant fibre virtuelle gParameterPtr->nSegment=%u\n", gParameterPtr->nSegment);

	//----------- rajoute une fibre virtuelle pour completer à un multiple de 2 * WG
	missingSegment = ((gParameterPtr->nSegment + 2 * WG - 1) / (2 * WG)) * 2 * WG - gParameterPtr->nSegment;
	printf("Missing segment is %d\n", missingSegment);
	iFiber = gParameterPtr->nFiber - 1;
	AllocateOneFiberLib(gFiber, iFiber, missingSegment);
	gFiber[iFiber].status = STATUS_VIRTUAL;

	//---------- recompte les segments
	gParameterPtr->nSegment = 0;
	for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++)
		gParameterPtr->nSegment += gFiber[iFiber].n;
	printf("nSegment apres fibre virtuelle gParameterPtr->nSegment=%u\n", gParameterPtr->nSegment);

	printf("SetFibersParametersElastica finished\n");
}

//---------------------------------------------------------------------
void SetFibersInitialPositionsElastica(void) {

	int				i, j, k, iFiber = 0, isc, i1;
	double			sc, lx, ly, l, deltac, iDouble, alpha, beta, s, s1;

	printf("SetFibersInitialPositionsElastica started\n");
	
	//--------------- fiber 0
	iFiber = 0;
	for (i = 0; i < gFiber[iFiber].n; i++) {
		gFiber[iFiber].xt[i][0] = (double)(i) - 30.;
		gFiber[iFiber].xt[i][1] = 0.;
		gFiber[iFiber].xt[i][2] = 2.01;
	}
	iFiber = 1;
	for (i = 0; i < gFiber[iFiber].n; i++) {
		gFiber[iFiber].xt[i][0] = -15.;
		gFiber[iFiber].xt[i][1] = (((double)i) - 0.5) * gFiber[iFiber].l0;
		gFiber[iFiber].xt[i][2] = 0.;
	}
	iFiber = 2;
	for (i = 0; i < gFiber[iFiber].n; i++) {
		gFiber[iFiber].xt[i][0] = 15.;
		gFiber[iFiber].xt[i][1] = (((double)i) - 0.5) * gFiber[iFiber].l0;
		gFiber[iFiber].xt[i][2] = 0.;
	}
	for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].m1[i][0] = 1.;
			gFiber[iFiber].m1[i][1] = 0.;
			gFiber[iFiber].m1[i][2] = 0.;
			gFiber[iFiber].m1_bar[i][0] = 0.;
			gFiber[iFiber].m1_bar[i][1] = 1.;
			gFiber[iFiber].m1_bar[i][2] = 0.;
			gFiber[iFiber].kappa1_bar[i] = 0.;
			gFiber[iFiber].kappa2_bar[i] = 0.;
		}
	}

	CalculateLeLib(gFiber, gParameterPtr);

	//------------- au repos -----------------------------
	for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			for (j = 0; j < 3; j++) {
				gFiber[iFiber].thetat[i][j] = 0.;
				gFiber[iFiber].xtm[i][j] = gFiber[iFiber].xt[i][j];
				gFiber[iFiber].thetatm[i][j] = gFiber[iFiber].thetat[i][j];
			}
		}
	}
}