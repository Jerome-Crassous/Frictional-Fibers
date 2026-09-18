#include	"main.h"
#define L	15

extern FILE* g_FilePtr;
extern parameter* gParameterPtr;
extern fiber* gFiber;
extern contact* gContact;

void S2Position(int t, double* position);

//---------------------------------------------------------------------
void SetFibersParameters(void) {

	unsigned int				iFiber, missingSegment;

	printf("SetFibersParametersTennis started\n");

	//---------------------- traversantes horizontale
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++)
		AllocateOneFiberLib(gFiber, iFiber, 301);

	//-------------------------------------------------
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		gFiber[iFiber].status = STATUS_FREE;
		gFiber[iFiber].radius = gParameterPtr->r;
		gFiber[iFiber].l0 = 1.;
		gFiber[iFiber].k0 = 1.;
		gFiber[iFiber].masse = 1.;
		gFiber[iFiber].bending = gParameterPtr->bending;
		gFiber[iFiber].c = gParameterPtr->c;
		gFiber[iFiber].j = gFiber[iFiber].masse;
	}

	//---------- compte les segments
	gParameterPtr->nSegment = 0;
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++)
		gParameterPtr->nSegment += gFiber[iFiber].n;
	printf("nSegment avant fibre virtuelle gParameterPtr->nSegment=%u\n", gParameterPtr->nSegment);

	//----------- rajoute une fibre virtuelle pour completer à un multiple de WG
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

	printf("SetFibersParametersTennis finished\n");
}

//-----------------------------------------------------------
void SetFibersInitialPositions(void) {

	int						iFiber, i, k;
	double					ez[3] = { 0.,0.,1. };
	double					d=10., shift;
	int						nHalf, nHalfHalf, iMiddle;
	double					t, dt = 1.e-3, position[3];

	printf("SetFibersInitialPositions start\n");
	printf("zzzzzzzzzzz %d\n", gParameterPtr->nFiber);
	nHalf = (gParameterPtr->nFiber - 1) / 2;
	nHalfHalf = (nHalf - 1) / 2;
	printf("nHalf %d - nHalffHalf %d\n", nHalf, nHalfHalf);
	
	//--------------- positions
	for (iFiber = 0; iFiber < nHalf; iFiber++) {
		printf("%d\n", iFiber);
		iMiddle = (gFiber[iFiber].n - 1) / 2.;
		t = 0.;

		//----------------- position
		for (i = iMiddle; i < gFiber[iFiber].n - 1; i++) {
			S2Position(t, gFiber[iFiber].xt[i]);
			for (;;) {
				t += dt;
				S2Position(t, gFiber[iFiber].xt[i + 1]);
				if (Distance(gFiber[iFiber].xt[i], gFiber[iFiber].xt[i + 1]) > 1.)
					break;
			}
		}

		//----------------- symetrique
		for (i = iMiddle - 1; i >= 0; i--) {
			gFiber[iFiber].xt[i][0] = -gFiber[iFiber].xt[2 * iMiddle - i][0];
			gFiber[iFiber].xt[i][1] = gFiber[iFiber].xt[2 * iMiddle - i][1];
			gFiber[iFiber].xt[i][2] = gFiber[iFiber].xt[2 * iMiddle - i][2];
		}

		//----------------- symetrique
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].xt[i][1] += ((double)(iFiber - nHalfHalf)) * L;
			gFiber[iFiber].xt[i][2] *= pow(-1, (iFiber - nHalfHalf));
		}
	}
	printf("aaaaaaaaaa %d\n", gParameterPtr->nFiber);

	//--------------- positions
	for (iFiber = nHalf; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].xt[i][0] = ((double)(iFiber - nHalf - nHalfHalf)) * L;
			gFiber[iFiber].xt[i][1] = (double)(i - 150);
			gFiber[iFiber].xt[i][2] = 0.;
		}
	}

	//------------ thetat
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].thetat[i][0] = 0.;
			gFiber[iFiber].thetat[i][1] = 0.;
			gFiber[iFiber].thetat[i][2] = 0.;
		}
	}

	//---------------- m1 et m1_bar
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		CalculateLeLib(&(gFiber[iFiber]));
		gFiber[iFiber].m1_bar[0][0] = 0.;
		gFiber[iFiber].m1_bar[0][1] = 0.5;
		gFiber[iFiber].m1_bar[0][2] = 0.5;
		Calculate_Relaxed_M1_Bar(&(gFiber[iFiber]));
		Calculate_M1_From_M1Bar_And_Theta(&(gFiber[iFiber]));
	}

	//----------- au repos
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].kappa1_bar[i] = 0.;
			gFiber[iFiber].kappa2_bar[i] = 0.;
			for (k = 0; k < 3; k++) {
				gFiber[iFiber].xtm[i][k] = gFiber[iFiber].xt[i][k];
				gFiber[iFiber].thetatm[i][k] = gFiber[iFiber].thetat[i][k];
			}
		}
	}

	printf("SetFibersInitialPositions finished\n");
}

//---------------------------------------------
void S2Position(int t, double* position) {

	double				amp = 1.;

	position[0] = t;
	position[1] = 0.;
	position[2] = amp * cos(2. * PI * t / (2. * L));
}