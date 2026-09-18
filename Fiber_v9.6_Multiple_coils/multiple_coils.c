#include	"main.h"

extern FILE				* g_FilePtr;
extern parameter		* gParameterPtr;
extern fiber			* gFiber;
extern contact			* gContact;

void S2PositionHelix(double t, double* x);
void SetOneFiberInitialPositions(unsigned int iFiber);

//---------------------------------------------------------------------
void SetFibersParameters(void) {

	unsigned int				iFiber, missingSegment;

	printf("SetFibersParameters started\n");

	//---------------------- traversantes horizontale
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		AllocateOneFiberLib(gFiber, iFiber, 48);
		gFiber[iFiber].status = STATUS_FREE;
		gFiber[iFiber].radius = gParameterPtr->r;
		gFiber[iFiber].l0 = 1.;
		gFiber[iFiber].k0 = 1.;
		gFiber[iFiber].masse = 1.;
		gFiber[iFiber].bending = 1.;
		gFiber[iFiber].c = gFiber[iFiber].bending / 1.5;
		gFiber[iFiber].j = gFiber[iFiber].masse;
	}

	//---------- compte les segments
	gParameterPtr->nSegment = 0;
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++)
		gParameterPtr->nSegment += gFiber[iFiber].n;
	printf("nSegment avant fibre virtuelle gParameterPtr->nSegment=%u\n", gParameterPtr->nSegment);

	//----------- rajoute une fibre virtuelle pour completer à un multiple de WG
	missingSegment = ((gParameterPtr->nSegment + WG - 1) / (WG)) * WG - gParameterPtr->nSegment;
	printf("Missing segment is %d\n", missingSegment);
	iFiber = gParameterPtr->nFiber - 1;
	AllocateOneFiberLib(gFiber, iFiber, missingSegment);
	gFiber[iFiber].status = STATUS_VIRTUAL;

	//---------- recompte les segments
	gParameterPtr->nSegment = 0;
	for (iFiber = 0; iFiber < gParameterPtr->nFiber; iFiber++)
		gParameterPtr->nSegment += gFiber[iFiber].n;
	printf("nSegment apres fibre virtuelle gParameterPtr->nSegment=%u\n", gParameterPtr->nSegment);

	printf("SetFibersParameter finished\n");
}

//-----------------------------------------------------------
void SetFibersInitialPositions(void) {

	unsigned int				iFiber;
	
	printf("SetFibersInitialPositions started\n");
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++)
		SetOneFiberInitialPositions(iFiber);

	printf("SetFibersInitialPositions finished\n");
}

//-----------------------------------------------------------
void SetOneFiberInitialPositions(unsigned int iFiber) {

	unsigned int			i, k;
	double					t = 0, dt = 1.e-4, dT[3], m1[3], m2[3];

	//-------------- positions
	S2PositionHelix(t, gFiber[iFiber].xt[0]);
	for (i = 1; i < gFiber[iFiber].n; i++) {
		for (;;) {
			t += dt;
			S2PositionHelix(t, gFiber[iFiber].xt[i]);
			if (Distance(gFiber[iFiber].xt[i - 1], gFiber[iFiber].xt[i]) > gFiber[iFiber].l0)
				break;
		}
	}

	//------------- theta
	for (i = 0; i < gFiber[iFiber].n; i++) {
		gFiber[iFiber].thetat[i][0] = 0.;
		gFiber[iFiber].thetat[i][1] = 0.;
		gFiber[iFiber].thetat[i][2] = 0.;
	}

	//---------------- m1 et m1_bar
	CalculateLeLib(&(gFiber[iFiber]));
	gFiber[iFiber].m1_bar[0][0] = 0.;
	gFiber[iFiber].m1_bar[0][1] = 0.5;
	gFiber[iFiber].m1_bar[0][2] = 0.5;
	Calculate_Relaxed_M1_Bar(&(gFiber[iFiber]));
	Calculate_M1_From_M1Bar_And_Theta(&(gFiber[iFiber]));

	//----------------- kappa
	gFiber[iFiber].kappa1_bar[0] = 0.;
	gFiber[iFiber].kappa2_bar[0] = 0.;
	gFiber[iFiber].kappa1_bar[gFiber[iFiber].n - 1] = 0.;
	gFiber[iFiber].kappa2_bar[gFiber[iFiber].n - 1] = 0.;
	for (i = 1; i < gFiber[iFiber].n - 1; i++) {
		for (k = 0; k < 3; k++) {
			dT[k] = (gFiber[iFiber].e[i][k] - gFiber[iFiber].e[i - 1][k]) / gFiber[iFiber].l0;
			m1[k] = (gFiber[iFiber].m1_bar[i][k] + gFiber[iFiber].m1_bar[i - 1][k]) / 2.;
		}
		CrossProduct(gFiber[iFiber].e[i], m1, m2);
		gFiber[iFiber].kappa1_bar[i] = DotProduct(m1, dT);
		gFiber[iFiber].kappa2_bar[i] = DotProduct(m2, dT);
	}

	//-------------- au repos
	for (i = 0; i < gFiber[iFiber].n; i++) {
		for (k = 0; k < 3; k++) {
			gFiber[iFiber].xtm[i][k] = gFiber[iFiber].xt[i][k];
			gFiber[iFiber].thetatm[i][k] = gFiber[iFiber].thetat[i][k];
		}
	}

	//-------------- shifté
	double shift[3];
	shift[0] = .8 * (MyRand() - 0.5) * gParameterPtr->lx;
	shift[1] = .8 * (MyRand() - 0.5) * gParameterPtr->ly;
	shift[2] = 0.7 * (MyRand() - 0.5) * gParameterPtr->lz;

	for (i = 0; i < gFiber[iFiber].n; i++) {
		for (k = 0; k < 3; k++) {
			gFiber[iFiber].xt[i][k] += shift[k];
			gFiber[iFiber].xtm[i][k] += shift[k];
		}
	}
}

//------------------------------------------------------
void S2PositionHelix(double t, double* x) {

	x[0] = gParameterPtr->RHelix * sin(t);
	x[1] = gParameterPtr->RHelix * cos(t);
	x[2] = gParameterPtr->pitch * t;
}