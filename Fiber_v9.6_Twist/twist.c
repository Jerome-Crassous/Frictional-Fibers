#include	"main.h"

extern FILE				* g_FilePtr;
extern parameter		* gParameterPtr;
extern fiber			* gFiber;
extern contact			* gContact;

extern cl_command_queue command_queue;
extern cl_double		*param_Double_hostPtr;
extern cl_mem			param_Uint_devPtr, param_Double_devPtr;

//---------------------------------------------------------------------
void SetFibersParameters(void) {

	unsigned int				iFiber, missingSegment;

	printf("SetFibersParameters started\n");

	//--------------------------------------------------------
	iFiber = 0;
	AllocateOneFiberLib(gFiber, iFiber, 20);
	gFiber[iFiber].status = STATUS_FREE;
	gFiber[iFiber].radius = 1.;
	gFiber[iFiber].l0 = gParameterPtr->ls;
	gFiber[iFiber].k0 = 1.;
	gFiber[iFiber].masse = 1.;
	gFiber[iFiber].bending = gParameterPtr->bending;
	gFiber[iFiber].c = gParameterPtr->c;
	gFiber[iFiber].j = gFiber[iFiber].masse;

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

//---------------------------------------------------------------------
void SetFibersInitialPositions(void) {

	int				i, k, iFiber, n;
	double			theta, x0, y0, z0, aux1, aux2, sMiddle, ll, bis, z1, z2;
	double					T[3], dT[3], Darboux[3], aux, m1[3], m2[3];

	printf("SetFibersInitialPositions started\n");

	//------------- fibre 0
	iFiber = 0;
	
	for (i = 0; i < gFiber[0].n; i++) {
		gFiber[iFiber].xt[i][0] = 0.;
		gFiber[iFiber].xt[i][1] = 0.;
		gFiber[iFiber].xt[i][2] = -((double)(i)) * gFiber[0].l0;

		gFiber[iFiber].thetat[i][0] = 0.; 
		gFiber[iFiber].thetat[i][1] = 0.; 
		gFiber[iFiber].thetat[i][2] = 0.;

		gFiber[iFiber].kappa1_bar[i] = 0.;
		gFiber[iFiber].kappa2_bar[i] = 0.;
	}
	CalculateLeLib(&(gFiber[iFiber]));
	gFiber[iFiber].m1_bar[0][0] = 0.;
	gFiber[iFiber].m1_bar[0][1] = 1.;
	gFiber[iFiber].m1_bar[0][2] = 0.;
	Calculate_Relaxed_M1_Bar(&(gFiber[iFiber]));
	Calculate_M1_From_M1Bar_And_Theta(&(gFiber[iFiber]));

	//-------------- au repos
	for(iFiber=0;iFiber<gParameterPtr->nFiber-1;iFiber++){
		for (i = 0; i < gFiber[iFiber].n; i++) {
			for (k = 0; k < 3; k++) {
				gFiber[iFiber].xtm[i][k] = gFiber[iFiber].xt[i][k];
				gFiber[iFiber].thetatm[i][k] = gFiber[iFiber].thetat[i][k];
			}
		}
	}
	printf("SetFibersInitialPositions finished\n");
}