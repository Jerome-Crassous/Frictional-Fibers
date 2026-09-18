#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

//---------------------------------------------------------------------
void AllocateOneFiberLib(fiber* theFiberPtr, unsigned int iFiber, unsigned int npt) {

	unsigned int				i, failed = 0;

	//-------------------------------------------------------
	theFiberPtr[iFiber].n = npt;

	//--------------------------------------------------------
	theFiberPtr[iFiber].xtm = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].xtm == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].xtm[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].xtm[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].xt = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].xt == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].xt[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].xt[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].thetatm = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].thetatm == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].thetatm[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].thetatm[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].thetat = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].thetat == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].thetat[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].thetat[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].f = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].f == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].f[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].f[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].e = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].e == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].e[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].e[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].moment = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].moment == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].moment[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].moment[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].m1_bar = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].m1_bar == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].m1_bar[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].m1_bar[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].m1 = (double**)GlobalAlloc(0, theFiberPtr[iFiber].n * sizeof(double*));
	if (theFiberPtr[iFiber].m1 == NULL)
		failed = 1;
	else {
		for (i = 0; i < theFiberPtr[iFiber].n; i++) {
			theFiberPtr[iFiber].m1[i] = (double*)GlobalAlloc(GPTR, 3 * sizeof(double));
			if (theFiberPtr[iFiber].m1[i] == NULL)
				failed = 1;
		}
	}

	//--------------------------------------------------------
	theFiberPtr[iFiber].flag = (unsigned int*)GlobalAlloc(GPTR, theFiberPtr[iFiber].n * sizeof(unsigned int));
	if (theFiberPtr[iFiber].flag == NULL) failed = 1;
	theFiberPtr[iFiber].ltm = (double*)GlobalAlloc(GPTR, theFiberPtr[iFiber].n * sizeof(double));
	if (theFiberPtr[iFiber].ltm == NULL) failed = 1;
	theFiberPtr[iFiber].lt = (double*)GlobalAlloc(GPTR, theFiberPtr[iFiber].n * sizeof(double));
	if (theFiberPtr[iFiber].lt == NULL) failed = 1;
	theFiberPtr[iFiber].kappa1_bar = (double*)GlobalAlloc(GPTR, theFiberPtr[iFiber].n * sizeof(double));
	if (theFiberPtr[iFiber].kappa1_bar == NULL) failed = 1;
	theFiberPtr[iFiber].kappa2_bar = (double*)GlobalAlloc(GPTR, theFiberPtr[iFiber].n * sizeof(double));
	if (theFiberPtr[iFiber].kappa2_bar == NULL) failed = 1;

	if (failed != 0)
		printf("Allocation Fiber failed\n");

	theFiberPtr[iFiber].n = npt;
}