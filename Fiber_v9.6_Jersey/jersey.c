#include	"main.h"

extern FILE				* g_FilePtr;
extern parameter		* gParameterPtr;
extern fiber			* gFiber;
extern contact			* gContact;

//---------------------------------------------------------------------
void SetFibersParameters(void) {

	unsigned int				iFiber, missingSegment;

	printf("SetFibersParameters started\n");

	//---------------------- traversantes horizontale
	for (iFiber = 0; iFiber < GN_CELL_Y - 1; iFiber++)
		AllocateOneFiberLib(gFiber, iFiber, GN_CELL_X * NLOOP + 1);
	printf("aa\n");
	//---------------------- coin haut/gauche
	AllocateOneFiberLib(gFiber, GN_CELL_Y-1, NLOOP / 4 + 1);

	//---------------------- fibres basses ---------
	for (iFiber = GN_CELL_Y; iFiber < GN_CELL_X + GN_CELL_Y; iFiber++)
		AllocateOneFiberLib(gFiber, iFiber, NLOOP / 2 + 1);

	//---------------------- fibres hautes ---------
	for (iFiber = GN_CELL_X + GN_CELL_Y; iFiber <= 2 * GN_CELL_X + GN_CELL_Y - 2; iFiber++)
		AllocateOneFiberLib(gFiber, iFiber, NLOOP / 2 + 1);

	//---------------------- coin haut/droit
	AllocateOneFiberLib(gFiber, 2 * GN_CELL_X + GN_CELL_Y -1, NLOOP / 4 + 1);



	//-------------------------------------------------
	for (iFiber = 0; iFiber <= 2 * GN_CELL_X + GN_CELL_Y; iFiber++) {
		gFiber[iFiber].status = STATUS_FREE;
		gFiber[iFiber].radius = gParameterPtr->rTarget;
		gFiber[iFiber].l0 = 1.;
		gFiber[iFiber].k0 = 1.;
		gFiber[iFiber].masse = 1.;
		gFiber[iFiber].bending = gParameterPtr->bending * gParameterPtr->rTarget * gParameterPtr->rTarget / 4.;
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

//---------------------------------------------------------------------
void SetFibersInitialPositions(void) {

	int			iFiber, i, index, is, k;
	double		** loop_raw, ** loop, shift[3], s, ds;
	FILE		* filePtr;

	printf("SetFibersInitialPositions started\n");

	//------------ lecture de la loop raw
	filePtr = fopen("x_theta_lw70_lc70.txt", "r");
	if (!filePtr)
		printf("opening file failed in SetFiberJersey\n");

	loop_raw = (double**)malloc(NLOOP_RAW * sizeof(double*));
	for (i = 0; i < NLOOP_RAW; i++)
		loop_raw[i] = (double*)malloc(3 * sizeof(double));

	for (i = 0; i < NLOOP_RAW; i++) {
		fscanf(filePtr, "%lf %lf %lf", &(loop_raw[i][0]), &(loop_raw[i][1]), &(loop_raw[i][2]));
		loop_raw[i][0] += LC_LW_RAW / 2;
		loop_raw[i][1] = LC_LW_RAW - loop_raw[i][1] - LC_LW_RAW / 2.;
	}
	fclose(filePtr);
	printf("Loading loop ended\n");


	//------------- fabrication de la loop renormalisée
	loop = (double**)malloc(NLOOP * sizeof(double*));
	for (i = 0; i < NLOOP; i++)
		loop[i] = (double*)malloc(3 * sizeof(double));

	for (i = 0; i < NLOOP; i++) {
		s = 0.99999 * ((double)(i)) / SCALE;
		is = IMin((int)s, NLOOP_RAW - 1);
		ds = s - (double)is;
		for (k = 0; k < 3; k++) {
			loop[i][k] = (1. - ds) * loop_raw[is][k] + ds * loop_raw[is + 1][k];
			loop[i][k] *= SCALE;
		}
	}
	printf("Rescaling loop ended\n");

	//----------------------- fibres traversantes
	for (iFiber = 0; iFiber < GN_CELL_Y - 1; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].xt[i][0] = +((double)(i / NLOOP)) * LC_LW + loop[i % NLOOP][0];
			gFiber[iFiber].xt[i][1] = +((double)iFiber) * LC_LW + loop[i % NLOOP][1];
			gFiber[iFiber].xt[i][2] = +loop[i % NLOOP][2];
		}
	}
	printf("Setting fiber step #1 ended\n");

	//------------------ coin haut gauche ---------
	iFiber = GN_CELL_Y-1;
	for (i = 0; i < gFiber[iFiber].n; i++) {
		gFiber[iFiber].xt[i][0] = loop[i % NLOOP][0];
		gFiber[iFiber].xt[i][1] = +((double)GN_CELL_Y - 1) * LC_LW + loop[i % NLOOP][1];
		gFiber[iFiber].xt[i][2] = +loop[i % NLOOP][2];
	}
	printf("Setting fiber step #2 ended\n");

	//------------------ fibres basses ------------
	for (iFiber = GN_CELL_Y; iFiber < GN_CELL_X + GN_CELL_Y; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].xt[i][0] = ((double)(iFiber - GN_CELL_Y)) * LC_LW + loop[i + NLOOP / 4][0];
			gFiber[iFiber].xt[i][1] = -LC_LW + loop[i + NLOOP / 4][1];
			gFiber[iFiber].xt[i][2] = loop[i + NLOOP / 4][2];
		}
	}
	printf("Setting fiber step #3 ended\n");

	//------------------ fibres hautes --------
	for (iFiber = GN_CELL_X + GN_CELL_Y; iFiber <= 2 * GN_CELL_X + GN_CELL_Y - 2; iFiber++) {
		for (i = 0; i < NLOOP / 4; i++) {
			index = i + 3 * NLOOP / 4;
			gFiber[iFiber].xt[i][0] = LC_LW * ((double)(iFiber - (GN_CELL_X + GN_CELL_Y))) + loop[index][0];
			gFiber[iFiber].xt[i][1] = LC_LW * ((double)GN_CELL_Y - 1) + loop[index][1];
			gFiber[iFiber].xt[i][2] = loop[index][2];
		}
		for (; i < gFiber[iFiber].n; i++) {
			index = i + -NLOOP / 4;
			gFiber[iFiber].xt[i][0] = LC_LW * ((double)(iFiber - (GN_CELL_X + GN_CELL_Y) + 1)) + loop[index][0];
			gFiber[iFiber].xt[i][1] = LC_LW * ((double)GN_CELL_Y - 1) + loop[index][1];
			gFiber[iFiber].xt[i][2] = loop[index][2];
		}
	}
	printf("Setting fiber step #4 ended\n");

	//------------------ coin haut droit ---------
	iFiber = 2 * GN_CELL_X + GN_CELL_Y-1;
	for (i = 0; i < NLOOP / 4; i++) {
		index = i + 3 * NLOOP / 4;
		gFiber[iFiber].xt[i][0] = LC_LW * ((double)(iFiber - (GN_CELL_X + GN_CELL_Y))) + loop[index][0];
		gFiber[iFiber].xt[i][1] = LC_LW * ((double)GN_CELL_Y - 1) + loop[index][1];
		gFiber[iFiber].xt[i][2] = loop[index][2];
	}
	for (; i < gFiber[iFiber].n; i++) {
		index = i - NLOOP / 4;
		gFiber[iFiber].xt[i][0] = LC_LW * ((double)(iFiber - (GN_CELL_X + GN_CELL_Y) + 1)) + loop[index][0];
		gFiber[iFiber].xt[i][1] = LC_LW * ((double)GN_CELL_Y - 1) + loop[index][1];
		gFiber[iFiber].xt[i][2] = loop[index][2];
	}
	printf("Setting fiber step #5 ended\n");

	//----------------------recentrage
	shift[0] = -LC_LW * ((double)(GN_CELL_X) / 2);
	shift[1] = -LC_LW * ((double)(GN_CELL_Y) / 2);
	shift[2] = 0;

	for (iFiber = 0; iFiber < gParameterPtr->nFiber -1; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			for (k = 0; k < 3; k++)
				gFiber[iFiber].xt[i][k] += shift[k];
		}
	}
	printf("recentering ended\n");

	//-------------------- theta, kappa
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			for (k = 0; k < 3; k++)
				gFiber[iFiber].thetat[i][k] = 0.;
			gFiber[iFiber].kappa1_bar[i] = 0.;
			gFiber[iFiber].kappa2_bar[i] = 0.;
		}
	}
	printf("set theta kappa ended\n");

	//-------------------- m1_bar, m1
	for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
		CalculateLeLib(&(gFiber[iFiber]));
		gFiber[iFiber].m1_bar[0][0] = 0.;
		gFiber[iFiber].m1_bar[0][1] = 0.;
		gFiber[iFiber].m1_bar[0][2] = 1.;
		Calculate_Relaxed_M1_Bar(&(gFiber[iFiber]));
		Calculate_M1_From_M1Bar_And_Theta(&(gFiber[iFiber]));
	}
	printf("set m1_bar m1 kappa ended\n");

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