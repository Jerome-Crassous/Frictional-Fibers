#include "main.h"

parameter* gParameterPtr;
FILE* g_FilePtr;
fiber* gFiber;
contact* gContact;

//----------------------------------------------------------------------
int main(int argc, char** argv) {

	int					i, ok;
	double				trash, param1, param2, param3, sMin, sMax;

	//------------ init globale parameter
	gParameterPtr = (parameter*)malloc(sizeof(parameter));

	gParameterPtr->fx = 0.0;
	gParameterPtr->fy = 0.0;
	gParameterPtr->fix = 0.0;
	gParameterPtr->sc = 0.05;		// sc/l		0 < sc/l < 1/4
	gParameterPtr->thetac = 0.;
	gParameterPtr->lx = 2.770833e-01;		// lx/l		sc > lx/4
	gParameterPtr->ly = 3.256906e-01;		// ly/l		sc < l/4 - ly/2		/**/

	gParameterPtr->flag0 = 2;

	gParameterPtr->iter0 = 10 * KI;		//low eta - high eta
	gParameterPtr->iter1 = 150 * MI;	//start harmonic compression 
	gParameterPtr->iter3 = 150 * KI;	//start compression
	gParameterPtr->iterStop = MI;

	gParameterPtr->iter = 0;
	gParameterPtr->flag_Reset_Contact = 1;
	gParameterPtr->periodic = PERIODIC_NO;
	gParameterPtr->periodicVisu = PERIODIC_VISU_NO; 
	gParameterPtr->line = 64;
	gParameterPtr->nFiber = 2;
	gParameterPtr->nSegment = 0;
	gParameterPtr->nContactMax = (32*WG);
	gParameterPtr->epsStar = 0.1;

	gParameterPtr->r = 0.2;
	gParameterPtr->rTarget = 2. / 30;
	gParameterPtr->R = 20.;
	gParameterPtr->RTarget = 100.;
	gParameterPtr->h = 150.;
	gParameterPtr->hTarget = 14.;
	gParameterPtr->dhdt = 4.e-5;
	gParameterPtr->RHelix = 8.;
	gParameterPtr->pitch = 0.1;

	gParameterPtr->dt = 0.1;
	gParameterPtr->kn = 1.;
	gParameterPtr->kt = 0.5;
	gParameterPtr->lambda = 1.e-3;
	gParameterPtr->lambda_internal = 2.8;
	gParameterPtr->lambda_contact_n = 1.;
	gParameterPtr->lambda_contact_t = 0.;
	gParameterPtr->bending = 1.;
	gParameterPtr->c = gParameterPtr->bending;
	gParameterPtr->mu = 0.4;
	gParameterPtr->muTarget = gParameterPtr->mu;

	//------------- int globales file
	g_FilePtr = fopen("logFile.tmp", "w");

	//-------------- init. des fibres et billes
	printf("Initialisation fibre start\n");
	gFiber = (fiber*)malloc((gParameterPtr->nFiber) * sizeof(fiber));
	gContact = (contact*)malloc((gParameterPtr->nContactMax) * sizeof(contact));

	SetFibersParameters();
	SetFibersInitialPositions();
	SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
	printf("Initialisation fibre end\n");

    //-------------- init. device
	printf("Initialisation device start\n");
	InitDevice();
    Create_Device_Ptrs();
    Create_Host_Ptrs();
    Create_Kernels();
	printf("Initialisation device end\n");

    //-------------- copy etat initial host-> device
    Copy_Utils_host2dev();      //à faire une fois suffit
	Copy_Fiber_host2dev();		//initialisation
    Copy_Param_host2dev();		//idem

	//-------------- auvegarde config initiale
	SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
	printf("Initializations finished\n");

    //-------------- main loop
    for (;;)
        DoOneIteration();

    return 0;
}