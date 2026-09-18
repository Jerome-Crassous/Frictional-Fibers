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
	gParameterPtr->flag_Reset_Contact = 1;

	gParameterPtr->low_High_Eta = 10 * KI;		//low eta - high eta
	gParameterPtr->iter1 = 150 * MI;	//start harmonic compression 
	gParameterPtr->iter3 = 150 * KI;	//start compression
	gParameterPtr->low_High_Mu = 25 * KI;		//start mu
	gParameterPtr->iter6 = 5 * MI;		//start cylinder expansion

	gParameterPtr->periodic = PERIODIC_NO;
	gParameterPtr->periodicVisu = PERIODIC_VISU_NO;

	gParameterPtr->flag_Reset_Contact = 1;
	gParameterPtr->iter = 0;
	gParameterPtr->periodic = 0;
	gParameterPtr->nSegment = 0;
	gParameterPtr->nContactMax = (512*WG);
	gParameterPtr->epsStar = 0.1;

	gParameterPtr->bound = 100.;
	gParameterPtr->r = 0.5;
	gParameterPtr->rTarget = 0.5;

	gParameterPtr->line = 10;
	gParameterPtr->nFiber = 31;
	gParameterPtr->hTarget = 1.;

	gParameterPtr->dt = 0.1;
	gParameterPtr->kn = 1.;
	gParameterPtr->kt = 0.5;
	gParameterPtr->lambda = 0.001;
	gParameterPtr->lambda_internal = 2.8;
	gParameterPtr->lambda_contact_n = 1.;
	gParameterPtr->lambda_contact_t = 0.;
	gParameterPtr->bending = 0.001;
	gParameterPtr->c = gParameterPtr->bending;
	gParameterPtr->mu = 0.5;
	gParameterPtr->muTarget = gParameterPtr->mu;

	//------------- overload des arguments
	printf("Nombre d’arguments passes au programme : %d\n", argc);
	if (argc != 1) {
		for (i = 1; i < argc; i++)
			printf("argv[%d] : %s\n", i, argv[i]);
		ok = sscanf(argv[1], "%u", &(gParameterPtr->line));
		ok = sscanf(argv[2], "%u", &(gParameterPtr->flag0));
	}

	//---------- 15 contacts/fibre * 50
//	gParameterPtr->nContactMax = 15 * gParameterPtr->nFiber * 50;
//	gParameterPtr->nContactMax = ((gParameterPtr->nContactMax-1) / (WG*WG) +1)*WG*WG;
//	printf("gParameterPtr->nContactMax = %u\n", gParameterPtr->nContactMax);

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

    //-------------- copy etat initial
    Copy_Utils_host2dev();      //à faire une fois suffit
	Copy_Fiber_host2dev();		//initialisation
    Copy_Param_host2dev();		//idem

	SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
	printf("Initializations finished\n");
//	exit(33);

    //-------------- main loop
    for (;;)
        DoOneIteration();

    return 0;
}