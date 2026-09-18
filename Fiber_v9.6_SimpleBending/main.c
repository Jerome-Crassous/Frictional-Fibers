#include "main.h"

parameter* gParameterPtr;
FILE* g_FilePtr, *g_File2Ptr;
fiber* gFiber;
contact* gContact;

//----------------------------------------------------------------------
int main(int argc, char** argv) {

	int         i, j, ok;

	//------------ init globale parameter
	gParameterPtr = (parameter*)malloc(sizeof(parameter));	
	
	gParameterPtr->iterStop = 300 * KI;			//end

	gParameterPtr->iter = 0;
	gParameterPtr->flag_Reset_Contact = 1;
	gParameterPtr->periodic = PERIODIC_NO;
	gParameterPtr->periodicVisu = PERIODIC_VISU_NO;

	gParameterPtr->line = 1000;
	gParameterPtr->nFiber = 2;
	gParameterPtr->nContactMax = (32 * WG);
	gParameterPtr->epsStar = 10.;

	gParameterPtr->ls = 1.5;
	gParameterPtr->r = 1.;

	gParameterPtr->dt = 0.1;
	gParameterPtr->kn = 1.;
	gParameterPtr->kt = .5;
	gParameterPtr->lambda = 1.e-3;
	gParameterPtr->lambda_internal = 2.8;
	gParameterPtr->lambda_contact_n = 1.;
	gParameterPtr->lambda_contact_t = 0.;
	gParameterPtr->bending = 0.1;
	gParameterPtr->c = gParameterPtr->bending;
	gParameterPtr->mu = 0.5;
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

    //-------------- copy etat initial
    Copy_Utils_host2dev();      //à faire une fois suffit
	Copy_Fiber_host2dev();		//initialisation
    Copy_Param_host2dev();		//idem

	SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
	printf("Initializations finished\n");
//	exit(44);

    //-------------- main loop
    for (;;)
        DoOneIteration();

    return 0;
}