#include "main.h"

parameter* gParameterPtr;
FILE* g_FilePtr, *g_File2Ptr;
fiber* gFiber;
contact* gContact;

//----------------------------------------------------------------------
int main(int argc, char** argv) {

	//------------ init globale parameter
	gParameterPtr = (parameter*)malloc(sizeof(parameter));	

	gParameterPtr->flag_Reset_Contact = 1;
	gParameterPtr->iter = 0;
	gParameterPtr->low_High_Eta = 10 * KI;	//low eta - high eta
	gParameterPtr->rRate = 10 * KI;			//rate radius expansion	

	gParameterPtr->periodic = PERIODIC_NO;
	gParameterPtr->periodicVisu = PERIODIC_VISU_NO;
	gParameterPtr->line = 502;
	gParameterPtr->nContactMax =  (WG * WG);
	gParameterPtr->epsStar = 1.;

	gParameterPtr->nFiber = 2000;

	gParameterPtr->lx = 150.;
	gParameterPtr->ly = 150.;
	gParameterPtr->lz = 200;

	gParameterPtr->r = .1;
	gParameterPtr->rTarget = .5;
	gParameterPtr->RHelix = 4.;
	gParameterPtr->pitch = 3.;

	gParameterPtr->bending = 1.; 
	gParameterPtr->c = gParameterPtr->bending;
	gParameterPtr->mu = 0.5;
	gParameterPtr->muTarget = gParameterPtr->mu;

	gParameterPtr->dt = 0.1;
	gParameterPtr->kn = 1.;
	gParameterPtr->kt = .5;
	gParameterPtr->lambda = 1.e-4; 
	gParameterPtr->lambda_internal = 2.8;
	gParameterPtr->lambda_contact_n = 1.;
	gParameterPtr->lambda_contact_t = 0.;

	//------------- int globales file
	g_FilePtr = fopen("logFile.tmp", "w");

	//-------------- init. des fibres OU lecture configs 
	#ifdef MODE_DROP
		printf("Initialisation fibre start\n");
		gFiber = (fiber*)malloc((gParameterPtr->nFiber) * sizeof(fiber));
		SetFibersParameters();
		SetFibersInitialPositions();
		printf("Initialisation fibre end\n");
		gParameterPtr->flag_Phase = 0;	
		gParameterPtr->iterStop = 1000 * KI;
	#endif

	#ifdef MODE_EXTRACT
		ReadConfigLib(&gFiber, gParameterPtr, "", -1);
		gParameterPtr->flag_Phase = 1;
		gParameterPtr->iterStop = 4000 * KI;
		gParameterPtr->r = gParameterPtr->rTarget;
		int iFiber;
		for (iFiber = 0; iFiber < gParameterPtr->nFiber - 1; iFiber++) {
			if (gFiber[iFiber].xt[0][0] * gFiber[iFiber].xt[0][0] +
				gFiber[iFiber].xt[0][1] * gFiber[iFiber].xt[0][1] < 100.) {
				fprintf(g_FilePtr, "%d %e %e %e\n", iFiber,
					gFiber[iFiber].xt[0][0], gFiber[iFiber].xt[0][1], gFiber[iFiber].xt[0][2]);
			}
			gFiber[iFiber].bending *= 10.;
		}
		//exit(2);
	#endif

	//--------- allocation contact et forcage reset contact
	printf("gParameterPtr->nContactMax %u\n", gParameterPtr->nContactMax);	
	gContact = (contact*)malloc((gParameterPtr->nContactMax) * sizeof(contact));

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

	//-------------- initialisation des temps
	MyTimeLib();
//	exit(44);

    //-------------- main loop
    for (;;)
        DoOneIteration();

    return 0;
}