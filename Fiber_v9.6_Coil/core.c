#include	"main.h"

extern FILE* g_FilePtr;
extern parameter* gParameterPtr;
extern fiber* gFiber;
extern cl_double* param_Float_hostPtr;
extern contact* gContact;

void PrintContact(contact* theContactPtr, parameter* theParameterPtr, unsigned int index);

//-----------------------------------------------------------------
void DoOneIteration(void) {

	unsigned int			iFiber, i, HalfArc;
	static int				inited = 0;
	double					eBending, eStreching, eTwist;

	//------- force reset + elastic + visqueux
	Execute_Kernel_Calculate_le(gParameterPtr);
	Execute_Kernel_Bending_Force(gParameterPtr);
	Execute_Kernel_Twist_Force(gParameterPtr);

	//------- force specific
	Execute_Kernel_Specific_Force(gParameterPtr);

	//------- tout contact
	Execute_Kernel_PossibleContact(gParameterPtr);
	Execute_Kernel_Contact(gParameterPtr);
	Execute_Kernel_RemoveDoubleContact(gParameterPtr);
	Execute_Kernel_CalculateContactForce(gParameterPtr);
	Execute_Kernel_AddContactForce(gParameterPtr);
	Execute_Kernel_ShiftContact(gParameterPtr);

	//-------- integrate
	Execute_kernel_Integrate_and_Shift(gParameterPtr);

	//-------- geometry
	Execute_kernel_Compute_m1bar(gParameterPtr);
	Execute_kernel_Specific_Position(gParameterPtr);
	Execute_kernel_Compute_m1(gParameterPtr);

	//-------- update MaxDisplacement
	Execute_Kernel_Max_Displacement(gParameterPtr);
	Execute_Kernel_Unwarp(gParameterPtr);

	//---------- arret si trop proche de la traversée
	if ((gParameterPtr->iter > 10*KI) && (Execute_Kernel_Max_Delta(gParameterPtr) > 0.3))
		exit(44);

	//------------ taches millenaires
	if (gParameterPtr->iter % (10 * KI) == 0) {

		//----------- transfer
		Copy_Fiber_dev2host();
		SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
		
		eBending = BendingEnergyLib(gFiber, gParameterPtr);
		eStreching = TractionEnergyLib(gFiber, gParameterPtr);
		eTwist = TwistEnergyLib(gFiber, gParameterPtr);	
		HalfArc =  PI *
			sqrt(gParameterPtr->RHelix* gParameterPtr->RHelix 
				+ gParameterPtr->pitch * gParameterPtr->pitch);
		//---------- force
		printf("iter %u - dx %e - %e %e %e\n",
			gParameterPtr->iter,
			gFiber[0].xt[gFiber[0].n - 1 - HalfArc][2] - gFiber[0].xt[HalfArc][2],
			eBending, eStreching, eTwist);

		fprintf(g_FilePtr, "%u %e %e %e %e\n",
			gParameterPtr->iter,
			gFiber[0].xt[gFiber[0].n - 1 - HalfArc][2] - gFiber[0].xt[HalfArc][2],
			eBending, eStreching, eTwist);
	}

	//------------ arret
	if (gParameterPtr->iter == gParameterPtr->iterStop)
		exit(55);

	//------------ post-operation
	fflush(g_FilePtr);
	gParameterPtr->iter++;
}

//-----------------------------------------------------------
void PrintContact(contact* theContactPtr, parameter* theParameterPtr, unsigned int index) {

	unsigned int		iContact, k;
	double				s1, s2;
	FILE* filePtr;
	char				fileName[256];

	sprintf(fileName, "C:\\data\\test\\contact_%d_%d.txt", gParameterPtr->line, index);
	filePtr = fopen(fileName, "w");
	fprintf(filePtr, "fiber1 fiber2 s1 s2 Cx Cy Cz nx ny nz fcx fcy fcz\n");
	if (!filePtr)
		printf("creation %s failed\n", fileName);

	for (iContact = 0; iContact < theParameterPtr->nContact; iContact++) {
		fprintf(filePtr, "%d %d %e %e %e %e %e %e %e %e %e %e %e\n",
			theContactPtr[iContact].fiber1,
			theContactPtr[iContact].fiber2,
			(((double)theContactPtr[iContact].node1) + theContactPtr[iContact].s1)
			* gFiber[theContactPtr[iContact].fiber1].l0,
			(((double)theContactPtr[iContact].node2) + theContactPtr[iContact].s2)
			* gFiber[theContactPtr[iContact].fiber2].l0,
			theContactPtr[iContact].C_t[0], theContactPtr[iContact].C_t[1], theContactPtr[iContact].C_t[2],
			theContactPtr[iContact].n_t[0], theContactPtr[iContact].n_t[1], theContactPtr[iContact].n_t[2],
			theContactPtr[iContact].fc[0], theContactPtr[iContact].fc[1], theContactPtr[iContact].fc[2]);
		fprintf(g_FilePtr, "%e ", theContactPtr[iContact].fc[2]);
	}
	fprintf(g_FilePtr, "\n");
	fclose(filePtr);
	return;
}