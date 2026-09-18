#include	"main.h"

extern FILE* g_FilePtr;
extern parameter* gParameterPtr;
extern fiber* gFiber;
extern contact* gContact;

//-----------------------------------------------------------------
void DoOneIteration(void) {

	double				eBending, eKinetic;
	
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

	//-------- taches millénaires
	if (gParameterPtr->iter % (KI) == 0) {
		Copy_Fiber_dev2host();
		SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
		SaveReducedConfigLib(gFiber, gParameterPtr, "c:\\data\\multiple_coils\\", gParameterPtr->iter/KI);
		eBending = BendingEnergyLib(gFiber, gParameterPtr);
		eKinetic = KineticEnergyLib(gFiber, gParameterPtr);
		Execute_Kernel_Dissipation(gParameterPtr);
		printf("%u %e %e %e %e \n", gParameterPtr->iter, MyTimeLib(), eBending, eKinetic, eBending + eKinetic);
		fprintf(g_FilePtr, "%u %e %e %e %e \n", gParameterPtr->iter, MyTimeLib(), eBending, eKinetic, eBending + eKinetic);
	}

	if (gParameterPtr->iter == gParameterPtr->iterStop) {
		static FILE* summaryFilePtr;
		summaryFilePtr = fopen("summary.txt", "a");
		fprintf(summaryFilePtr, "%u %e \n", gParameterPtr->nContactMax, MyTimeLib());
		#ifdef MODE_DROP
			SaveConfigLib(gFiber, gParameterPtr,"", -1);
		#endif // MODE_DROP
		exit(87523);
	}

	//------------ post-operation
	fflush(g_FilePtr);	
	gParameterPtr->iter++;
}