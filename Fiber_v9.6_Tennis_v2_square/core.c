#include	"main.h"

extern FILE* g_FilePtr;
extern parameter* gParameterPtr;
extern fiber* gFiber;
extern cl_double* param_Float_hostPtr;
extern contact* gContact;

void PrintContact(contact* theContactPtr, parameter* theParameterPtr, unsigned int index);
void PrintFiber0(void);

//-----------------------------------------------------------------
void DoOneIteration(void) {

	unsigned int			iFiber, i, k, isc, iTest, iCenter;
	double					timeStart, theta, lx, ly, l, force[3], dx, fInt0[3], fInt1[3];
	double					fxx, fyy, fint, t[3], n[3], sc;					
	char					fileName[256];
	static FILE* filePtr;
	static int inited = 0, reset = 0, direction;
	static double			scMin, scMax, dgsc, dgscMin;
	double					eBending, eStreching, eTwist, eKinetic;

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

	//------------ arret
	if (gParameterPtr->iter % (KI) == 0) {
		Copy_Fiber_dev2host();
		Copy_Contact_dev2host();
		SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);
	//	getchar();

		iCenter = (gParameterPtr->line +1) * 10;
		iTest = iCenter + 10;
		printf("%d %e %e %e %\n",
			gParameterPtr->iter,gFiber[0].xt[iCenter][2], gFiber[0].xt[iTest][2],
			DAbs(gFiber[0].xt[iTest][2] - gFiber[0].xtm[iTest][2]));

		fprintf(g_FilePtr,"%d %e %e %e %\n",
			gParameterPtr->iter, gFiber[0].xt[iCenter][2], gFiber[0].xt[iTest][2],
			DAbs(gFiber[0].xt[iTest][2] - gFiber[0].xtm[iTest][2]));


		if (gParameterPtr->iter % (10*KI) == 0) {
			PrintFiber0();
			PrintContact(gContact, gParameterPtr, 1);
		}
	}

	if (gParameterPtr->iter == (50 * 100 * KI))
		exit(41);

	//------------ post-operation
	fflush(g_FilePtr);
	gParameterPtr->iter++;
}

//-----------------------------------------------------------
void PrintFiber0(void) {

	unsigned int		i, k, mu, delta;
	double				s1, s2;
	FILE* filePtr;
	char				fileName[256];

	delta =  (gFiber[0].radius - gFiber[0].xt[210][2]);
	mu = 100 * gParameterPtr->mu;
	printf("### %d %e %d\n", gParameterPtr->iter,
		gFiber[0].radius - gFiber[0].xt[210][2], delta);

	sprintf(fileName, "profil_%d_%d.txt", mu, delta);
	filePtr = fopen(fileName, "w");
	fprintf(filePtr, "xx_%d_%d zz_%d_%d zzRed_%d_%d\n", mu, delta, mu, delta, mu, delta);

	for(i=0;i<gFiber[0].n;i++)
		fprintf(filePtr, "%e %e %e\n",
			gFiber[0].xt[i][0],
			(gFiber[0].radius - gFiber[0].xt[i][2]),
			(gFiber[0].radius-gFiber[0].xt[i][2])/ (gFiber[0].radius - gFiber[0].xt[210][2]));
	fclose(filePtr);
	return;
}

//-----------------------------------------------------------
void PrintContact(contact* theContactPtr, parameter* theParameterPtr, unsigned int index) {

	unsigned int		iContact, k, delta, mu;
	double				s1, s2, fn, ft[3];
	FILE* filePtr;
	char				fileName[256];

	delta = (gFiber[0].radius - gFiber[0].xt[210][2]);
	mu = 100 * gParameterPtr->mu;
	printf("### %d %e %d\n", gParameterPtr->iter,
		gFiber[0].radius - gFiber[0].xt[210][2], delta);

	sprintf(fileName, "vontact_%d_%d.txt", mu, delta);
	filePtr = fopen(fileName, "w");
	fprintf(filePtr, "fiber1 fiber2 s1 s2 Cx Cy Cz nx ny nz fcx fcy fcz mu\n");
	if (!filePtr)
		printf("creation %s failed\n", fileName);

	for (iContact = 0; iContact < theParameterPtr->nContact; iContact++) {

		fn = DotProduct(theContactPtr[iContact].fc, theContactPtr[iContact].n_t);
		for (k = 0; k < 3; k++)
			ft[k] = theContactPtr[iContact].fc[k] - fn * theContactPtr[iContact].n_t[k];

		fprintf(filePtr, "%d %d %e %e %e %e %e %e %e %e %e %e %e %e\n",
			theContactPtr[iContact].fiber1,
			theContactPtr[iContact].fiber2,
			(((double)theContactPtr[iContact].node1) + theContactPtr[iContact].s1)
			* gFiber[theContactPtr[iContact].fiber1].l0,
			(((double)theContactPtr[iContact].node2) + theContactPtr[iContact].s2)
			* gFiber[theContactPtr[iContact].fiber2].l0,
			theContactPtr[iContact].C_t[0], theContactPtr[iContact].C_t[1], theContactPtr[iContact].C_t[2],
			theContactPtr[iContact].n_t[0], theContactPtr[iContact].n_t[1], theContactPtr[iContact].n_t[2],
			theContactPtr[iContact].fc[0], theContactPtr[iContact].fc[1], theContactPtr[iContact].fc[2],
			Norm(ft)/fn
		);
	}
	fclose(filePtr);
	return;
}


//-----------------------------------------------------------
unsigned int MeanContactPosition(contact* theContactPtr, parameter* theParameterPtr,
	unsigned int fiber1, unsigned int fiber2, double* position) {

	//return the mean position of the contact fiber1-fiber2
	unsigned int		iContact, k;
	double				deltaTot = 0;

	for (k = 0;k < 3;k++)
		position[k] = 0;

	for (iContact = 0; iContact < theParameterPtr->nContact; iContact++) {
		if ((theContactPtr[iContact].fiber1 == fiber1)
			&& (theContactPtr[iContact].fiber2 == fiber2)) {

			deltaTot += theContactPtr[iContact].deltat;
			for (k = 0;k < 3;k++)
				position[k] += theContactPtr[iContact].deltat * theContactPtr[iContact].C_t[k];
		}
	}

	if (deltaTot < 1.e-12)
		return 0;
	else {
		for (k = 0;k < 3;k++)
			position[k] /= deltaTot;
		return 1;
	}
}