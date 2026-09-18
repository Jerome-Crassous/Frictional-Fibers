#include	"main.h"

extern FILE* g_FilePtr;
extern parameter* gParameterPtr;
extern fiber* gFiber;
extern contact* gContact;

extern cl_command_queue command_queue;

extern cl_mem		param_Uint_devPtr, param_Double_devPtr;
extern cl_mem		cont_exist_devPtr, cont_deltat_devPtr;
extern cl_mem		possContSeg1_devPtr, possContSeg2_devPtr;
extern cl_mem		cont_s1_devPtr, cont_s2_devPtr;
extern cl_mem		n_possCont_in_wg_devPtr;

extern cl_uint		*param_Uint_hostPtr;
extern cl_double	*param_Double_hostPtr;
extern cl_uint		*cont_exist_hostPtr;
extern cl_float		*cont_deltat_hostPtr;
extern cl_uint		*possContSeg1_hostPtr, *possContSeg2_hostPtr;		
extern cl_float		*cont_s1_hostPtr, *cont_s2_hostPtr;
extern cl_uint		*n_possCont_in_wg_hostPtr;						//[NcontactMax/WG]

//-----------------------------------------------------------------
void DoOneIteration(void) {

	unsigned int			iFiber, i, n;
	double					eBending, eStreching, eTwist, eKinetic, eContactNormal, eContactTangential, dissipFrict;
	FILE* filePtr;

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
		Copy_Contact_dev2host();
		SaveReducedConfigLib(gFiber, gParameterPtr, NULL, -1);

		//---------- affichage
		eTwist = TwistEnergyLib(gFiber, gParameterPtr);
		printf("iter = %u theta = %e enrgy = %e\n",gParameterPtr->iter,gFiber[0].thetat[2][2], eTwist);
		fprintf(g_FilePtr, "%u %e %e\n", gParameterPtr->iter, gFiber[0].thetat[8][2], eTwist);
	}

	//--------- fin
	if (gParameterPtr->iter == gParameterPtr->iterStop) {
		filePtr = fopen("profile.txt", "w");
		for (i = 0; i < gFiber[0].n; i++)
			fprintf(filePtr, "%e %e\n", gFiber[0].xt[i][2], gFiber[0].thetat[i][2]);
		fclose(filePtr);
		exit(6);
	}



	//------------ post-operation
	fflush(g_FilePtr);	
	gParameterPtr->iter++;
}

//-----------------------------------------------------------
void PrintContactDelta(void) {

	unsigned int  		k, nC = 0;
	cl_uint				nContactMax = gParameterPtr->nContactMax;
	cl_int				ret;
	FILE				*filePtr;
	char				fileName[256];
	
	//--------------------------------------------------------------------	
	ret = clEnqueueReadBuffer(command_queue, cont_exist_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), cont_exist_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_deltat_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_float), cont_deltat_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");
	
	//--------------------------------------------------------------------	
	sprintf(fileName, "C:\\data\\packing\\contacts_%d_%d.txt", gParameterPtr->line, gParameterPtr->iter/KI);
	filePtr = fopen(fileName, "w");
	if (!filePtr) {
		printf("opening %s failed\n", fileName);
		return;
	}

	//--------------------------------------------------------------------	
	for (k = 0; k < nContactMax; k++) {
		if (cont_exist_hostPtr[k] != NOCONTACT) {
			fprintf(filePtr, "%e\n", cont_deltat_hostPtr[k]);
		}
	}
	fclose(filePtr);
	return;
}

//-----------------------------------------------------------
void DetailedContact(void) {

	unsigned int  			iFiber, k, nC = 0, iSeg1, iSeg2;
	cl_uint					nContactMax = gParameterPtr->nContactMax;
	cl_int					ret;
	float					sum = 0, delta;
	double					s1, s2;
	static FILE				* filePtr;
	static unsigned int		inited = 0;
	char					fileName[256];

	//-------------------------------------------------------------------
	sprintf(fileName, "C:\\data\\detailed_contact_%d_%d.txt", gParameterPtr->line, gParameterPtr->iter / KI);
	filePtr = fopen(fileName, "w");
	fprintf(filePtr, "fiber1 fiber2 s1 s2 delta\n");

	//-------------------------- ------------------------------------------	
	ret = clEnqueueReadBuffer(command_queue, cont_exist_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), cont_exist_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_deltat_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_float), cont_deltat_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, possContSeg1_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), possContSeg1_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, possContSeg2_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), possContSeg2_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_s1_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_float), cont_s1_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_s2_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_float), cont_s2_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");

	//--------------------------------------------------------------------	
	for (k = 0; k < nContactMax; k++) {
		if (cont_exist_hostPtr[k] != NOCONTACT) {
			iSeg1 = possContSeg1_hostPtr[k];
			iSeg2 = possContSeg2_hostPtr[k];
			s1 = ((double)(possContSeg1_hostPtr[k]%9) + cont_s1_hostPtr[k]);
			s2 = ((double)(possContSeg2_hostPtr[k]%9) + cont_s2_hostPtr[k]);
			delta = cont_deltat_hostPtr[k];
			fprintf(filePtr, "%u %u %e %e %e\n", iSeg1/9, iSeg2/9, s1, s2, delta);
		}
	}
	fclose(filePtr);
	fprintf(filePtr, "\n");
	fflush(filePtr);
	
	return;
}

//---------------------------------------------------------------------
double PrintBendingEnergy(fiber* theFiberPtr, parameter* theParameterPtr) {

	unsigned int		i, iFiber;
	double				aux = 0, omega, energy = 0;
	static FILE* filePtr;
	char					fileName[256];


	//-------------------------------------------------------------------
	sprintf(fileName, "C:\\data\\bending_energy_%d_%d.txt", gParameterPtr->line, gParameterPtr->iter / KI);
	filePtr = fopen(fileName, "w");
	fprintf(filePtr, "iFiber_eb eb\n");

	//-------------------------------------------------------------------
	for (iFiber = 0; iFiber < theParameterPtr->nFiber; iFiber++) {
		if (theFiberPtr[iFiber].status == STATUS_FREE) {
			aux = 0;
			for (i = 1; i < theFiberPtr[iFiber].n - 1; i++) {
				omega = KappaLib(theFiberPtr[iFiber], i, 1) - theFiberPtr[iFiber].kappa1_bar[i];
				aux += omega * omega;
				omega = KappaLib(theFiberPtr[iFiber], i, 2) - theFiberPtr[iFiber].kappa2_bar[i];
				aux += omega * omega;
			}
			aux *= theFiberPtr[iFiber].bending;
			fprintf(filePtr, "%d %e\n", iFiber, aux/2);
		}
		energy += aux;
	}
	fclose(filePtr);
	return energy / 2;
}

//--------------------------------------------------------------------
void Check_PossibleContact(void) {

	int			        nContactMax = gParameterPtr->nContactMax;
	unsigned int		iGroup, iContact, maxPossibleContacts = 0, iGroupMaxPossibleContacts;
	size_t              global_Size[1];
	size_t              local_Size[1];
	cl_int              ret;
	cl_float			integratedDispacement, epsStar;
	char				fileName[256];

	//-----------------------------------------
	static unsigned int	inited = 0;
	static FILE* filePtr;
	if (!inited) {
		sprintf(fileName, "Check_PossibleContact_%u.txt", gParameterPtr->line);
		filePtr = fopen(fileName, "w");
		inited = 1;
	}

	//---------------- lit le nombre de potentiel contact par group et determine le nombre max de potentiel contact	
	ret = clEnqueueReadBuffer(command_queue, n_possCont_in_wg_devPtr, CL_TRUE, 0, nContactMax / WG * sizeof(cl_uint), n_possCont_in_wg_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");

	maxPossibleContacts = 0;
	iGroupMaxPossibleContacts = 0;
	for (iGroup = 0; iGroup < nContactMax / WG; iGroup++) {
		if (n_possCont_in_wg_hostPtr[iGroup] > maxPossibleContacts) {
			maxPossibleContacts = n_possCont_in_wg_hostPtr[iGroup];
			iGroupMaxPossibleContacts = iGroup;
		}
	}

	fprintf(filePtr,"%u %u %u\n",gParameterPtr->iter, maxPossibleContacts, iGroupMaxPossibleContacts);
	fflush(filePtr);
}

