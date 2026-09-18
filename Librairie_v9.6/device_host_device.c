#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdlib.h>
#include <stdio.h>
#include <CL/cl.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

void check_error(cl_int err, const char* operation, char* filename, int line);
#define checkError(E, S) check_error(E, S, __FILE__, __LINE__)

#define __CL_ENABLE_EXCEPTIONS

extern unsigned int* gIter;
extern fiber* gFiber;
extern contact* gContact;
extern parameter* gParameterPtr;
extern cl_context					context;
extern cl_command_queue				command_queue;
extern FILE* g_FilePtr;

//----------------------------------- parameters
cl_uint* param_UInt_hostPtr;
cl_double* param_Double_hostPtr;
cl_mem			param_UInt_devPtr;
cl_mem			param_Double_devPtr;

//----------------------------------- fibres
cl_uint* fiber_n_hostPtr, * fiber_status_hostPtr;											//Nf
cl_float* fiber_radius_hostPtr, * fiber_l0_hostPtr, * fiber_k0_hostPtr;						//Nf
cl_float* fiber_mass_hostPtr, * fiber_b_hostPtr, * fiber_c_hostPtr, * fiber_j0_hostPtr;		//Nf

cl_uint* fiber_flag_hostPtr;																//Nseg
cl_float* fiber_lt_hostPtr, * fiber_ltm_hostPtr;											//Nseg
cl_float* fiber_kappa1_bar_hostPtr, * fiber_kappa2_bar_hostPtr;								//Nseg

cl_double* fiber_xt_hostPtr, * fiber_xtm_hostPtr;								//3*Nseg
cl_double* fiber_thetat_hostPtr, * fiber_thetatm_hostPtr;						//3*Nseg
cl_float* fiber_f_hostPtr, * fiber_e_hostPtr, *fiber_moment_hostPtr;			//3*Nseg
cl_float* fiber_m1_bar_hostPtr, * fiber_m1_hostPtr;								//3*Nseg

cl_mem			fiber_n_devPtr, fiber_status_devPtr;
cl_mem			fiber_radius_devPtr, fiber_l0_devPtr, fiber_k0_devPtr;
cl_mem			fiber_mass_devPtr, fiber_b_devPtr, fiber_c_devPtr, fiber_j0_devPtr;

cl_mem			fiber_flag_devPtr;
cl_mem			fiber_lt_devPtr, fiber_ltm_devPtr;
cl_mem			fiber_kappa1_bar_devPtr, fiber_kappa2_bar_devPtr;

cl_mem			fiber_xt_devPtr, fiber_xtm_devPtr;
cl_mem			fiber_thetat_devPtr, fiber_thetatm_devPtr;
cl_mem			fiber_f_devPtr, fiber_e_devPtr, fiber_moment_devPtr;
cl_mem			fiber_m1_bar_devPtr, fiber_m1_devPtr;
cl_mem			fiber_dissip_visc_global_devPtr, fiber_dissip_visc_stretch_devPtr, fiber_Wop_devPtr;

//---------------------------------- util segment to i
cl_uint			* iFiberFromSegment_hostPtr, * iFromSegment_hostPtr;	//Nseg
cl_mem			iFiberFromSegment_devPtr, iFromSegment_devPtr;

//----------------------------------- possContacts
cl_uint			*n_possCont_in_wg_hostPtr;								//[NcontactMax/WG]
cl_uint			*possContSeg1_hostPtr, *possContSeg2_hostPtr;			//[NcontactMax]
cl_uint			*oldPossContSeg1_hostPtr, *oldPossContSeg2_hostPtr;
cl_uint			*iShiftCont_hostPtr;
cl_mem			n_possCont_in_wg_devPtr;
cl_mem			possContSeg1_devPtr, possContSeg2_devPtr;
cl_mem			oldPossContSeg1_devPtr, oldPossContSeg2_devPtr;
cl_mem			iShiftCont_devPtr;

//----------------------------------- contacts			
cl_uint			*cont_exist_hostPtr;				//[NcontactMax]
cl_float		*cont_s1_hostPtr, *cont_s2_hostPtr;
cl_double		*cont_deltat_hostPtr, *cont_deltatm_hostPtr;
cl_double		*cont_W1_hostPtr, *cont_W2_hostPtr;
cl_float		*cont_Ct_hostPtr;					//[3*NcontactMax]
cl_float		*cont_nt_hostPtr;
cl_float		*cont_ut_hostPtr, * cont_utm_hostPtr;
cl_float		*cont_fc_hostPtr;

cl_mem			cont_exist_devPtr, * old_cont_exist_devPtr;
cl_mem			cont_s1_devPtr, cont_s2_devPtr;
cl_mem			cont_deltat_devPtr, cont_deltatm_devPtr, old_cont_deltatm_devPtr;
cl_mem			cont_W1_devPtr, cont_W2_devPtr;
cl_mem			cont_dissip_ft_devPtr, cont_dissip_fn_devPtr;
cl_mem			cont_Ct_devPtr;
cl_mem			cont_nt_devPtr;
cl_mem			cont_ut_devPtr, cont_utm_devPtr, old_cont_utm_devPtr;
cl_mem			cont_fc_devPtr;

//----------------------------------- force per segment	
cl_uint* nForce_per_Segment_hostPtr, * nMoment_per_Segment_hostPtr;		//NSEG
cl_float* force_per_Segment_hostPtr, * moment_per_Segment_hostPtr;		//3*NSEG*NFPS
cl_mem			nForce_per_Segment_devPtr, nMoment_per_Segment_devPtr;
cl_mem			force_per_Segment_devPtr, moment_per_Segment_devPtr;

//------------------------ reductions --- only on device
cl_mem			reduction_devPtr, partialCount_devPtr;

//--------------------------------------------------------------------
void Create_Host_Ptrs(void) {

	cl_uint				nFiber = gParameterPtr->nFiber;
	cl_uint				nSegment = gParameterPtr->nSegment;
	cl_uint				nContactMax = gParameterPtr->nContactMax;

	//---------- parametres
	param_UInt_hostPtr = (cl_uint*)malloc(NPARAM * sizeof(cl_uint));
	param_Double_hostPtr = (cl_double*)malloc(NPARAM * sizeof(cl_double));

	//---------- fibres
	fiber_n_hostPtr = (cl_uint*)malloc(nFiber * sizeof(cl_uint));
	fiber_status_hostPtr = (cl_uint*)malloc(nFiber * sizeof(cl_uint));
	fiber_radius_hostPtr = (cl_float*)malloc(nFiber * sizeof(cl_float));
	fiber_l0_hostPtr = (cl_float*)malloc(nFiber * sizeof(cl_float));
	fiber_k0_hostPtr = (cl_float*)malloc(nFiber * sizeof(cl_float));
	fiber_mass_hostPtr = (cl_float*)malloc(nFiber * sizeof(cl_float));
	fiber_b_hostPtr = (cl_float*)malloc(nFiber * sizeof(cl_float));
	fiber_c_hostPtr = (cl_float*)malloc(nFiber * sizeof(cl_float));
	fiber_j0_hostPtr = (cl_float*)malloc(nFiber * sizeof(cl_float));

	fiber_flag_hostPtr = (cl_uint*)malloc( nSegment * sizeof(cl_uint));
	fiber_lt_hostPtr = (cl_float*)malloc(nSegment * sizeof(cl_float));
	fiber_ltm_hostPtr = (cl_float*)malloc(nSegment * sizeof(cl_float));
	fiber_kappa1_bar_hostPtr = (cl_float*)malloc(nSegment * sizeof(cl_float));
	fiber_kappa2_bar_hostPtr = (cl_float*)malloc(nSegment * sizeof(cl_float));

	fiber_xt_hostPtr = (cl_double*)malloc(3 * nSegment * sizeof(cl_double));
	fiber_xtm_hostPtr = (cl_double*)malloc(3 * nSegment * sizeof(cl_double));
	fiber_thetat_hostPtr = (cl_double*)malloc(3 * nSegment * sizeof(cl_double));
	fiber_thetatm_hostPtr = (cl_double*)malloc(3 * nSegment * sizeof(cl_double));
	fiber_f_hostPtr = (cl_float*)malloc(3 * nSegment * sizeof(cl_float));
	fiber_e_hostPtr = (cl_float*)malloc(3 * nSegment * sizeof(cl_float));
	fiber_moment_hostPtr = (cl_float*)malloc(3 * nSegment * sizeof(cl_float));
	fiber_m1_bar_hostPtr = (cl_float*)malloc(3 * nSegment * sizeof(cl_float));
	fiber_m1_hostPtr = (cl_float*)malloc(3 * nSegment * sizeof(cl_float));

	//---------- utils
	iFiberFromSegment_hostPtr = (cl_uint*)malloc(nSegment * sizeof(cl_uint));
	iFromSegment_hostPtr = (cl_uint*)malloc(nSegment * sizeof(cl_uint));

	//----------- possContacts
	n_possCont_in_wg_hostPtr = (cl_uint*)malloc(nContactMax / WG * sizeof(cl_uint));
	possContSeg1_hostPtr = (cl_uint*)malloc(nContactMax * sizeof(cl_uint));
	possContSeg2_hostPtr = (cl_uint*)malloc(nContactMax * sizeof(cl_uint));
	oldPossContSeg1_hostPtr = (cl_uint*)malloc(nContactMax * sizeof(cl_uint));
	oldPossContSeg2_hostPtr = (cl_uint*)malloc(nContactMax * sizeof(cl_uint));
	iShiftCont_hostPtr = (cl_uint*)malloc(nContactMax * sizeof(cl_uint));

	//---------- contact		
	cont_exist_hostPtr = (cl_uint*)malloc(nContactMax * sizeof(cl_uint));			//[NcontactMax]
	cont_s1_hostPtr = (cl_float*)malloc(nContactMax * sizeof(cl_float));
	cont_s2_hostPtr = (cl_float*)malloc(nContactMax * sizeof(cl_float));
	cont_deltat_hostPtr = (cl_double*)malloc(nContactMax * sizeof(cl_double));
	cont_deltatm_hostPtr = (cl_double*)malloc(nContactMax * sizeof(cl_double));
	cont_W1_hostPtr = (cl_double*)malloc(nContactMax * sizeof(cl_double));
	cont_W2_hostPtr = (cl_double*)malloc(nContactMax * sizeof(cl_double));
	cont_Ct_hostPtr = (cl_float*)malloc(3 * nContactMax * sizeof(cl_float));		//[3*NcontactMax]
	cont_nt_hostPtr = (cl_float*)malloc(3 * nContactMax * sizeof(cl_float));
	cont_ut_hostPtr = (cl_float*)malloc(3 * nContactMax * sizeof(cl_float));
	cont_utm_hostPtr = (cl_float*)malloc(3 * nContactMax * sizeof(cl_float));
	cont_fc_hostPtr = (cl_float*)malloc(3 * nContactMax * sizeof(cl_float));


	//------------ force per segment
	nForce_per_Segment_hostPtr = (cl_uint*)malloc(nSegment * sizeof(cl_uint));
	nMoment_per_Segment_hostPtr = (cl_uint*)malloc(nSegment * sizeof(cl_uint));
	force_per_Segment_hostPtr = (cl_float*)malloc(3 * nSegment * NFPS * sizeof(cl_float));
	moment_per_Segment_hostPtr = (cl_float*)malloc(3 * nSegment * NFPS * sizeof(cl_float));
}

//--------------------------------------------------------------------
void Create_Device_Ptrs(void) {

	cl_uint				nFiber = gParameterPtr->nFiber;
	cl_uint				nSegment = gParameterPtr->nSegment;
	cl_uint				nContactMax = gParameterPtr->nContactMax;
	cl_int				ret;

	//---------- parametres
	param_UInt_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, NPARAM * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	param_Double_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, NPARAM * sizeof(cl_double), NULL, &ret);checkError(ret, "clCreateBuffer");

	//---------- fibres
	fiber_n_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_status_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_radius_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_l0_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_k0_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_mass_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_b_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_c_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_j0_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nFiber * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");

	fiber_flag_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_uint), NULL, &ret); checkError(ret, "clCreateBuffer");
	fiber_lt_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_ltm_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_kappa1_bar_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_kappa2_bar_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");

	fiber_xt_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_double), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_xtm_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_double), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_thetat_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_double), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_thetatm_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_double), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_f_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_e_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_moment_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_m1_bar_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_m1_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	fiber_dissip_visc_global_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	fiber_dissip_visc_stretch_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	fiber_Wop_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");

	//------------- utils
	iFiberFromSegment_devPtr = clCreateBuffer(context, CL_MEM_READ_ONLY, nSegment * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	iFromSegment_devPtr = clCreateBuffer(context, CL_MEM_READ_ONLY, nSegment * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");

	//---------------- possContacts
	n_possCont_in_wg_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, (nContactMax / WG) * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	possContSeg1_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	oldPossContSeg1_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	possContSeg2_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	oldPossContSeg2_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	iShiftCont_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");

	//---------------- contacts
	cont_exist_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_uint), NULL, &ret); checkError(ret, "clCreateBuffer");
	old_cont_exist_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_uint), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_s1_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_s2_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_deltat_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_deltatm_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	old_cont_deltatm_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_Ct_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_nt_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_ut_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_utm_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	old_cont_utm_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_fc_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nContactMax * sizeof(cl_float), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_dissip_ft_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_dissip_fn_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_W1_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");
	cont_W2_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nContactMax * sizeof(cl_double), NULL, &ret); checkError(ret, "clCreateBuffer");

	//--------------- force per segment
	nForce_per_Segment_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	nMoment_per_Segment_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, nSegment * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
	force_per_Segment_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * NFPS * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	moment_per_Segment_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, 3 * nSegment * NFPS * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");

	//---------- reductions
	reduction_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, (64 * WG * WG / SIZE_REDUCTION) * sizeof(cl_float), NULL, &ret);checkError(ret, "clCreateBuffer");
	partialCount_devPtr = clCreateBuffer(context, CL_MEM_READ_WRITE, (nContactMax / SIZE_REDUCTION) * sizeof(cl_uint), NULL, &ret);checkError(ret, "clCreateBuffer");
}

//--------------------------------------------------------------------
void Copy_Utils_host2dev(void) {

	cl_uint				i, iFiber;
	cl_uint				nFiber = gParameterPtr->nFiber;
	cl_uint				nSegment = gParameterPtr->nSegment;
	cl_int				ret;

	//-------------------------------------------------------------
	for (iFiber = 0; iFiber < nFiber; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			iFromSegment_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)] = i;
			iFiberFromSegment_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)] = iFiber;
		}
	}

	//----------------- copy host -> devive
	ret = clEnqueueWriteBuffer(command_queue, iFromSegment_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_uint), iFromSegment_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, iFiberFromSegment_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_uint), iFiberFromSegment_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
}

//--------------------------------------------------------------------
void Copy_Param_host2dev(void) {

	cl_int					ret;

	param_UInt_hostPtr[0] = gParameterPtr->iter;
	param_UInt_hostPtr[1] = gParameterPtr->line;
	param_UInt_hostPtr[2] = gParameterPtr->nFiber;
	param_UInt_hostPtr[3] = gParameterPtr->nSegment;
	param_UInt_hostPtr[4] = gParameterPtr->nContactMax;
	param_UInt_hostPtr[5] = gParameterPtr->nContact; //semble ne pas être utilisé. A remplacer par param_UInt_hostPtr[7]
	param_UInt_hostPtr[6] = gParameterPtr->periodic;
	//	param_UInt_hostPtr[7] Nombre de contacts;

	param_UInt_hostPtr[20] = gParameterPtr->iter0;
	param_UInt_hostPtr[21] = gParameterPtr->iter1;
	param_UInt_hostPtr[22] = gParameterPtr->low_High_Eta;
	param_UInt_hostPtr[23] = gParameterPtr->iter3;
	param_UInt_hostPtr[24] = gParameterPtr->rRate;
	param_UInt_hostPtr[25] = gParameterPtr->low_High_Mu;
	param_UInt_hostPtr[26] = gParameterPtr->iter6;
	param_UInt_hostPtr[27] = gParameterPtr->iter7;
	param_UInt_hostPtr[28] = gParameterPtr->iter8;
	param_UInt_hostPtr[29] = gParameterPtr->iterStop;

	param_UInt_hostPtr[30] = gParameterPtr->flag0;
	param_UInt_hostPtr[31] = gParameterPtr->flag1;
	param_UInt_hostPtr[32] = gParameterPtr->flag2;
	param_UInt_hostPtr[33] = gParameterPtr->flag3;
	param_UInt_hostPtr[34] = gParameterPtr->flag4;
	param_UInt_hostPtr[35] = gParameterPtr->flag5;
	param_UInt_hostPtr[36] = gParameterPtr->flag6;
	param_UInt_hostPtr[37] = gParameterPtr->flag7;
	param_UInt_hostPtr[38] = gParameterPtr->flag_Phase;
	param_UInt_hostPtr[39] = gParameterPtr->flag_Reset_Contact;

	param_UInt_hostPtr[40] = gParameterPtr->NN;

	//param_UInt_hostPtr[NPARAM - 2]	//reservé nombre de parametre pour dissipation
	//param_UInt_hostPtr[NPARAM - 1]	//reservé numero slot pour sortie dissipation

	//--------- double
	param_Double_hostPtr[0] = gParameterPtr->dt;
	param_Double_hostPtr[1] = gParameterPtr->kn;
	param_Double_hostPtr[2] = gParameterPtr->kt;
	param_Double_hostPtr[3] = gParameterPtr->lambda;
	param_Double_hostPtr[4] = gParameterPtr->lambda_internal;
	param_Double_hostPtr[5] = gParameterPtr->lambda_contact_n;
	param_Double_hostPtr[6] = gParameterPtr->lambda_contact_t;
	param_Double_hostPtr[7] = gParameterPtr->mu;
	param_Double_hostPtr[10] = gParameterPtr->muTarget;
	param_Double_hostPtr[11] = gParameterPtr->lambdaTarget;

	param_Double_hostPtr[15] = gParameterPtr->traction;
	param_Double_hostPtr[16] = gParameterPtr->pressure;
	param_Double_hostPtr[17] = gParameterPtr->pressureTarget;
	param_Double_hostPtr[18] = gParameterPtr->d;
	param_Double_hostPtr[19] = gParameterPtr->dTarget;
	param_Double_hostPtr[20] = gParameterPtr->r;
	param_Double_hostPtr[21] = gParameterPtr->rTarget;
	param_Double_hostPtr[22] = gParameterPtr->R;
	param_Double_hostPtr[23] = gParameterPtr->RTarget;
	param_Double_hostPtr[24] = gParameterPtr->h;
	param_Double_hostPtr[25] = gParameterPtr->hTarget;
	param_Double_hostPtr[26] = gParameterPtr->dhdt;
	param_Double_hostPtr[27] = gParameterPtr->RHelix;
	param_Double_hostPtr[28] = gParameterPtr->pitch;

	param_Double_hostPtr[30] = gParameterPtr->epsStar;
	//	param_Double_hostPtr[31]	//maxDisplacement de l'Iteration
	//	param_Double_hostPtr[32]	//integratedDispacement
	param_Double_hostPtr[33] = 0;	//dissipation visqueuse globale
	param_Double_hostPtr[34] = 0;	//dissipation visqueuse stretch
	param_Double_hostPtr[35] = 0;	//dissipation visqueuse fn
	param_Double_hostPtr[36] = 0;	//dissipation friction solide
	param_Double_hostPtr[37] = 0;	//W operateur

	param_Double_hostPtr[40] = gParameterPtr->fix;
	param_Double_hostPtr[41] = gParameterPtr->fx;
	param_Double_hostPtr[42] = gParameterPtr->fy;
	param_Double_hostPtr[43] = gParameterPtr->fz;
	param_Double_hostPtr[47] = gParameterPtr->sc;
	param_Double_hostPtr[48] = gParameterPtr->thetac;

	param_Double_hostPtr[50] = gParameterPtr->lx;
	param_Double_hostPtr[51] = gParameterPtr->ly;
	param_Double_hostPtr[52] = gParameterPtr->lz;
	param_Double_hostPtr[53] = gParameterPtr->lxTarget;
	param_Double_hostPtr[54] = gParameterPtr->lyTarget;
	param_Double_hostPtr[55] = gParameterPtr->lzTarget;
	param_Double_hostPtr[56] = gParameterPtr->nCellX;
	param_Double_hostPtr[57] = gParameterPtr->nCellY;

	param_Double_hostPtr[60] = gParameterPtr->bound;
	param_Double_hostPtr[61] = gParameterPtr->boundX;
	param_Double_hostPtr[62] = gParameterPtr->boundY;
	param_Double_hostPtr[63] = gParameterPtr->boundZ;
	param_Double_hostPtr[64] = gParameterPtr->alpha;
	param_Double_hostPtr[65] = gParameterPtr->alphaX;
	param_Double_hostPtr[66] = gParameterPtr->alphaY;
	param_Double_hostPtr[67] = gParameterPtr->alphaZ;
	param_Double_hostPtr[68] = gParameterPtr->theta;
	param_Double_hostPtr[69] = gParameterPtr->beta;

	param_Double_hostPtr[70] = gParameterPtr->delta;

	//----------------- copy host -> devive
	ret = clEnqueueWriteBuffer(command_queue, param_UInt_devPtr, CL_TRUE, 0, NPARAM * sizeof(cl_uint), param_UInt_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, param_Double_devPtr, CL_TRUE, 0, NPARAM * sizeof(cl_double), param_Double_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
}

//--------------------------------------------------------------------
void Copy_Param_dev2host(void) {

	cl_int			ret;

	//----------------- copy device -> host
	ret = clEnqueueReadBuffer(command_queue, param_UInt_devPtr, CL_TRUE, 0, NPARAM * sizeof(cl_uint), param_UInt_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, param_Double_devPtr, CL_TRUE, 0, NPARAM * sizeof(cl_double), param_Double_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
}

//--------------------------------------------------------------------
void Write_One_UInt_Param_host2dev(unsigned int slot, unsigned int value) {

	cl_int					ret;
	if (slot == 7)
		printf("############ Warning Write_One_UInt_Param_host2dev ############\n");
	ret = clEnqueueWriteBuffer(command_queue, param_UInt_devPtr, CL_TRUE, slot * sizeof(cl_uint), sizeof(cl_uint), &value, 0, NULL, NULL); checkError(ret, "clEnqueueReadBuffer");
}

//--------------------------------------------------------------------
void Write_One_Double_Param_host2dev(unsigned int slot, double value) {

	cl_int					ret;
	if ((slot == 31) || (slot == 32))
		printf("############ Warning Write_One_Double_Param_host2dev ############\n");
	ret = clEnqueueWriteBuffer(command_queue, param_Double_devPtr, CL_TRUE, slot * sizeof(cl_double), sizeof(cl_double), &value, 0, NULL, NULL); checkError(ret, "clEnqueueReadBuffer");
}

//--------------------------------------------------------------------
unsigned int Read_One_UInt_Param_dev2host(unsigned int slot) {

	cl_int					ret;
	unsigned int value;

	ret = clEnqueueReadBuffer(command_queue, param_UInt_devPtr, CL_TRUE, slot * sizeof(cl_uint), sizeof(cl_uint), &value, 0, NULL, NULL); checkError(ret, "clEnqueueReadBuffer");
	return value;
}

//--------------------------------------------------------------------
double Read_One_Double_Param_dev2host(unsigned int slot) {

	cl_int					ret;
	double					value;

	ret = clEnqueueReadBuffer(command_queue, param_Double_devPtr, CL_TRUE, slot * sizeof(cl_double), sizeof(cl_double), &value, 0, NULL, NULL); checkError(ret, "clEnqueueReadBuffer");
	return value;
}

//--------------------------------------------------------------------
void Copy_Fiber_host2dev(void) {

	cl_uint				i, k, iFiber;
	cl_uint				nFiber = gParameterPtr->nFiber;
	cl_uint				nSegment = gParameterPtr->nSegment;
	cl_int				ret;

	//----------------- wrap gFiber -> host
	for (iFiber = 0; iFiber < nFiber; iFiber++) {
		fiber_n_hostPtr[iFiber] = gFiber[iFiber].n;
		fiber_status_hostPtr[iFiber] = gFiber[iFiber].status;
		fiber_radius_hostPtr[iFiber] = gFiber[iFiber].radius;
		fiber_mass_hostPtr[iFiber] = gFiber[iFiber].masse;
		fiber_l0_hostPtr[iFiber] = gFiber[iFiber].l0;
		fiber_k0_hostPtr[iFiber] = gFiber[iFiber].k0;
		fiber_b_hostPtr[iFiber] = gFiber[iFiber].bending;
		fiber_c_hostPtr[iFiber] = gFiber[iFiber].c;
		fiber_j0_hostPtr[iFiber] = gFiber[iFiber].j;
		for (i = 0; i < gFiber[iFiber].n; i++) {
			fiber_flag_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)] = gFiber[iFiber].flag[i];
			fiber_lt_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)] = gFiber[iFiber].lt[i];
			fiber_ltm_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)] = gFiber[iFiber].ltm[i];
			fiber_kappa1_bar_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)] = gFiber[iFiber].kappa1_bar[i];
			fiber_kappa2_bar_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)] = gFiber[iFiber].kappa2_bar[i];
			for (k = 0; k < 3; k++) {
				fiber_xt_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].xt[i][k];
				fiber_xtm_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].xtm[i][k];
				fiber_thetat_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].thetat[i][k];
				fiber_thetatm_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].thetatm[i][k];
				fiber_f_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].f[i][k];
				fiber_e_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].e[i][k];
				fiber_moment_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].moment[i][k];
				fiber_m1_bar_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].m1_bar[i][k];
				fiber_m1_hostPtr[3 * NSegment(gFiber, gParameterPtr, iFiber, i) + k] = gFiber[iFiber].m1[i][k];
			}
		}
	}

	//----------------- copy host -> devive
	ret = clEnqueueWriteBuffer(command_queue, fiber_n_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_uint), fiber_n_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_status_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_uint), fiber_status_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_radius_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_radius_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_mass_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_mass_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_l0_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_l0_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_k0_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_k0_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_b_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_b_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_c_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_c_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_j0_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_j0_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");

	ret = clEnqueueWriteBuffer(command_queue, fiber_flag_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_uint), fiber_flag_hostPtr, 0, NULL, NULL); checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_lt_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_lt_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_ltm_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_ltm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_kappa1_bar_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_kappa1_bar_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_kappa2_bar_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_kappa2_bar_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");

	ret = clEnqueueWriteBuffer(command_queue, fiber_xt_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_xt_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_xtm_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_xtm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_thetat_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_thetat_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_thetatm_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_thetatm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_f_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_f_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_e_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_e_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_moment_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_moment_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_m1_bar_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_m1_bar_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, fiber_m1_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_m1_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
}

//--------------------------------------------------------------------
void Copy_Fiber_dev2host(void) {

	cl_uint				i, j, iFiber;
	cl_uint				nFiber = gParameterPtr->nFiber;
	cl_uint				nSegment = gParameterPtr->nSegment;
	cl_int				ret;

	//----------------- copy device -> host
	ret = clEnqueueReadBuffer(command_queue, fiber_n_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_uint), fiber_n_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_status_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_uint), fiber_status_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_radius_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_radius_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_mass_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_mass_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_l0_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_l0_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_k0_devPtr, CL_TRUE, 0, nFiber * sizeof(cl_float), fiber_k0_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	ret = clEnqueueReadBuffer(command_queue, fiber_flag_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_uint), fiber_flag_hostPtr, 0, NULL, NULL); checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_lt_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_lt_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_ltm_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_ltm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_kappa1_bar_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_kappa1_bar_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_kappa2_bar_devPtr, CL_TRUE, 0, nSegment * sizeof(cl_float), fiber_kappa2_bar_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	ret = clEnqueueReadBuffer(command_queue, fiber_xt_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_xt_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_xtm_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_xtm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_thetat_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_thetat_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_thetatm_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_thetatm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_f_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_f_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_e_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_e_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_m1_bar_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_m1_bar_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, fiber_m1_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_float), fiber_m1_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	//----------------- wrap host -> gFiber 
	for (iFiber = 0; iFiber < nFiber; iFiber++) {
		gFiber[iFiber].n = fiber_n_hostPtr[iFiber];
		gFiber[iFiber].status = fiber_status_hostPtr[iFiber];
		gFiber[iFiber].radius = fiber_radius_hostPtr[iFiber];
		gFiber[iFiber].masse = fiber_mass_hostPtr[iFiber];
		gFiber[iFiber].l0 = fiber_l0_hostPtr[iFiber];
		gFiber[iFiber].k0 = fiber_k0_hostPtr[iFiber];
		for (i = 0; i < gFiber[iFiber].n; i++) {
			gFiber[iFiber].flag[i] = fiber_flag_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)];
			gFiber[iFiber].lt[i] = fiber_lt_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)];
			gFiber[iFiber].ltm[i] = fiber_ltm_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i)];
			for (j = 0;j < 3;j++) {
				gFiber[iFiber].moment[i][j] = fiber_moment_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].xt[i][j] = fiber_xt_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].xtm[i][j] = fiber_xtm_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].thetat[i][j] = fiber_thetat_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].thetatm[i][j] = fiber_thetatm_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].f[i][j] = fiber_f_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].e[i][j] = fiber_e_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].m1_bar[i][j] = fiber_m1_bar_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
				gFiber[iFiber].m1[i][j] = fiber_m1_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
			}
		}
	}
}

//--------------------------------------------------------------------
void Copy_FiberPosition_dev2host(void) {

	cl_uint				i, j, iFiber;
	cl_uint				nFiber = gParameterPtr->nFiber;
	cl_uint				nSegment = gParameterPtr->nSegment;
	cl_int				ret;

	//----------------- copy device -> host
	ret = clEnqueueReadBuffer(command_queue, fiber_xt_devPtr, CL_TRUE, 0, 3 * nSegment * sizeof(cl_double), fiber_xt_hostPtr, 0, NULL, NULL); checkError(ret, "clEnqueueReadBuffer");

	//----------------- wrap host -> gFiber 
	for (iFiber = 0; iFiber < nFiber; iFiber++) {
		for (i = 0; i < gFiber[iFiber].n; i++) {
			for (j = 0; j < 3; j++) {
				gFiber[iFiber].xt[i][j] = fiber_xt_hostPtr[NSegment(gFiber, gParameterPtr, iFiber, i) * 3 + j];
			}
		}
	}
}

//--------------------------------------------------------------------
void Copy_Contact_dev2host(void) {

	unsigned int			j, k, nC;
	cl_uint					nContactMax = gParameterPtr->nContactMax;
	cl_int					ret;

	//--------------------------------------------------------------------	
	ret = clEnqueueReadBuffer(command_queue, cont_exist_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), cont_exist_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_s1_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_float), cont_s1_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_s2_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_float), cont_s2_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_deltat_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_double), cont_deltat_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_deltatm_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_double), cont_deltatm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_W1_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_double), cont_W1_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_W2_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_double), cont_W2_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_Ct_devPtr, CL_TRUE, 0, nContactMax * 3 * sizeof(cl_float), cont_Ct_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_nt_devPtr, CL_TRUE, 0, nContactMax * 3 * sizeof(cl_float), cont_nt_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_ut_devPtr, CL_TRUE, 0, nContactMax * 3 * sizeof(cl_float), cont_ut_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_utm_devPtr, CL_TRUE, 0, nContactMax * 3 * sizeof(cl_float), cont_utm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_fc_devPtr, CL_TRUE, 0, nContactMax * 3 * sizeof(cl_float), cont_fc_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	//--------------------------------------------------------------------	
	ret = clEnqueueReadBuffer(command_queue, possContSeg1_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), possContSeg1_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, possContSeg2_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), possContSeg2_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	//-------------------------------------------------------------------
	nC = 0;
	for (k = 0; k < nContactMax; k++) {
		if (cont_exist_hostPtr[k] == CONTACT) {
			gContact[nC].exist = cont_exist_hostPtr[k];
			gContact[nC].fiber1 = iFiberFromSegment_hostPtr[possContSeg1_hostPtr[k]];
			gContact[nC].fiber2 = iFiberFromSegment_hostPtr[possContSeg2_hostPtr[k]];
			gContact[nC].node1 = iFromSegment_hostPtr[possContSeg1_hostPtr[k]];
			gContact[nC].node2 = iFromSegment_hostPtr[possContSeg2_hostPtr[k]];
			gContact[nC].s1 = cont_s1_hostPtr[k];
			gContact[nC].s2 = cont_s2_hostPtr[k];
			gContact[nC].deltat = cont_deltat_hostPtr[k];
			gContact[nC].deltatm = cont_deltatm_hostPtr[k];
			gContact[nC].W1 = cont_W1_hostPtr[k];
			gContact[nC].W2 = cont_W2_hostPtr[k];
			for (j = 0; j < 3; j++) {
				gContact[nC].C_t[j] = cont_Ct_hostPtr[3 * k + j];
				gContact[nC].u_t[j] = cont_ut_hostPtr[3 * k + j];
				gContact[nC].u_tm[j] = cont_utm_hostPtr[3 * k + j];
				gContact[nC].n_t[j] = cont_nt_hostPtr[3 * k + j];
				gContact[nC].fc[j] = cont_fc_hostPtr[3 * k + j];
			}
			nC++;
		}
	}
	gParameterPtr->nContact = nC;
}

//--------------------------------------------------------------------
void Print_Possible_dev2host(unsigned int flag) {

	unsigned int			k;
	cl_uint					nContactMax = gParameterPtr->nContactMax;
	cl_int					ret;
	FILE* filePtr;
	char					fileName[256];

	//--------------------------------------------------------------------
	sprintf(fileName, "c://data//possContact_%d_%d.txt", gParameterPtr->iter, flag);
	filePtr = fopen(fileName, "w");

	//--------------------------------------------------------------------	
	ret = clEnqueueReadBuffer(command_queue, possContSeg1_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), possContSeg1_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, possContSeg2_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), possContSeg2_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, oldPossContSeg1_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), oldPossContSeg1_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, oldPossContSeg2_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), oldPossContSeg2_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_exist_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), cont_exist_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");
	ret = clEnqueueReadBuffer(command_queue, cont_utm_devPtr, CL_TRUE, 0, nContactMax * 3 * sizeof(cl_float), cont_utm_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	for (k = 0; k < WG; k++) {
		if (cont_exist_hostPtr[k] != NOCONTACT) {
			fprintf(filePtr, "%d %d - %d %d - %d %e\n",
				oldPossContSeg1_hostPtr[k], oldPossContSeg2_hostPtr[k], possContSeg1_hostPtr[k], possContSeg2_hostPtr[k],
				cont_exist_hostPtr[k],
				cont_utm_hostPtr[3 * k + 2]);
		}
	}
	fclose(filePtr);
}

//--------------------------------------------------------------------
unsigned int  NPossContact_dev(void) {

	unsigned int   		k, nPossC = 0;
	cl_uint				nContactMax = gParameterPtr->nContactMax;
	cl_int				ret;

	//--------------------------------------------------------------------	
	ret = clEnqueueReadBuffer(command_queue, n_possCont_in_wg_devPtr, CL_TRUE, 0, nContactMax / WG * sizeof(cl_uint), n_possCont_in_wg_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	for (k = 0; k < nContactMax / WG; k++)
		nPossC += n_possCont_in_wg_hostPtr[k];
	return nPossC;
}

//--------------------------------------------------------------------
unsigned int   NContact_dev(void) {

	unsigned int  		k, nC = 0;
	cl_uint				nContactMax = gParameterPtr->nContactMax;
	cl_int				ret;

	//--------------------------------------------------------------------	
	ret = clEnqueueReadBuffer(command_queue, cont_exist_devPtr, CL_TRUE, 0, nContactMax * sizeof(cl_uint), cont_exist_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueReadBuffer");

	for (k = 0; k < nContactMax; k++) {
		if (cont_exist_hostPtr[k] != NOCONTACT) {
			nC++;
		}
	}
	return nC;
}

//--------------------------------------------------------------------
void ResetDissipationBuffers(void) {

	cl_double			zero = 0.;
	cl_int				ret;

	ret = clEnqueueFillBuffer(command_queue, fiber_dissip_visc_global_devPtr, &zero, sizeof(cl_double), 0, gParameterPtr->nSegment * sizeof(cl_double), 0, NULL, NULL); checkError(ret, "clEnqueueFillBuffer");
	ret = clEnqueueFillBuffer(command_queue, fiber_dissip_visc_stretch_devPtr, &zero, sizeof(cl_double), 0, gParameterPtr->nSegment * sizeof(cl_double), 0, NULL, NULL); checkError(ret, "clEnqueueFillBuffer");
	ret = clEnqueueFillBuffer(command_queue, cont_dissip_ft_devPtr, &zero, sizeof(cl_double), 0, gParameterPtr->nContactMax * sizeof(cl_double), 0, NULL, NULL); checkError(ret, "clEnqueueFillBuffer");
	ret = clEnqueueFillBuffer(command_queue, cont_dissip_fn_devPtr, &zero, sizeof(cl_double), 0, gParameterPtr->nContactMax * sizeof(cl_double), 0, NULL, NULL); checkError(ret, "clEnqueueFillBuffer");
	ret = clEnqueueFillBuffer(command_queue, fiber_Wop_devPtr, &zero, sizeof(cl_double), 0, gParameterPtr->nSegment * sizeof(cl_double), 0, NULL, NULL); checkError(ret, "clEnqueueFillBuffer");
}

//--------------------------------------------------------------------
void ResetW12Buffers(void) {

	cl_double			zero = 0.;
	cl_int				ret;

	ret = clEnqueueFillBuffer(command_queue, cont_W1_devPtr, &zero, sizeof(cl_double), 0, gParameterPtr->nContactMax * sizeof(cl_double), 0, NULL, NULL); checkError(ret, "clEnqueueFillBuffer");
	ret = clEnqueueFillBuffer(command_queue, cont_W2_devPtr, &zero, sizeof(cl_double), 0, gParameterPtr->nContactMax * sizeof(cl_double), 0, NULL, NULL); checkError(ret, "clEnqueueFillBuffer");
}