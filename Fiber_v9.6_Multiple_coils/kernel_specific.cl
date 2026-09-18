
//-----------------------------------------------------
__attribute__((alyays_inline)) uint FastCheck(uint fiber1, uint i1, uint fiber2, uint i2) {
	if ((fiber1 == fiber2) && (i2 < i1 + 5)) return SKIP_CONTACT_DETECTION;
	return 0;
}

//-----------------------------------------------------------
__kernel void kernel_Specific_OneTime(
	__global uint* param_Uint,
	__global double* param_Double) {

	double dr;	
	param_Uint[0]++;

	//--------- gr
	if (param_Double[20] < param_Double[21]) {
		dr = param_Double[21] / convert_double(param_Uint[24]);
		param_Double[20] += dr;		//r+=dr
		param_Double[32] += dr;		//integratedDisplacement+=dr
	}
}	

//-----------------------------------------------------------
__kernel void kernel_Specific_Position(
	__global uint * param_UInt, __global double* param_Double,
	__global uint * n,
	__global float* radius, __global float* l0,
	__global uint * iFiberFromSegment, __global uint * iFromSegment,
	__global uint* flag, 
	__global float* e,
	__global double* xtm, __global double* xt,
	__global double* thetatm, __global double* thetat,
	__global float* m1_bar, __global float* m1) {
	
	uint						iSegment, iFiber, i, k;
	double						lx, ly;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	if(i==0)
		radius[iFiber]=param_Double[20];

	if ((param_UInt[38]!=0) && (iFiber == 1206) && (i == 0))
		xt[3 * iSegment + 2] = xtm[3 * iSegment + 2] + 1.e-4;
}

//-----------------------------------------------------------
__kernel void kernel_Specific_Force(
	__global uint* param_UInt, __global double* param_Double,
	__global uint* n,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global uint* flag, 
	__global double* xtm, __global double* xt,
	__global double* thetatm, __global double* thetat,
	__global float* f, __global float* moment, __global double* Wop) {

	uint				iSegment,iFiber,i,k;
	double3				bound;
	float				gravity;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	//----------- force en z plaques
	bound[0] = param_Double[50]/2.;
	bound[1] = param_Double[51]/2.;
	bound[2] = param_Double[52]/2.;

	for(k=0;k<3;k++){
		if((xt[3 * iSegment + k] > bound[k] && (k!=2)))
			f[3 * iSegment + k] -= 1.*(xt[3 * iSegment + k]-bound[k]);
		if(xt[3 * iSegment + k] < -bound[k])
			f[3 * iSegment + k] -= 1.*(xt[3 * iSegment + k]+bound[k]);
	}

	//---------- gravité
	if (param_Double[0]< 500*KI)
		gravity = 1.e-6;
	else
		gravity = 1.e-7;

	f[3 * iSegment + 2] -= gravity;

	return;
}