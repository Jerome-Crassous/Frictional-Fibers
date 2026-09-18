//-----------------------------------------------------
//------- specific kernel must be inserted here
//-----------------------------------------------------
__attribute__((always_inline)) uint FastCheck(uint fiber1, uint i1, uint fiber2, uint i2) {

	if ((fiber1 == fiber2) && (i2 < i1 + 5)) return SKIP_CONTACT_DETECTION;
	return 0;
}

//-----------------------------------------------------------
__kernel void kernel_Specific_OneTime(
	__global uint* param_UInt,
	__global double* param_Double) {

	unsigned int			k, iSegment;
		
	//-----------------------------------------
	param_UInt[0]++;

	//--------- lambda
	if(param_UInt[0]<1000)
		param_Double[3] = 1.e-1;
	else
		param_Double[3] = 1.e-1;
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
	double delta;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	if ((i ==0)||(i==n[iFiber]-1)){
		xt[3 * iSegment + 0] = xtm[3 * iSegment + 0];
		xt[3 * iSegment + 1] = xtm[3 * iSegment + 1];
		xt[3 * iSegment + 2] = xtm[3 * iSegment + 2];
	}

	delta = 1.e-4 * convert_double(param_UInt[0]);
	if(delta>10.) 
		delta = 10.;
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
	
	uint						iFiber, iSegment, i, k;
	float						tension=0.01;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];	

	 if (i == 0){
		if (iFiber < 15) 
			f[3 * iSegment + 0] -= tension;
		else if (iFiber < 30)
			f[3 * iSegment + 1] -= tension;
	}
	else if (i == 300){
		if (iFiber < 15) 
			f[3 * iSegment + 0] += tension;
		else if (iFiber < 30)
			f[3 * iSegment + 1] += tension;
	}
}