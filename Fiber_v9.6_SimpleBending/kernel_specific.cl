//-----------------------------------------------------
__attribute__((always_inline)) uint FastCheck(uint fiber1, uint i1, uint fiber2, uint i2) {
	return SKIP_CONTACT_DETECTION;
}

//-----------------------------------------------------------
__kernel void kernel_Specific_OneTime(
	__global uint* param_Uint,
	__global double* param_Double) {
	
	param_Uint[0]++;
}	

//-----------------------------------------------------------
__kernel void kernel_Specific_Position(
	__global uint* param_UInt, __global double* param_Double,
	__global uint* n,
	__global float* radius, __global float* l0,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global uint* flag,
	__global float* e,
	__global double* xtm, __global double* xt,
	__global double* thetatm, __global double* thetat,
	__global float* m1_bar, __global float* m1) {


	uint						iSegment, iFiber, i, j, k;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	if ((iFiber == 0) && (i<=1)) {
		for (k = 0; k < 3; k++) {
			xt[3 * iSegment + k] = xtm[3 * iSegment + k];
			thetat[3 * iSegment + k] = thetatm[3 * iSegment + k];
		}
	}
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

	uint				iSegment,iFiber,i;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	if ((iFiber == 0) && (i + 1 == n[iFiber]))
		f[3 * iSegment + 2] -= 1.e-4;
}