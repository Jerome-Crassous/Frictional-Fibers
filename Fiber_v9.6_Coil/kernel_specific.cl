//-----------------------------------------------------
__attribute__((always_inline)) uint FastCheck(uint fiber1, uint i1, uint fiber2, uint i2) {
	return SKIP_CONTACT_DETECTION;
}

//-----------------------------------------------------------
__kernel void kernel_Specific_OneTime(__global uint* param_Int, __global float* param_Float) {

	param_Int[0]++;
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

	uint						iSegment, iFiber, i, j, k, isc, dx, flag0;
	float						l, sc, deltac;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];
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

	uint				iSegment, iFiber, i, flag0;
	float				force, delta, r, r2, arc;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	arc = 2. * PI * 
		sqrt(param_Double[27] * param_Double[27] + param_Double[28] * param_Double[28]);

	force =  convert_float(param_UInt[0]) * 1.e-10 / arc;
	if (i < convert_uint(arc) ) f[3 * iSegment + 2] -= force;
	if ((i + convert_uint(arc)) > n[iFiber]) f[3 * iSegment + 2] += force;
}