
//-----------------------------------------------------
__attribute__((alyays_inline)) uint FastCheck(uint fiber1, uint i1, uint fiber2, uint i2) {
	if ((fiber1 == fiber2) && (i2 < i1 + 5)) return SKIP_CONTACT_DETECTION;
	return 0;
}

//-----------------------------------------------------------
__kernel void kernel_Specific_OneTime(
	__global uint* param_Uint,
	__global double* param_Double) {

	param_Uint[0]++;

	//--------- lambda
	if (param_Uint[0] < param_Uint[22])
		param_Double[3] = 0.1;
	else
		param_Double[3] = 0.01;

	//--------- lx ly
	if (param_Double[50] < param_Double[53])
		param_Double[50] += .001;

	if (param_Double[51] < param_Double[54])
		param_Double[51] += .001;
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

	uint						iSegment, iFiber, i, k;
	double						lx, ly;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	lx=param_Double[50];
	ly=param_Double[51];

	#define GN_CELL_X			16
	#define GN_CELL_Y			16

	//--------------------- traversantes
	if (iFiber+1< GN_CELL_Y){
		if (i==0){
			xt[3 * iSegment + 0] = -(GN_CELL_X/2)*lx;
			xt[3 * iSegment + 1] = (-GN_CELL_Y/2+(double)(iFiber)+0.214)*ly;
			xt[3 * iSegment + 2] = 0.;
		}
		else if (i+1==n[iFiber]){
			xt[3 * iSegment + 0] = +(GN_CELL_X/2)*lx;
			xt[3 * iSegment + 1] = (-GN_CELL_Y/2+(double)(iFiber)+0.214)*ly;
			xt[3 * iSegment + 2] = 0.;
		}
	}
	//---------------------- haut gauche
	else if (iFiber +1 == GN_CELL_Y){
		if (i==0){
			xt[3 * iSegment + 0] = -(GN_CELL_X/2)*lx;
			xt[3 * iSegment + 1] = (-GN_CELL_Y/2+(double)(iFiber)+0.214)*ly;
			xt[3 * iSegment + 2] = 0.;
		}
		else if (i+1==n[iFiber]){
			xt[3 * iSegment + 0] = ((double)(iFiber+1)-GN_CELL_Y - GN_CELL_X / 2. + 1. / 4.) * lx;
			xt[3 * iSegment + 1] = (GN_CELL_Y / 2 ) * ly;
			xt[3 * iSegment + 2] = 0.;
		}
	}
	//---------------------- bas
	else if ((iFiber +1 > GN_CELL_Y) && (iFiber + 1 <= GN_CELL_Y + GN_CELL_X)) {
		if (i == 0) {
			xt[3 * iSegment + 0] = ((double)(iFiber+1)-GN_CELL_Y - GN_CELL_X / 2. - 3. / 4.) * lx;
			xt[3 * iSegment + 1] = -(GN_CELL_Y / 2) * ly;
			xt[3 * iSegment + 2] = 0.;
		}
		else if (i + 1 == n[iFiber]) {
			xt[3 * iSegment + 0] = ((double)(iFiber+1)-GN_CELL_Y - GN_CELL_X / 2. - 1. / 4.) * lx;
			xt[3 * iSegment + 1] = -(GN_CELL_Y / 2) * ly;
			xt[3 * iSegment + 2] = 0.;
		}
	}
	//--------------------- haut centre
	else if (iFiber+1 < GN_CELL_Y + 2 * GN_CELL_X) {
		if (i == 0) {
			xt[3 * iSegment + 0] = ((double)(iFiber+1)-GN_CELL_Y - GN_CELL_X - GN_CELL_X / 2. - 1. / 4.) * lx;
			xt[3 * iSegment + 1] = +(GN_CELL_Y / 2) * ly;
			xt[3 * iSegment + 2] = 0.;
		}
		else if (i + 1 == n[iFiber]) {
			xt[3 * iSegment + 0] = ((double)(iFiber+1)-GN_CELL_Y - GN_CELL_X - GN_CELL_X / 2. + 1. / 4.) * lx;
			xt[3 * iSegment + 1] = +(GN_CELL_Y / 2) * ly;
			xt[3 * iSegment + 2] = 0.;
		}
	}
	//--------------------- haut droit
	else if (iFiber +1 == GN_CELL_Y + 2 * GN_CELL_X) {
		if (i == 0) {
			xt[3 * iSegment + 0] = ((double)(iFiber+1)-GN_CELL_Y - GN_CELL_X - GN_CELL_X / 2. - 1. / 4.) * lx;
			xt[3 * iSegment + 1] = +(GN_CELL_Y / 2) * ly;
			xt[3 * iSegment + 2] = 0.;
		}
		else if (i + 1 == n[iFiber]) {
			xt[3 * iSegment + 0] = +(GN_CELL_X / 2) * lx;
			xt[3 * iSegment + 1] = (GN_CELL_Y / 2 - 1 + 0.214) * ly;
			xt[3 * iSegment + 2] = 0.;
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
	double				bound;

	iSegment = (uint)(get_global_id(0));
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];

	#define GN_CELL_X			16
	#define GN_CELL_Y			16

	//----------- force en z plaques
	bound= param_Double[50]*GN_CELL_X/2.;
	if(xt[3 * iSegment + 0] > bound)
		f[3 * iSegment + 0] -= 1.*(xt[3 * iSegment + 0]-bound);

	if(xt[3 * iSegment + 0] < -bound)
		f[3 * iSegment + 0] -= 1.*(xt[3 * iSegment + 0]+bound);

	bound= param_Double[51]*GN_CELL_Y/2.;
	if(xt[3 * iSegment + 1] > bound)
		f[3 * iSegment + 1] -= 1.*(xt[3 * iSegment + 1]-bound);

	if(xt[3 * iSegment + 1] < -bound)
		f[3 * iSegment + 1] -= 1.*(xt[3 * iSegment + 1]+bound);

	return;		
}