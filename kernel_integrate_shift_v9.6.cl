//-----------------------------------------------------------
__kernel void kernel_Integrate_and_Shift(
	__global double* param_Double,
	__global uint* status, __global float* mass, __global float* j0,
	__global uint* iFiberFromSegment,
	__global double* xt, __global double* xtm,__global float* f,
	__global double* thetat, __global double* thetatm, __global float* moment) {

	uint						iSegment, iFiber, k;
	float						aux;

	//-----------------------------------------
	iSegment = get_global_id(0);
	iFiber = iFiberFromSegment[iSegment];
	if (status[iFiber] != STATUS_FREE) return;

	//--------- integrate & shift  of  xt & thetat
	for (k = 0; k < 3; k++) {
		aux = xt[3 * iSegment + k] - xtm[3 * iSegment + k];
		xtm[3* iSegment +k] = xt[3* iSegment +k];
		xt[3* iSegment +k] += aux + f[3* iSegment +k] * param_Double[0] * param_Double[0] / mass[iFiber];	

		aux = thetat[3 * iSegment + k] - thetatm[3 * iSegment + k];
		thetatm[3*iSegment+k] = thetat[3*iSegment+k];
		thetat[3*iSegment+k] += aux + moment[3*iSegment+k] * param_Double[0] * param_Double[0] / j0[iFiber];
	}
}

//-----------------------------------------------------------
__kernel void kernel_Compute_m1bar(
	__global double* param_Double,
	__global uint* n, __global uint* status,
	__global float* l0, __global float* mass, __global float* j0,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global float* lt,
	__global float* e, __global float* f,
	__global double* xt, __global double* xtm,
	__global double* thetat, __global double* thetatm,
	__global float* m1_bar, __global float* m1, __global float* moment) {

	uint						iSegment, iFiber, i, k;
	float3						et,etm, m1_bar_Loc;
	float						aux;

	//-----------------------------------------
	iSegment = get_global_id(0);
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];
	if (status[iFiber] != STATUS_FREE)		return;
	if (n[iFiber] == 1)	return;		//---- cas billes m1_bar sans importance

	//---------- loading	
	if (i +1 < n[iFiber]){	//---- cas segment
		for (k = 0; k < 3; k++){
			et[k] = xt[3*(iSegment +1)+k] - xt[3* iSegment +k];
			etm[k] = xtm[3*(iSegment +1)+k] - xtm[3* iSegment +k];
			m1_bar_Loc[k] = m1_bar[3* iSegment +k];
		}
	}
	else if (i +1 == n[iFiber]){ //---- cas bille terminale chargement segment-1
		for (k = 0; k < 3; k++){
			et[k] = xt[3* iSegment +k] - xt[3* (iSegment-1) +k];
			etm[k] = xtm[3* iSegment +k] - xtm[3* (iSegment-1) +k];
			m1_bar_Loc[k] = m1_bar[3* (iSegment-1) +k];
		}
	}
	else //----------- cas erreur
		printf("what in kernel_Compute_m1bar\n");

	et = normalize(et);
	etm = normalize(etm);
	m1_bar_Loc = parallel_transport(m1_bar_Loc, etm, et);		

	aux = dot(m1_bar_Loc, et);
	for (k = 0; k < 3; k++)
		m1_bar_Loc[k] -= aux * et[k];
	m1_bar_Loc = normalize(m1_bar_Loc);

	//------ sortie m1_bar
	for (k = 0; k < 3; k++)
		m1_bar[3 * iSegment + k] = m1_bar_Loc[k];
}

//-----------------------------------------------------------
__kernel void kernel_Compute_m1(
	__global double* param_Double,
	__global uint* n, __global uint* status,
	__global float* l0, __global float* mass, __global float* j0,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global float* lt,
	__global float* e, __global float* f,
	__global double* xt, __global double* xtm,
	__global double* thetat, __global double* thetatm,
	__global float* m1_bar, __global float* m1, __global float* moment) {

	uint						iSegment, iFiber, i, k;
	float3						et,etm, m1_Loc, vect, dTheta;
	float						aux;

	//-----------------------------------------
	iSegment = get_global_id(0);
	iFiber = iFiberFromSegment[iSegment];
	i = iFromSegment[iSegment];
	if (status[iFiber] != STATUS_FREE)		return;

	//-------- cas billes seules
	if (n[iFiber] == 1){
		dTheta = convert_float3(vload3(iSegment, thetat) - vload3(iSegment, thetatm));
		
		//--- rotate m1
		vect = convert_float3(vload3(iSegment, m1));
		vect += cross(dTheta, vect);
		vect = normalize(vect);
		vstore3(vect, iSegment, m1);

		//--- rotate m1
		vect = convert_float3(vload3(iSegment, e));
		vect += cross(dTheta, vect);
		vect = normalize(vect);
		vstore3(vect, iSegment, e);

		return;
	}

	//--------- cas non billes seules
	if (i +1 < n[iFiber]){	//---- cas segment ou bille normale
		for (k = 0; k < 3; k++){
			et[k] = xt[3*(iSegment +1)+k] - xt[3* iSegment +k];
			m1_Loc[k] = m1_bar[3* iSegment +k];
		}
	}
	else if (i +1 == n[iFiber]){ //---- cas bille terminale chargement segment-1
		for (k = 0; k < 3; k++){
			et[k] = xt[3*(iSegment)+k] - xt[3* (iSegment-1) +k];
			m1_Loc[k] = m1_bar[3* (iSegment-1) +k];
		}
	}
	else //----------- cas erreur
		printf("what in kernel_Compute_m1\n");

	et = normalize(et);
	m1_Loc = rotate(m1_Loc, et, thetat[3 * iSegment + 2]);
	for (k = 0; k < 3; k++)
		m1[3 * iSegment + k] = m1_Loc[k];
}