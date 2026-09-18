//-----------------------------------------------------------
__kernel void kernel_Calculate_le(
	__global double* param_Double, __global uint* n, __global uint* status,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global float* k0, __global float* l0,
	__global float* lt, __global float* ltm,
	__global float* e, __global float* f, __global double* xt, __global double* xtm,
	__global double* thetat, __global double* thetatm, __global float* moment,
	__global double* dissip_visc_global, __global double* dissip_visc_stretch){

	uint					iSegment, i, k, iFiber;
	float					dt, lambda, lambda_internal;
	float3					aux, force_Loc=0, moment_Loc=0.;
	double					llt,lltm,daux;

	//-----------------------
	iSegment = get_global_id(0); 
	iFiber = iFiberFromSegment[iSegment];
	if (status[iFiber] == STATUS_VIRTUAL) return;
	i = iFromSegment[iSegment];

	//------------ lt -- ltm -- e
	if (i + 1 < n[iFiber]) {
		aux = convert_float3(vload3(iSegment + 1, xtm) - vload3(iSegment, xtm));
		ltm[iSegment] = length(aux);

		aux = convert_float3(vload3(iSegment + 1, xt) - vload3(iSegment, xt));
		lt[iSegment] = length(aux);

		aux = normalize(aux);
		vstore3(aux, iSegment, e);
	}

	//------------- traction + visqueux
	dt = param_Double[0];
	lambda = param_Double[3];
	lambda_internal = param_Double[4];

	//------------ force et dissipation visqueuse
	force_Loc = -lambda * convert_float3(vload3(iSegment, xt) - vload3(iSegment, xtm)) / dt;
	dissip_visc_global[iSegment] += param_Double[3] *
		(dot(vload3(iSegment, xt) - vload3(iSegment, xtm), vload3(iSegment, xt) - vload3(iSegment, xtm))) / param_Double[0];

	//----------- force et disspation stretch
	if (i > 0) {
		lltm = length(vload3(iSegment, xtm) - vload3(iSegment - 1, xtm));
		llt = length(vload3(iSegment, xt) - vload3(iSegment - 1, xt));
		force_Loc -= vload3(iSegment - 1, e) * 
			(k0[iFiber] * (((float)llt) / l0[iFiber] - 1.0f) + lambda_internal * ((float)(llt - lltm)) / dt);
	}
	if (i + 1 < (n[iFiber])) {
		lltm = length(vload3(iSegment + 1, xtm) - vload3(iSegment, xtm));
		llt = length(vload3(iSegment + 1, xt) - vload3(iSegment, xt));
		force_Loc += vload3(iSegment, e) * (k0[iFiber] * (((float)llt) / l0[iFiber] - 1.0f) + lambda_internal * ((float)(llt - lltm)) / dt);
		dissip_visc_stretch[iSegment] += lambda_internal * (llt - lltm) * (llt - lltm) / dt;
	}
	vstore3(force_Loc, iSegment, f);

	//------------ moment
	moment_Loc = - lambda * convert_float3(vload3(iSegment, thetat) - vload3(iSegment, thetatm)) / dt;
	dissip_visc_global[iSegment] += param_Double[3] * (thetat[3* iSegment+2]- thetatm[3 * iSegment + 2]) 
		* (thetat[3 * iSegment + 2] - thetatm[3 * iSegment + 2]) / param_Double[0];

	vstore3(moment_Loc, iSegment, moment);
}

//-----------------------------------------------------------
__kernel void kernel_BendingForce(
	__global uint* param_UInt, __global double* param_Double, __global uint* n, __global uint* status,
	__global float* bending,	__global float* l0,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global float* kappa1_bar, __global float* kappa2_bar,
	__global double* xt, __global double* xtm, 
	__global float* e, __global float* f, __global float* moment, __global float* m1) {

	uint		iFiber, iSegment, i, k;
	int			shiftNode, shiftSegment;
	float		kappa1_red, kappa2_red;
	float3		forceLoc =0.;
	float3		momentLoc = 0., momentLocAux=0.;
	float3		dTds = 0., m1Loc = 0., m2Loc = 0.;
	float3		e_im2 = 0., e_im1 = 0., e_i = 0., e_ip1 = 0.;
	float3		m_im2 = 0., m_im1 = 0., m_i = 0., m_ip1 = 0.;
	double3		xt_im2 = 0., xt_im1 = 0., xt_i = 0., xt_ip1 = 0., xt_ip2 = 0.;
	double3		xtm_im2 = 0., xtm_im1 = 0., xtm_i = 0., xtm_ip1 = 0., xtm_ip2 = 0.;

	//-----------------------------------------
	iSegment = get_global_id(0);
	iFiber = iFiberFromSegment[iSegment];
 	if (n[iFiber]==1) return;	
	if (status[iFiber] != STATUS_FREE) return;
	i = iFromSegment[ iSegment];

	//------------------------------------------
	if ((i >= 2) && (i + 1 <= n[iFiber])) {
		e_im2 = vload3(iSegment - 2, e);
		m_im2 = vload3(iSegment - 2, m1);
		xt_im2 = vload3(iSegment - 2, xt);
		xtm_im2 = vload3(iSegment - 2, xtm);
	}
	if (i >= 1) {
		e_im1 = vload3(iSegment - 1, e);
		m_im1 = vload3(iSegment - 1, m1);
		xt_im1 = vload3(iSegment - 1, xt);
		xtm_im1 = vload3(iSegment - 1, xtm);
	}
	if (i + 1 <= n[iFiber]) {
		e_i = vload3(iSegment, e);
		m_i = vload3(iSegment, m1);
		xt_i = vload3(iSegment, xt);
		xtm_i = vload3(iSegment, xtm);
	}	
	if (i + 2 <= n[iFiber]) {
		e_ip1 = vload3(iSegment + 1, e);
		m_ip1 = vload3(iSegment + 1, m1);
		xt_ip1 = vload3(iSegment + 1, xt);
		xtm_ip1 = vload3(iSegment + 1, xtm);
	}
	if (i + 3 <= n[iFiber]) {
		xt_ip2 = vload3(iSegment + 2, xt);
		xtm_ip2 = vload3(iSegment + 2, xtm);
	}

	//------------------------------------------
	if ((i >= 2) && (i + 1 <= n[iFiber])) {
		dTds = convert_float3(xt_im2 - 2 * xt_im1 + xt_i) / l0[iFiber] / l0[iFiber];
		m1Loc = 0.5f * (m_im2 + m_im1);
		m2Loc = 0.5f * (cross(e_im2, m_im2) + cross(e_im1, m_im1));
		kappa1_red = dot(dTds, m1Loc) - kappa1_bar[iSegment - 1];
		kappa2_red = dot(dTds, m2Loc) - kappa2_bar[iSegment - 1];

		forceLoc -= (kappa1_red * m1Loc + kappa2_red * m2Loc) * bending[iFiber] / l0[iFiber];
		momentLocAux = (kappa1_red * cross(dTds, m_im1) + kappa2_red * cross(dTds, cross(e_im1, m_im1))) * 0.5f * bending[iFiber]* l0[iFiber];
		forceLoc -= cross(e_im1, momentLocAux) / l0[iFiber];
	}

	//------------------------------------------
	if ((i >= 1) && (i + 2 <= n[iFiber])) {
		dTds = convert_float3(xt_im1 - 2 * xt_i + xt_ip1) / l0[iFiber] / l0[iFiber];
		m1Loc = 0.5f * (m_im1 + m_i);
		m2Loc = 0.5f * (cross(e_im1, m_im1) + cross(e_i, m_i));
		kappa1_red = dot(dTds, m1Loc) - kappa1_bar[iSegment - 0];
		kappa2_red = dot(dTds, m2Loc) - kappa2_bar[iSegment - 0];

		forceLoc += 2 * (kappa1_red * m1Loc + kappa2_red * m2Loc) * bending[iFiber] / l0[iFiber];

		momentLocAux = (kappa1_red * cross(dTds, m_i) + kappa2_red * cross(dTds, cross(e_i, m_i))) * 0.5f * bending[iFiber]* l0[iFiber];
		momentLoc += momentLocAux;
		forceLoc += cross(e_i, momentLocAux) / l0[iFiber];

		momentLocAux = (kappa1_red * cross(dTds, m_im1) + kappa2_red * cross(dTds, cross(e_im1, m_im1))) * 0.5f * bending[iFiber]* l0[iFiber];
		forceLoc -= cross(e_im1, momentLocAux) / l0[iFiber];
}

	//------------------------------------------
	if ((i >= 0) && (i + 3 <= n[iFiber])) {
		dTds = convert_float3(xt_i - 2 * xt_ip1 + xt_ip2) / l0[iFiber] / l0[iFiber];
		m1Loc = 0.5f * (m_i + m_ip1);
		m2Loc = 0.5f * (cross(e_i, m_i) + cross(e_ip1, m_ip1));
		kappa1_red = dot(dTds, m1Loc) - kappa1_bar[iSegment + 1];
		kappa2_red = dot(dTds, m2Loc) - kappa2_bar[iSegment + 1];

		forceLoc -= (kappa1_red * m1Loc + kappa2_red * m2Loc) * bending[iFiber] / l0[iFiber];

		momentLocAux = (kappa1_red * cross(dTds, m_i) + kappa2_red * cross(dTds, cross(e_i, m_i))) * 0.5f * bending[iFiber] * l0[iFiber];
		momentLoc += momentLocAux;
		forceLoc += cross(e_i, momentLocAux) / l0[iFiber];
	}

	//------------------------------------------
	forceLoc += vload3(iSegment, f);
	vstore3(forceLoc, iSegment, f);
	moment[3 * iSegment + 2] += dot(momentLoc, e_i);
}

//-----------------------------------------------------------
__kernel void kernel_TwistForce(
	__global uint* param_UInt, __global uint* n, __global uint* status,
	__global float* c, __global float* l0,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global float* e, __global float* f, __global float* m1, __global float* moment) {

	uint		iFiber, iSegment, j, i;
	float		beta, a;
	float3		vect = 0., force=0.;
	float3		e_im2 = 0., e_im1 = 0., e_i = 0., e_ip1 = 0.;
	float3		m_im2 = 0., m_im1 = 0., m_i = 0., m_ip1 = 0.;
	float3		moment_sym_im1 = 0., moment_sym_i = 0., moment_sym_ip1 = 0.;

	//-----------------------------------------	
	iSegment = get_global_id(0);
	iFiber = iFiberFromSegment[iSegment];
	if (n[iFiber] == 1) return;
	if (status[iFiber] != STATUS_FREE) return;
	i = iFromSegment[iSegment];

	//------------------------------------------
	if (i >= 2) {
		e_im2 = vload3(iSegment - 2, e);
		m_im2 = vload3(iSegment - 2, m1);
	}
	if (i >= 1) {
		e_im1 = vload3(iSegment - 1, e);
		m_im1 = vload3(iSegment - 1, m1);
	}
	if (i + 2 <= n[iFiber]) {
		e_i = vload3(iSegment, e);
		m_i = vload3(iSegment, m1);
	}	
	if (i + 3 <= n[iFiber]) {
		e_ip1 = vload3(iSegment + 1, e);
		m_ip1 = vload3(iSegment + 1, m1);
	}

	//------------------------------------------
	if ((i >= 2) && (i + 1 <= n[iFiber])) {
		vect = e_im2 + e_im1;
		beta = dot(vect, cross(m_im1,m_im2));
		a = -c[iFiber] / 4. / l0[iFiber] * beta * dot(vect, cross(cross(e_im1, m_im1),m_im2));
		moment_sym_im1 = 0.25f * a * (1.f + dot(e_im2,e_im1)) * (e_im2+e_im1);
	}

	//------------------------------------------
	if ((i >= 1) && (i + 2 <= n[iFiber])) {
		vect = e_im1 + e_i;
		beta = dot(vect, cross(m_i, m_im1));
		a = -c[iFiber] / 4. / l0[iFiber] * beta * dot(vect, cross(cross(e_i, m_i), m_im1));
		moment_sym_i = 0.25f * a * (1.f + dot(e_im1, e_i)) * (e_im1 + e_i);
	}

	//------------------------------------------
	if ((i >= 0) && (i + 3 <= n[iFiber])) {
		vect = e_i + e_ip1;
		beta = dot(vect, cross(m_ip1, m_i));
		a = -c[iFiber] / 4. / l0[iFiber] * beta * dot(vect, cross(cross(e_ip1, m_ip1), m_i));
		moment_sym_ip1 = 0.25f * a * (1.f + dot(e_i, e_ip1)) * (e_i + e_ip1);
	}

	moment[3 * iSegment + 2] += dot(moment_sym_i, e_i) - dot(moment_sym_ip1, e_i);

	//------------------------------------------
	if ((i >= 2) && (i + 1 <= n[iFiber]))
		force -= cross(e_im1, moment_sym_im1);
	if ((i >= 1) && (i + 2 <= n[iFiber])) {
		force += cross(e_im1, moment_sym_i);
		force += cross(e_i, moment_sym_i);
	}
	if ((i >= 0) && (i + 3 <= n[iFiber]))
		force -= cross(e_i, moment_sym_ip1);

	force += vload3(iSegment, f);
	vstore3(force, iSegment, f);
}