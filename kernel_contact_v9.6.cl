
//-----------------------------------------------------------
__kernel void kernel_Contact(
	__global uint* param_UInt, __global double* param_Double,												//------ parameter
	__global uint* fiber_n, __global float* radius, __global double* xt, __global double* xtm,		//------ fibre
	__global uint* iFiberFromSegment, __global uint* iFromSegment,				//------ utils
	__global uint* possContSeg1, __global uint* possContSeg2,					//------ possCont	
	__global uint* iShiftCont,
	__global uint* cont_exist,													//------ cont
	__global float* cont_s1, __global float* cont_s2,
	__global double* cont_deltat, __global double* cont_deltatm,
	__global float* cont_Ct,
	__global float* cont_nt,
	__global float* cont_ut, __global float* cont_utm) {

	uint						iSeg1, iSeg2, iFiber1, iFiber2, i1, i2, n1, n2;
	uint						k, found = 0, iCont, type;
	float3						shift = { 0.,0.,0. };
	double3						a, b, c, n, nAux, CC, dC;
	double3						am, bm, cm;
	double						aa, bb, cc, ab, bc, ac, aux = 0, s1, s2, dist, delta, r12;
	double						aam, bbm, ccm, abm, bcm, acm, auxm = 0, s1m, s2m, distm, deltam;
	double						deltaCC = -1000., deltaCS = -1000., deltaSC = -1000., deltaSS = -1000.;
	double						deltaCCm = -1000., deltaCSm = -1000., deltaSCm = -1000., deltaSSm = -1000.;

	//---------------------------------------------
	iCont = (uint)(get_global_id(0));
	if (cont_exist[iCont] == NOPOSSCONTACT)
							 
												
		return;

	//---------- cherche si il y a un possible contact
	iSeg1 = possContSeg1[iCont];
	iSeg2 = possContSeg2[iCont];
	iFiber1 = iFiberFromSegment[iSeg1];
	iFiber2 = iFiberFromSegment[iSeg2];
	i1 = iFromSegment[iSeg1];
	i2 = iFromSegment[iSeg2];
	n1 = fiber_n[iFiber1];
	n2 = fiber_n[iFiber2];

	ComputeShift(param_UInt[6],iShiftCont[iCont], param_Double[50], param_Double[51], param_Double[52], &shift);

	//--------------- distances
	for (k = 0; k < 3; k++) {
		a[k] = xt[3 * iSeg2 + k] + shift[k] - xt[3 * iSeg1 + k];
		am[k] = xtm[3 * iSeg2 + k] + shift[k] - xtm[3 * iSeg1 + k];
		if (i1 + 2 <= n1) {
			b[k] = (xt[3 * (iSeg1 + 1) + k] - xt[3 * iSeg1 + k]);
			bm[k] = (xtm[3 * (iSeg1 + 1) + k] - xtm[3 * iSeg1 + k]);
		}
		else {
			b[k] = 0;
			bm[k] = 0;
		}
		if (i2 + 2 <= n2) {
			c[k] = (xt[3 * (iSeg2 + 1) + k] - xt[3 * iSeg2 + k]);
			cm[k] = (xtm[3 * (iSeg2 + 1) + k] - xtm[3 * iSeg2 + k]);
		}
		else {
			c[k] = 0;
			cm[k] = 0;
		}
	}
	aa = dot(a, a), bb = dot(b, b), cc = dot(c, c);
	ab = dot(a, b), ac = dot(a, c), bc = dot(b, c);
	aam = dot(am, am), bbm = dot(bm, bm), ccm = dot(cm, cm);
	abm = dot(am, bm), acm = dot(am, cm), bcm = dot(bm, cm);
	r12 = radius[iFiber1] + radius[iFiber2];

	//------------------------ SS
	dist = sqrt(aa);
	deltaSS = r12 - dist;
	distm = sqrt(aam);
	deltaSSm = r12 - distm;

	//------------------------ CS
	if (i1 + 2 <= n1) {
		s1 = ab / bb;
		s1m = abm / bbm;
		if ((s1 > 0) && (s1 < 1)) {
			dist = sqrt(fmax(0, (aa - 2. * ab * s1 + bb * s1 * s1)));
			deltaCS = r12 - dist;
			distm = sqrt(fmax(0, (aam - 2. * abm * s1m + bbm * s1m * s1m)));
			deltaCSm = r12 - distm;
		}
	}

	//------------------------ SC
	if (i2 + 2 <= n2) {
		s2 = -ac / cc;
		s2m = -acm / ccm;
		if ((s2 > 0) && (s2 < 1)) {
			dist = sqrt(fmax(0, (aa + 2 * ac * s2 + cc * s2 * s2)));
			deltaSC = r12 - dist;
			distm = sqrt(fmax(0, (aam + 2 * acm * s2m + ccm * s2m * s2m)));
			deltaSCm = r12 - distm;
		}
	}

	//------------------------ CC
	if ((i1 + 2 <= n1) && (i2 + 2 <= n2)) {
		aux = bb * cc - bc * bc;
		auxm = bbm * ccm - bcm * bcm;
		if (aux < EPS) 	//----------- segments paralleles
			;
		else {
			s1 = (cc * ab - ac * bc) / aux;
			s2 = -(bb * ac - ab * bc) / aux;
			s1m = (ccm * abm - acm * bcm) / auxm;
			s2m = -(bbm * acm - abm * bcm) / auxm;
			if ((s1 > 0) && (s1 < 1) && (s2 > 0) && (s2 < 1)) {
				dist = sqrt(fmax(0, (aa - 2. * ab * s1 + 2 * ac * s2 - 2 * bc * s1 * s2 + bb * s1 * s1 + cc * s2 * s2)));
				deltaCC = r12 - dist;
				distm = sqrt(fmax(0, (aam - 2. * abm * s1m + 2 * acm * s2m - 2 * bcm * s1m * s2m + bbm * s1m * s1m + ccm * s2m * s2m)));
				deltaCCm = r12 - distm;
			}
		}
	}

	//--------------------- cherche le plus grand delta	
	if ((deltaSS > 0) && (deltaSS > deltaCS) && (deltaSS > deltaSC) && (deltaSS > deltaCC))
		s1 = 0, s2 = 0, delta = deltaSS, deltam = deltaSSm, found = 1;
	else if ((deltaCS > 0) && (deltaCS > deltaSC) && (deltaCS > deltaCC))
		s1 = ab / bb, s2 = 0, delta = deltaCS, deltam = deltaCSm, found = 1;
	else if ((deltaSC > 0) && (deltaSC > deltaCC))
		s1 = 0, s2 = -ac / cc, delta = deltaSC, deltam = deltaSCm, found = 1;
	else if (deltaCC > 0) {
		aux = bb * cc - bc * bc;
		if (aux < EPS) 	//----------- segments paralleles
			;
		else {
			s1 = (cc * ab - ac * bc) / aux;
			s2 = -(bb * ac - ab * bc) / aux;
			delta = deltaCC, deltam = deltaCCm, found = 1;
		}
	}
	else {
		cont_exist[iCont] = NOCONTACT;
		return;
	}

	//------------------ Fill contact
	//------------------ normale --- point de contact
	for (k = 0; k < 3; k++) {
		n[k] = a[k] + s2 * c[k] - s1 * b[k];
		nAux[k] = n[k];
	}
	if (length(n) < EPS)
		printf("what kernel_Contact %d %d %e %e %e %e %e %e?\n", iSeg1, iSeg2, length(n), s1, s2, a[0], a[1], a[2]);
	else if (length(n) > EPS) {
		n = normalize(n);
		if (dot(n, nAux) < 0.)
			printf("what kernel_Contact %e\n", dot(n, nAux));
	}
	for (k = 0; k < 3; k++)
		CC[k] = (xt[3 * iSeg1 + k] + s1 * b[k] + radius[iFiber1] * n[k]
			+ xt[3 * iSeg2 + k] + shift[k] + s2 * c[k] - radius[iFiber2] * n[k]) / 2.;

												  
									 
					  
	 
					

	//-------------- remplit le contact	
						  
	cont_s1[iCont] = s1;
	cont_s2[iCont] = s2;
	cont_deltat[iCont] = delta;
							  

	for (k = 0; k < 3; k++) {
		cont_nt[3 * iCont + k] = n[k];
		cont_Ct[3 * iCont + k] = CC[k];
	}

	if (cont_exist[iCont] == NOCONTACT) { //if contact did not exist in previous iteration initialise utm
		cont_deltatm[iCont] = delta;
		for (k = 0; k < 3; k++) {
			cont_utm[3 * iCont + k] = 0.;
		}
	}
	cont_exist[iCont] = CONTACT;
}

//-----------------------------------------------------------
__kernel void kernel_RemoveDoubleContact(
	__global uint* param_UInt,											//------ parameter
	__global uint* possContSeg1, __global uint* possContSeg2,			//------ possCont					
	__global uint* cont_exist,											//------ cont
	__global double* cont_deltat, __global float* cont_Ct) {

	uint						k, iCont, iContP;
	uint						iSeg1, iSeg2, iSeg1P, iSeg2P;
	float						d2;

	//--------------------------------------------------------
/*	iCont = (uint)(get_global_id(0));
	if (possContSeg1[iCont] == UINT_MAX)	return;
	if (possContSeg2[iCont] == UINT_MAX)	return;
	if ((cont_exist[iCont] == NOCONTACT) || (cont_exist[iCont] == OLDCONTACT)) return;

	//------- il faut chercher quels contacts ?
	for (iContP = 0; iContP < param_UInt[4]; iContP++) {

		if ((cont_exist[iContP] == NOCONTACT) || (cont_exist[iContP] == OLDCONTACT)) return;
		if (possContSeg1[iContP] == UINT_MAX)	return;
		if (possContSeg2[iContP] == UINT_MAX)	return;

		iSeg1 = possContSeg1[iCont];
		iSeg2 = possContSeg2[iCont];
		iSeg1P = possContSeg1[iContP];
		iSeg2P = possContSeg2[iContP];

		if ((abs(iSeg1 - iSeg1P) * abs(iSeg2 - iSeg2P) > 1) || (iCont == iContP)) return;

		d2 = 0;
		for (k = 0; k < 3; k++)
			d2 += (cont_Ct[3 * iCont + k] - cont_Ct[3 * iContP + k])
			* (cont_Ct[3 * iCont + k] - cont_Ct[3 * iContP + k]);

		if (d2 < EPS_DEP * EPS_DEP) {	//---c'est un double contact
			if (cont_deltat[iContP] < cont_deltat[iCont])
				cont_exist[iContP] = OLDCONTACT;		//-------- race fail possible ?
			else
				cont_exist[iCont] = OLDCONTACT;
			return;
		}
	}*/
}

//-----------------------------------------------------------
__kernel void kernel_CalculateContactForce(
	__global uint* param_UInt, __global double* param_Double,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global uint* fiber_n,
	__global float* fiber_radius, __global float* fiber_l0,
	__global double* fiber_xt, __global double* fiber_xtm,
	__global double* fiber_thetat, __global double* fiber_thetatm,
	__global float* fiber_e,
	__global uint* possContSeg1, __global uint* possContSeg2, __global uint* iShiftCont,
	__global uint* cont_exist,
	__global float* cont_s1, __global float* cont_s2,
	__global double* cont_deltat, __global double* cont_deltatm,
	__global float* cont_Ct, __global float* cont_nt,
	__global float* cont_ut, __global float* cont_utm, __global float* cont_fc,
	__global uint* nForce_per_Segment, __global float* force_per_segment,
	__global uint* nMoment_per_Segment, __global float* moment_per_segment,
	__global double* dissip_ft, __global double* dissip_fn,
	__global double* cont_W1, __global double* cont_W2) {

	uint						iSeg1, iSeg2, iFiber1, iFiber2, i1, i2, n1, n2;
	uint						iCont, index, k;
	float						s1, s2, fn, fnElast, alpha, beta;
	float						dt, kn, kt, lambda_n, mu, omegaz;
	float						dotP1, dotP2, dotP1_ft, dotP2_ft, dotP1_fn, dotP2_fn;
	float3						e1Loc, e2Loc, utLoc, ntLoc, vi, vip1, aux, bis, Omega;
	float3						fcLoc, ftLoc, fnLoc;
	float3						dr;
	float3						fAux, sliding_Loc, shift;
	volatile __global uint* pointer;

	//--------------------------------------------------------
	iCont = (uint)get_global_id(0);

											 
											 
	if (cont_exist[iCont] != CONTACT)
		return;

	//--------------------------------------------------------
	iSeg1 = possContSeg1[iCont];
	iSeg2 = possContSeg2[iCont];
	iFiber1 = iFiberFromSegment[iSeg1];
	iFiber2 = iFiberFromSegment[iSeg2];
	i1 = iFromSegment[iSeg1];
	i2 = iFromSegment[iSeg2];
	n1 = fiber_n[iFiber1];
	n2 = fiber_n[iFiber2];
	s1 = cont_s1[iCont];
	s2 = cont_s2[iCont];

	ComputeShift(param_UInt[6], iShiftCont[iCont], param_Double[50], param_Double[51], param_Double[52], &shift);

	//--------------------------------------------------------
	utLoc = vload3(iCont, cont_utm);
	ntLoc = vload3(iCont, cont_nt);
	e1Loc = vload3(iSeg1, fiber_e);
	e2Loc = vload3(iSeg2, fiber_e);
	dt = param_Double[0];
	kn = param_Double[1];
	kt = param_Double[2];
	lambda_n = param_Double[5];
	mu = param_Double[7];

	//--------------------------------------------------------------
	//----------------- deplacement point coincidant fibre 2
	//--------------------------------------------------------------
	//---------------- vi et vip1
	vi = convert_float3(vload3(iSeg2, fiber_xt) - vload3(iSeg2, fiber_xtm)) / dt;
	vip1 = convert_float((i2 + 1) < n2) * convert_float3(vload3(iSeg2 + 1, fiber_xt) - vload3(iSeg2 + 1, fiber_xtm)) / dt;

	//----------------- calcule Omega
	if (n2 == 1)		//--- cas bille, omega est le omega de la bille
		Omega = convert_float3(vload3(iSeg2, fiber_thetat) - vload3(iSeg2, fiber_thetatm)) / dt;
	else if (i2 + 1 < n2) {				//--- cas segment normal, calculer le omega complet
		Omega = cross(e2Loc, vip1 - vi) / fiber_l0[iFiber2];	//---composante transverse
		omegaz = (fiber_thetat[3 * iSeg2 + 2] - fiber_thetatm[3 * iSeg2 + 2]) / dt;
		Omega += omegaz * e2Loc; //---composante longitudinale
	}
	else if (i2 + 1 == n2) {//------ cas bout de fibre
		omegaz = (fiber_thetat[3 * (iSeg2 - 1) + 2] - fiber_thetatm[3 * (iSeg2 - 1) + 2]) / dt;
		Omega = omegaz * vload3(iSeg2 - 1, fiber_e);
	}
	else
		printf("######## error in kernel_CalculateContactForce\n");

	//-----------------  Omega x (rc-ri)
	aux = vload3(iCont, cont_Ct) - convert_float3(vload3(iSeg2, fiber_xt));
	bis = cross(Omega, aux);

	//----------------- update le deplacement
	sliding_Loc = (vi + bis) * dt; 

	//--------------------------------------------------------------
	//----------------- deplacement point coincidant fibre 1
	//--------------------------------------------------------------
	//---------------- vi et vip1
	vi = convert_float3(vload3(iSeg1, fiber_xt) - vload3(iSeg1, fiber_xtm)) / dt;
	vip1 = convert_float((i1 + 1) < n1) * convert_float3(vload3(iSeg1 + 1, fiber_xt) - vload3(iSeg1 + 1, fiber_xtm)) / dt;

	//----------------- calcule Omega
	if (n1 == 1)		//--- cas bille, omega est le omega de la bille
		Omega = convert_float3(vload3(iSeg1, fiber_thetat) - vload3(iSeg1, fiber_thetatm)) / dt;
	else if (i1 + 1 < n1) {				//--- cas segment normal, calculer le omega complet
		Omega = cross(e1Loc, vip1 - vi) / fiber_l0[iFiber1];	//---composante transverse
		omegaz = (fiber_thetat[3 * iSeg1 + 2] - fiber_thetatm[3 * iSeg1 + 2]) / dt;
		Omega += omegaz * e1Loc;	//---composante longitudinale
	}
	else if (i1 + 1 == n1) {	//------ cas bout de fibre
		omegaz = (fiber_thetat[3 * (iSeg1 - 1) + 2] - fiber_thetatm[3 * (iSeg1 - 1) + 2]) / dt;
		Omega = omegaz * vload3(iSeg1 - 1, fiber_e);
	}
	else
		printf("######## error in kernel_CalculateContactForce\n");

	//----------------- calcule Omega x (rc-ri)
	aux = vload3(iCont, cont_Ct) - convert_float3(vload3(iSeg1, fiber_xt));
	bis = cross(Omega, aux);

	//----------------- update le deplacement
	sliding_Loc -= (vi + bis) * dt;

	//--------------------------------------------------------------
	//-------- force normale
	//--------------------------------------------------------------	
	fn = kn * cont_deltat[iCont] + lambda_n * (cont_deltat[iCont] - cont_deltatm[iCont]) / dt;

	//--------------------------------------------------------------
	//-------- deplacement tangentiel seuillé
	//--------------------------------------------------------------
	utLoc += sliding_Loc;
	utLoc -= dot(utLoc, ntLoc) * ntLoc;
	beta = length(utLoc);

	if (beta < 1.e-8) 
		utLoc = 0.;
	else {
		if (kt * beta > mu * fn){
			alpha = mu * fn / (kt * beta);
			utLoc *= alpha;
		}
	}

	//--------- force totale , dissipation et dotP
	vstore3(utLoc, iCont, cont_ut);
	fcLoc = fn * ntLoc - kt * utLoc;
	vstore3(fcLoc, iCont, cont_fc);

	beta = kn * (cont_deltat[iCont] + cont_deltatm[iCont]) / 2. + lambda_n * (cont_deltat[iCont] - cont_deltatm[iCont]) / dt;
	fnLoc = beta * ntLoc;
	ftLoc = -kt * (utLoc + vload3(iCont, cont_utm)) / 2.f;

	dotP1 = dot(fcLoc, e1Loc) * fiber_radius[iFiber1] / fiber_l0[iFiber1];
	dotP1_ft = dot(ftLoc, e1Loc) * fiber_radius[iFiber1] / fiber_l0[iFiber1];
	dotP1_fn = dot(fnLoc, e1Loc) * fiber_radius[iFiber1] / fiber_l0[iFiber1];

	dotP2 = dot(fcLoc, e2Loc) * fiber_radius[iFiber2] / fiber_l0[iFiber2];
	dotP2_ft = dot(ftLoc, e2Loc) * fiber_radius[iFiber2] / fiber_l0[iFiber2];
	dotP2_fn = dot(fnLoc, e2Loc) * fiber_radius[iFiber2] / fiber_l0[iFiber2];
	
	//-------------- force sur 1 -------------------
	if (i1 + 1 < n1) {
		pointer = &(nForce_per_Segment[iSeg1]);
		index = atomic_inc(pointer);
		dr = convert_float3(vload3(iSeg1, fiber_xt) - vload3(iSeg1, fiber_xtm));
		fAux = (1.f - s1) * (-ftLoc) - dotP1_ft * ntLoc;
		dissip_ft[iCont] -= dot(fAux, dr);
		fAux = (1.f - s1) * (-fnLoc) - dotP1_fn * ntLoc;
		dissip_fn[iCont] -= dot(fAux, dr);
		fAux = (1.f - s1) * (-fcLoc) - dotP1 * ntLoc;
		cont_W1[iCont] -= dot(fAux, dr);
		vstore3(fAux, (iSeg1)*NFPS + index, force_per_segment);

		pointer = &(nForce_per_Segment[iSeg1 + 1]);
		index = atomic_inc(pointer);
		dr = convert_float3(vload3(iSeg1+1, fiber_xt) - vload3(iSeg1+1, fiber_xtm));
		fAux = s1 * (-ftLoc) + dotP1_ft * ntLoc;
		dissip_ft[iCont] -= dot(fAux, dr);
		fAux = s1 * (-fnLoc) + dotP1_fn * ntLoc;
		dissip_fn[iCont] -= dot(fAux, dr);
		fAux = s1 * (-fcLoc) + dotP1 * ntLoc;
		cont_W1[iCont] -= dot(fAux, dr);
		vstore3(fAux, (iSeg1+1)*NFPS + index, force_per_segment);
	}
	else if (i1 + 1 == n1) {
		pointer = &(nForce_per_Segment[iSeg1]);
		index = atomic_inc(pointer);
		dr = convert_float3(vload3(iSeg1, fiber_xt) - vload3(iSeg1, fiber_xtm));
		fAux = -ftLoc;
		dissip_ft[iCont] -= dot(fAux, dr);
		fAux = -fnLoc;
		dissip_fn[iCont] -= dot(fAux, dr);	
		fAux = -fcLoc;
		cont_W1[iCont] -= dot(fAux, dr);
		vstore3(fAux, (iSeg1) * NFPS + index, force_per_segment);
	}
	else
		printf("what ?\n");

	//-------------- force sur 2 -------------------
	if (i2 + 1 < n2) {
		pointer = &(nForce_per_Segment[iSeg2]);
		index = atomic_inc(pointer);
		dr = convert_float3(vload3(iSeg2, fiber_xt) - vload3(iSeg2, fiber_xtm));
		fAux = (1.f - s2) * ftLoc - dotP2_ft * ntLoc;
		dissip_ft[iCont] -= dot(fAux, dr);
		fAux = (1.f - s2) * fnLoc - dotP2_fn * ntLoc;
		dissip_fn[iCont] -= dot(fAux, dr);
		fAux = (1.f - s2) * fcLoc - dotP2 * ntLoc;
		cont_W2[iCont] += dot(fAux, dr);
		vstore3(fAux, iSeg2 * NFPS + index, force_per_segment);

		pointer = &(nForce_per_Segment[iSeg2 + 1]);
		index = atomic_inc(pointer);
		dr = convert_float3(vload3(iSeg2+1, fiber_xt) - vload3(iSeg2+1, fiber_xtm));
		fAux = s2 * ftLoc + dotP2_ft * ntLoc;
		dissip_ft[iCont] -= dot(fAux, dr);
		fAux = s2 * fnLoc + dotP2_fn * ntLoc;
		dissip_fn[iCont] -= dot(fAux, dr);
		fAux = s2 * fcLoc + dotP2 * ntLoc;
		cont_W2[iCont] += dot(fAux, dr);
		vstore3(s2 * fcLoc + dotP2 * ntLoc, (iSeg2 + 1) * NFPS + index, force_per_segment);
	}
	else if (i2 + 1 == n2) {
		pointer = &(nForce_per_Segment[iSeg2]);
		index = atomic_inc(pointer);
		dr = convert_float3(vload3(iSeg2, fiber_xt) - vload3(iSeg2, fiber_xtm));
		fAux = ftLoc;
		dissip_ft[iCont] -= dot(fAux, dr);
		fAux = fnLoc;
		dissip_fn[iCont] -= dot(fAux, dr);
		fAux = fcLoc;
		cont_W2[iCont] += dot(fAux, dr);
		vstore3(fAux, (iSeg2) * NFPS + index, force_per_segment);
	}
	else
		printf("what ?\n");

	//-------- moments
	aux = cross(ntLoc, fcLoc);

	//-------------------------------------------------
	pointer = &(nMoment_per_Segment[iSeg1]);
	index = atomic_inc(pointer);
	if (n1 == 1) {
		for (k = 0; k < 3; k++)
			moment_per_segment[3 * (iSeg1 * NFPS + index) + k] = -fiber_radius[iFiber1] * aux[k];
	}
	else {	//---fibre
		//for (k = 0; k < 3; k++)
		//	eLoc[k] = fiber_e[3 * iSeg1 + k];
		for (k = 0; k < 3; k++)		//seul 2 importe
			moment_per_segment[3 * (iSeg1 * NFPS + index) + k] = -fiber_radius[iFiber1] * dot(aux, e1Loc);
	}

	//-------------------------------------------------
	pointer = &(nMoment_per_Segment[iSeg2]);
	index = atomic_inc(pointer);
	if (n2 == 1) {
		for (k = 0; k < 3; k++)
			moment_per_segment[3 * (iSeg2 * NFPS + index) + k] = -fiber_radius[iFiber2] * aux[k];
	}
	else {	//---fibre
		//for (k = 0; k < 3; k++)
		//	eLoc[k] = fiber_e[3 * iSeg2 + k];
		for (k = 0; k < 3; k++)		//seul 2 importe
			moment_per_segment[3 * (iSeg2 * NFPS + index) + k] = -fiber_radius[iFiber2] * dot(aux, e2Loc);
	}
}

//-----------------------------------------------------------
__kernel void kernel_AddContactForce(									//------ parameter
	__global float* fiber_f, __global float* fiber_moment,
	__global uint* nForce_per_Segment, __global float* force_per_Segment,
	__global uint* nMoment_per_Segment, __global float* moment_per_Segment) {

	uint							index, k, iSeg;
	private float3					fLoc, mLoc;

	//--------------------------------------------------------
	iSeg = (uint)get_global_id(0);

	for (k = 0; k < 3; k++) {
		fLoc[k] = 0;
		mLoc[k] = 0;
	}

	for (index = 0; index < nForce_per_Segment[iSeg]; index++) {
		for (k = 0; k < 3; k++)
			fLoc[k] += force_per_Segment[3 * (iSeg * NFPS + index) + k];
	}

	for (index = 0; index < nMoment_per_Segment[iSeg]; index++) {
		for (k = 0; k < 3; k++)
			mLoc[k] += moment_per_Segment[3 * (iSeg * NFPS + index) + k];
	}

	for (k = 0; k < 3; k++) {
		fiber_f[3 * iSeg + k] += fLoc[k];
		fiber_moment[3 * iSeg + k] += mLoc[k];
	}
}

//-----------------------------------------------------------
__kernel void kernel_ShiftContact(
	__global uint* cont_exist,										//------ cont
	__global double* cont_deltat, __global double* cont_deltatm,
	__global float* cont_ut, __global float* cont_utm) {

	uint				k, iCont;

	//---------------------------------------------------------
	iCont = (uint)(get_global_id(0));

	//----- if contact exist, do shift
	if (cont_exist[iCont] == CONTACT) {
		cont_deltatm[iCont] = cont_deltat[iCont];
		for (k = 0; k < 3; k++) {
			cont_utm[3 * iCont + k] = cont_ut[3 * iCont + k];
		}
	}
	 
								
}