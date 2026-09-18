//-----------------------------------------------------------
__kernel void kernel_Refresh_and_Store_Poss_Contact(
	__global uint* param_UInt,
	__global uint* oldPossContSeg1, __global uint* possContSeg1,
	__global uint* oldPossContSeg2, __global uint* possContSeg2,
	__global uint* oldContExist, __global uint* contExist,
	__global double* old_cont_deltatm, __global double* cont_deltatm,
	__global float* old_cont_utm, __global float* cont_utm) {

	uint k, iCont;

	iCont = (uint)(get_global_id(0));

	if (param_UInt[0] != 0) {
		oldPossContSeg1[iCont] = possContSeg1[iCont];
		oldPossContSeg2[iCont] = possContSeg2[iCont];
		oldContExist[iCont] = contExist[iCont];
		old_cont_deltatm[iCont] = cont_deltatm[iCont];
		for (k = 0; k < 3; k++)
			old_cont_utm[3 * iCont + k] = cont_utm[3 * iCont + k];
	}
	else {
		oldPossContSeg1[iCont] = UINT_MAX;
		oldPossContSeg2[iCont] = UINT_MAX;
		oldContExist[iCont] = NOPOSSCONTACT;
	}

	possContSeg1[iCont] = UINT_MAX;
	possContSeg2[iCont] = UINT_MAX;
	contExist[iCont] = NOPOSSCONTACT;

	return;
}

//-----------------------------------------------------------
__kernel void kernel_Poss_Contact(
	__global uint* param_UInt, __global double* param_Double,
	__global uint* n, __global uint* status, __global float* radius, __global double* xt,
	__global uint* iFiberFromSegment, __global uint* iFromSegment,
	__global uint* n_possCont_in_wg,
	__global uint* possContSeg1, __global uint* possContSeg2,
	__global uint* iShiftCont, __global uint* cont_exist,
	__local float* x1, __local float* x2) {

	uint					iAux, iGroup, index, j, iShift;
	uint					iSeg1, iFiber1, i1, iSeg2, iFiber2, i2, t1, t2;
	float3					a = { 0.,0.,0. }, b = { 0.,0.,0. }, c = { 0.,0.,0. }, shift = { 0.,0.,0. };
	float					aa, bb, cc, ab, bc, ac, aux, s1, s2, dist, delta, l1, l2, epsStar;
	volatile __global uint* counterPtr = 0;

	//-------- retrieve segment
	iSeg1 = (uint)(get_global_id(0));
	iSeg2 = (uint)(get_global_id(1));

	//-------- read in local memory if necessary
	if (iSeg2 > iSeg1) {

		//---------------------------------------
		iFiber1 = iFiberFromSegment[iSeg1];
		iFiber2 = iFiberFromSegment[iSeg2];
		i1 = iFromSegment[iSeg1];
		i2 = iFromSegment[iSeg2];
		t1 = get_local_id(0);
		t2 = get_local_id(1);

		//---------------------------------------
		for (j = 0; j < 3; j++)		x1[3 * t1 + j] = xt[3 * iSeg1 + j];
		if (i1 + 2 <= n[iFiber1])
			for (j = 0; j < 3; j++)		x1[3 * (t1 + 1) + j] = xt[3 * (iSeg1 + 1) + j];
		for (j = 0; j < 3; j++)		x2[3 * t2 + j] = xt[3 * iSeg2 + j];
		if (i2 + 2 <= n[iFiber2])
			for (j = 0; j < 3; j++)		x2[3 * (t2 + 1) + j] = xt[3 * (iSeg2 + 1) + j];
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	//--------------------------------
	if (iSeg2 <= iSeg1)
		return;
	if ((status[iFiber1] == STATUS_VIRTUAL) || (status[iFiber2] == STATUS_VIRTUAL))
		return;
	if (FastCheck(iFiber1, i1, iFiber2, i2) == SKIP_CONTACT_DETECTION)
		return;

	//---------------------------------
	iShift = (uint)(get_global_id(2));
	//	if (iShift != 0) printf("Arghh kernel_Poss_Contact\n");
	ComputeShift(param_UInt[6], iShift, (float)param_Double[50], (float)param_Double[51], (float)param_Double[52], &shift);

	//-------------- retrieve group ID
	if (2 * (iSeg2 + 1) > param_UInt[3])
		iAux = param_UInt[3] - iSeg2 - 1;
	else
		iAux = iSeg2;

	iGroup = (2 * iAux * (param_UInt[4] / (WG * 32))) / (param_UInt[3] / 32);
	counterPtr = &(n_possCont_in_wg[iGroup]);

	//------------------------------------
	for (j = 0; j < 3; j++) {
		a[j] = x2[3 * t2 + j] + shift[j] - x1[3 * t1 + j];
		if (i1 + 2 <= n[iFiber1])
			b[j] = x1[3 * (t1 + 1) + j] - x1[3 * t1 + j];
		if (i2 + 2 <= n[iFiber2])
			c[j] = x2[3 * (t2 + 1) + j] - x2[3 * t2 + j];
	}
	aa = dot(a, a), bb = dot(b, b), cc = dot(c, c);
	ab = dot(a, b), ac = dot(a, c), bc = dot(b, c);
	l1 = sqrt(bb), l2 = sqrt(cc);
	epsStar = param_Double[30];

	//------------ test contact Cylinder-Cylinder
	if ((i1 + 2 <= n[iFiber1]) && (i2 + 2 <= n[iFiber2])) {
		aux = bb * cc - bc * bc;
		if (aux > EPS) {
			s1 = (cc * ab - ac * bc) / aux;
			s2 = -(bb * ac - ab * bc) / aux;
			if ((s1 > (-epsStar / l1)) && (s1 < (1 + epsStar / l1))
				&& (s2 > (-epsStar / l2)) && (s2 < (1 + epsStar / l2))) {
				dist = sqrt(fmax(0, (aa - 2. * ab * s1 + 2 * ac * s2 - 2 * bc * s1 * s2 + bb * s1 * s1 + cc * s2 * s2)));
				delta = radius[iFiber1] + radius[iFiber2] - dist;
				if (delta > -epsStar) {
					index = atomic_inc(counterPtr);
					possContSeg1[index + iGroup * WG] = iSeg1;
					possContSeg2[index + iGroup * WG] = iSeg2;
					iShiftCont[index + iGroup * WG] = iShift;
					cont_exist[index + iGroup * WG] = NOCONTACT;
					return;
				}
			}
		}
	}

	//------------ test contact Cylinder-Sphere
	if (i1 + 2 <= n[iFiber1]) {
		s1 = ab / bb;
		if ((s1 > (-epsStar / l1)) && (s1 < (1 + epsStar / l1))) {
			dist = sqrt(fmax(0, (aa - 2. * ab * s1 + bb * s1 * s1)));
			delta = radius[iFiber1] + radius[iFiber2] - dist;
			if (delta > -epsStar) {
				index = atomic_inc(counterPtr);
				possContSeg1[index + iGroup * WG] = iSeg1;
				possContSeg2[index + iGroup * WG] = iSeg2;
				iShiftCont[index + iGroup * WG] = iShift;
				cont_exist[index + iGroup * WG] = NOCONTACT;
				return;
			}
		}
	}

	//------------ test contact Sphere-Cylinder
	if (i2 + 2 <= n[iFiber2]) {
		s2 = -ac / cc;
		if ((s2 > (-epsStar / l2)) && (s2 < (1 + epsStar / l2))) {
			dist = sqrt(fmax(0, (aa + 2 * ac * s2 + cc * s2 * s2)));
			delta = radius[iFiber1] + radius[iFiber2] - dist;
			if (delta > -epsStar) {
				index = atomic_inc(counterPtr);
				possContSeg1[index + iGroup * WG] = iSeg1;
				possContSeg2[index + iGroup * WG] = iSeg2;
				iShiftCont[index + iGroup * WG] = iShift;
				cont_exist[index + iGroup * WG] = NOCONTACT;
				return;
			}
		}
	}

	//------------ test contact Sphere-Sphere
	dist = sqrt(fmax(0, (aa)));
	delta = radius[iFiber1] + radius[iFiber2] - dist;
	if (delta > -epsStar) {
		index = atomic_inc(counterPtr);
		possContSeg1[index + iGroup * WG] = iSeg1;
		possContSeg2[index + iGroup * WG] = iSeg2;
		iShiftCont[index + iGroup * WG] = iShift;
		cont_exist[index + iGroup * WG] = NOCONTACT;
		return;
	}
}

//-----------------------------------------------------------
__kernel void kernel_Align_Poss_Contact(
	__global uint* n_possCont_in_wg,
	__global uint* oldPossContactSeg1, __global uint* possContactSeg1,
	__global uint* oldPossContactSeg2, __global uint* possContactSeg2,
	__global uint* iShiftCont,
	__global uint* oldContExist, __global uint* contExist,
	__global double* old_cont_deltatm, __global double* cont_deltatm,
	__global float* old_cont_utm, __global float* cont_utm) {


	uint							k, iGroup, iContOld, iContNew, found;

	iGroup = (uint)(get_global_id(0));

	//---------------------------------------------------------------
	for (iContNew = 0; contExist[iContNew + iGroup * WG] != NOPOSSCONTACT; iContNew++) {
		for (iContOld = 0; oldContExist[iContOld + iGroup * WG] != NOPOSSCONTACT; iContOld++) {
			if ((oldPossContactSeg1[iContOld + iGroup * WG] == possContactSeg1[iContNew + iGroup * WG])
				&& (oldPossContactSeg2[iContOld + iGroup * WG] == possContactSeg2[iContNew + iGroup * WG])) {

				contExist[iContNew + iGroup * WG] = oldContExist[iContOld + iGroup * WG];
				cont_deltatm[iContNew + iGroup * WG] = old_cont_deltatm[iContOld + iGroup * WG];
				for (k = 0; k < 3; k++)
					cont_utm[3 * (iContNew + iGroup * WG) + k] = old_cont_utm[3 * (iContOld + iGroup * WG) + k];
				break;
			}
		}
	}

	//------------------- check de pas d'erreur
	for (iContOld = 0; oldContExist[iContOld + iGroup * WG] != NOPOSSCONTACT; iContOld++) {

		if (oldContExist[iContOld + iGroup * WG] == NOCONTACT) continue;

		found = 0;
		for (iContNew = 0; contExist[iContNew + iGroup * WG] != NOPOSSCONTACT; iContNew++) {
			if ((possContactSeg1[iContNew + iGroup * WG] == oldPossContactSeg1[iContOld + iGroup * WG])
				&& (possContactSeg2[iContNew + iGroup * WG] == oldPossContactSeg2[iContOld + iGroup * WG])) {
				found = 1;
				break;
			}
		}
		if (found == 0)
			printf("###### error in kernel_Align_Poss_Contact\n");
	}
}