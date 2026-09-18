//-----------------------------------------------------------
__kernel void kernel_Max_Displacement_Step1(
	__global uint* status,
	__global double* xt, __global double* xtm,
	__global float* partialSums, __local float* localSums) {

	uint	k;
	uint	local_id = get_local_id(0);
	uint	group_size = get_local_size(0);

	//---------- Copy from global memory to local memory
	localSums[local_id] = 0;
	if (status[get_global_id(0)] != STATUS_VIRTUAL) {
		for (k = 0; k < 3; k++)
			localSums[local_id] +=
			(xt[3 * get_global_id(0) + k] - xtm[3 * get_global_id(0) + k])
			* (xt[3 * get_global_id(0) + k] - xtm[3 * get_global_id(0) + k]);
		localSums[local_id] = sqrt(localSums[local_id]);
	}

	//----------- Loop for computing localSums
	for (uint stride = group_size / 2; stride > 0; stride /= 2) {
		// Waiting for each 2x2 addition into given workgroup
		barrier(CLK_LOCAL_MEM_FENCE);

		// Divide WorkGroup into 2 parts and add elements 2 by 2
		// between local_id and local_id + stride
		if (local_id < stride) {
			if (localSums[local_id + stride] > localSums[local_id])
				localSums[local_id] = localSums[local_id + stride];
		}
	}

	// Write result into partialSums[nWorkGroups]
	if (local_id == 0)
		partialSums[get_group_id(0)] = localSums[0];
}

//-----------------------------------------------------------
__kernel void kernel_Max_Displacement_Step2(
	__global uint* param_UInt, __global double* param_Double, __global float* reduction_devPtr) {

	uint			k;
	float			max_Deplacement = 0;

	for (k = 0; k < param_UInt[3] / (2 * SIZE_REDUCTION); k++) {
		if (reduction_devPtr[k] > max_Deplacement)
			max_Deplacement = reduction_devPtr[k];
	}

	param_Double[31] = max_Deplacement;
	param_Double[32] += 2. * max_Deplacement;
}

//-----------------------------------------------------------
__kernel void kernel_Max_Delta_Step1(
	__global uint* cont_exist,
	__global double* cont_deltat,
	__global float* partialSums, __local float* localSums) {

	uint	local_id = get_local_id(0);
	uint	group_size = get_local_size(0);

	//---------- Copy from global memory to local memory
	localSums[local_id] = 0;
	if (cont_exist[get_global_id(0)] == CONTACT) {
		localSums[local_id] = cont_deltat[get_global_id(0)];
		//printf("aa %d %e\n",local_id, localSums[local_id]);
	}

	//----------- Loop for computing localSums
	for (uint stride = group_size / 2; stride > 0; stride /= 2) {
		// Waiting for each 2x2 addition into given workgroup
		barrier(CLK_LOCAL_MEM_FENCE);

		// Divide WorkGroup into 2 parts and add elements 2 by 2
		// between local_id and local_id + stride
		if (local_id < stride) {
			if (localSums[local_id + stride] > localSums[local_id])
				localSums[local_id] = localSums[local_id + stride];
		}
	}

	// Write result into partialSums[nWorkGroups]
	if (local_id == 0) {
		partialSums[get_group_id(0)] = localSums[0];
		//	printf("bb %e\n",partialSums[get_global_id(0)]);
	}
}

//-----------------------------------------------------------
__kernel void kernel_Max_Delta_Step2(
	__global uint* param_UInt, __global double* param_Double, __global float* reduction_devPtr) {

	uint			k;
	float			max_Delta = 0;

	for (k = 0; k < param_UInt[4] / (2 * SIZE_REDUCTION); k++) {
		if (reduction_devPtr[k] > max_Delta) {
			max_Delta = reduction_devPtr[k];
		}
	}
	//	printf("reduction_devPtr %d %e\n",k, max_Delta);
	param_Double[81] = max_Delta;
}

//-----------------------------------------------------------
__kernel void kernel_Count_Contact_Step1(
	__global uint* cont_exist,
	__global uint* partialCounts, __local uint* localCounts) {

	uint	local_id = get_local_id(0);
	uint	group_size = get_local_size(0);

	//---------- Copy from global memory to local memory
	localCounts[local_id] = 0;
	if (cont_exist[get_global_id(0)] == NOCONTACT) {
		localCounts[local_id]++;
	}

	//----------- Loop for computing localSums
	for (uint stride = group_size / 2; stride > 0; stride /= 2) {
		// Waiting for each 2x2 addition into given workgroup
		barrier(CLK_LOCAL_MEM_FENCE);

		// Divide WorkGroup into 2 parts and add elements 2 by 2
		// between local_id and local_id + stride
		if (local_id < stride)
			localCounts[local_id] += localCounts[local_id + stride];
	}

	// Write result into partialCounts[nWorkGroups]
	if (local_id == 0)
		partialCounts[get_group_id(0)] = localCounts[0];
}

//-----------------------------------------------------------
__kernel void kernel_Count_Contact_Step2(
	__global uint* param_UInt, __global uint* partialCounts) {

	uint			k, countContact = 0;

	for (k = 0; k < param_UInt[4] / (2 * SIZE_REDUCTION); k++)
		countContact += partialCounts[k];

	param_UInt[7] = countContact;
}

//-----------------------------------------------------------
__kernel void kernel_Unwarp(
	__global uint* param_UInt, __global uint* status,
	__global uint* n, __global double* thetat) {

	uint					iSegment, iFiber, i, k = 2;
	double					aux;

	//-----------------------------------------
	iFiber = get_global_id(0);
	if (iFiber + 2 > param_UInt[2]) return;
	if (status[iFiber] != STATUS_FREE) return;

	for (i = 0; i < n[iFiber] - 1; i++) {
		iSegment = NSegment(n, iFiber, i);
		aux = thetat[3 * (iSegment + 1) + k] - thetat[3 * (iSegment)+k];

		if (fabs(aux - PI) < 0.1)
			thetat[3 * (iSegment + 1) + k] -= PI;
		else if (fabs(aux + PI) < 0.1)
			thetat[3 * (iSegment + 1) + k] += PI;
	}
}

//-----------------------------------------------------------
__kernel void kernel_Dissipation_Step1(
	__global double* data,
	__global float* partialCounts, __local float* localCounts) {

	uint	local_id = get_local_id(0);
	uint	group_size = get_local_size(0);

	localCounts[local_id] = data[get_global_id(0)];

	for (uint stride = group_size / 2; stride > 0; stride /= 2) {
		barrier(CLK_LOCAL_MEM_FENCE);
		if (local_id < stride)
			localCounts[local_id] += localCounts[local_id + stride];
	}

	if (local_id == 0)
		partialCounts[get_group_id(0)] = localCounts[0];
}

//-----------------------------------------------------------
__kernel void kernel_Dissipation_Step2(
	__global uint* param_UInt, __global double* param_Double, __global float* partialCounts) {

	uint			k;
	float			sum = 0;

	for (k = 0; k < param_UInt[NPARAM - 2] / SIZE_REDUCTION; k++)
		sum += partialCounts[k];

	param_Double[param_UInt[NPARAM - 1]] = sum;
}