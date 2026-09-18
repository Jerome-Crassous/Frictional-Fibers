#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <CL/cl.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

#define MAX_SOURCE_SIZE (0x050000)
#define MAX_FILE_SIZE (0x010000)
#define CHUNK_SIZE 1024
#define __CL_ENABLE_EXCEPTIONS

cl_device_id device_id;
cl_context context;
cl_command_queue command_queue;

cl_kernel kernel_Calculate_le, kernel_BendingForce, kernel_TwistForce;
cl_kernel kernel_Refresh_and_Store_Poss_Contact, kernel_Poss_Contact, kernel_Align_Poss_Contact;
cl_kernel kernel_Contact, kernel_RemoveDoubleContact, kernel_CalculateContactForce, kernel_AddContactForce, kernel_ShiftContact;
cl_kernel kernel_Integrate_and_Shift, kernel_Compute_m1bar, kernel_Compute_m1;
cl_kernel kernel_Max_Delta_Step1, kernel_Max_Delta_Step2;
cl_kernel kernel_Count_Contact_Step1, kernel_Count_Contact_Step2;
cl_kernel kernel_Unwarp;
cl_kernel kernel_Dissipation_Step1, kernel_Dissipation_Step2;
cl_kernel kernel_Max_Displacement_Step1, kernel_Max_Displacement_Step2;
cl_kernel kernel_Specific_OneTime, kernel_Specific_Position, kernel_Specific_Force;

//---------------------- hostPtr --------------------------------------
extern cl_uint* param_UInt_hostPtr;
extern cl_double* param_Double_hostPtr;
extern cl_uint* n_possCont_in_wg_hostPtr;						//[NcontactMax/WG]
extern cl_uint* possContID_hostPtr;							//[NcontactMax]
extern cl_uint* nForce_per_Segment_hostPtr;                     //NSegment
extern cl_uint* nMoment_per_Segment_hostPtr;                    //NSegment

//---------------------- devPtr --------------------------------------
extern cl_mem		param_UInt_devPtr, param_Double_devPtr;
extern cl_mem		fiber_n_devPtr, fiber_status_devPtr;
extern cl_mem		fiber_radius_devPtr, fiber_l0_devPtr, fiber_k0_devPtr;
extern cl_mem		fiber_mass_devPtr, fiber_b_devPtr, fiber_c_devPtr, fiber_j0_devPtr;
extern cl_mem		fiber_flag_devPtr; 
extern cl_mem		fiber_lt_devPtr, fiber_ltm_devPtr;
extern cl_mem		fiber_kappa1_bar_devPtr, fiber_kappa2_bar_devPtr;
extern cl_mem		fiber_xt_devPtr, fiber_xtm_devPtr;
extern cl_mem		fiber_thetat_devPtr, fiber_thetatm_devPtr;
extern cl_mem		fiber_f_devPtr, fiber_e_devPtr;
extern cl_mem		fiber_moment_devPtr;
extern cl_mem		fiber_m1_bar_devPtr, fiber_m1_devPtr;
extern cl_mem		fiber_dissip_visc_global_devPtr, fiber_dissip_visc_stretch_devPtr, fiber_Wop_devPtr;

extern cl_mem		iFiberFromSegment_devPtr, iFromSegment_devPtr;

extern cl_mem		n_possCont_in_wg_devPtr;
extern cl_mem		possContSeg1_devPtr, oldPossContSeg1_devPtr;
extern cl_mem		possContSeg2_devPtr, oldPossContSeg2_devPtr;
extern cl_mem		iShiftCont_devPtr;

extern cl_mem       cont_exist_devPtr, old_cont_exist_devPtr;
extern cl_mem       cont_s1_devPtr, cont_s2_devPtr;
extern cl_mem       cont_deltat_devPtr, cont_deltatm_devPtr, old_cont_deltatm_devPtr;
extern cl_mem       cont_Ct_devPtr;
extern cl_mem       cont_nt_devPtr;
extern cl_mem       cont_ut_devPtr, cont_utm_devPtr, old_cont_utm_devPtr, cont_fc_devPtr;
extern cl_mem		cont_dissip_ft_devPtr, cont_dissip_fn_devPtr;
extern cl_mem		cont_W1_devPtr, cont_W2_devPtr;
extern cl_mem		nForce_per_Segment_devPtr, nMoment_per_Segment_devPtr;      //NSEG
extern cl_mem		force_per_Segment_devPtr, moment_per_Segment_devPtr;        //NSEG*NFPS

extern cl_mem		reduction_devPtr, partialCount_devPtr;

void Execute_Kernel_Sum_One_Array_In_Slot(cl_mem dataBuffer, size_t size, unsigned int slot);

//-------------------------------------------------------
void InitDevice(void) {

	int         i, i0 = 0;
	cl_platform_id platform_id[10];
	cl_uint ret_num_devices;
	cl_uint num_platforms;
	cl_char string[WG * 10] = { 0 };
	cl_int ret;

	printf("input InitDevice\n");

	//-------------  Get platform and device information
	ret = clGetPlatformIDs(0, NULL, &num_platforms);
	checkError(ret, "clGetPlatformIDs");

	ret = clGetPlatformIDs(num_platforms, platform_id, NULL);
	checkError(ret, "clGetPlatformIDs");

	for (i = 0; i < num_platforms; i++) {
		ret = clGetPlatformInfo(platform_id[i], CL_PLATFORM_NAME, sizeof(string), &string, NULL);
		checkError(ret, "clGetPlatformInfo");
		printf("Platform %d: %s\n", i, string);
		if ((string[0] == 'N') && (string[1] == 'V'))
			i0 = i;
	}

	i = i0;

	ret = clGetPlatformInfo(platform_id[i], CL_PLATFORM_NAME, sizeof(string), &string, NULL);
	checkError(ret, "clGetPlatformInfo");
	printf("Platform %d: %s\n", i, string);

	// Create device
	device_id = NULL;
	ret = clGetDeviceIDs(platform_id[i], CL_DEVICE_TYPE_ALL, 1, &device_id, &ret_num_devices);
	checkError(ret, "clGetDeviceIDs");

	// Create an OpenCL context
	context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
	checkError(ret, "clCreateContext");

	// Create a command queue
	command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
	checkError(ret, "clCreateCommandQueue");

	// Get device name
	ret = clGetDeviceInfo(device_id, CL_DEVICE_OPENCL_C_VERSION, sizeof(string), &string, NULL);
	checkError(ret, "Getting device name");
	printf("\t\tCL_DEVICE_OPENCL_C_VERSION : %s\n", string);

	return;
}

//-------------------------------------------------------
void Create_Kernels(void) {

	FILE* fp;
	char* source_str;
	size_t read = 0, source_size = 0;

	cl_int			ret;
	cl_program program;
	unsigned int        j;

	printf("input Create_Kernels\n");

	source_str = (char*)GlobalAlloc(0, MAX_SOURCE_SIZE);

	//-------------  Load the kernel source code into the array source_str
	fp = fopen("..\\fiberLib_Common_Macros_v9.6.h", "r");

	if (!fp) {
		fprintf(stderr, "Failed to load fiberLib_Common_Macros_v9.6.h\n");
		exit(1);
	}
	while ((read = fread(source_str + source_size, 1, CHUNK_SIZE, fp)) > 0)
		source_size += read;

	fclose(fp);
	printf("la source fiberLib_Common_Macros_v9.6.h est %llu %llu\n", read, source_size);

	//-------------
	fp = fopen("..\\functions_v9.6.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load functions_v9.6.cl.\n");
		exit(1);
	}
	read = fread_s(source_str + source_size, MAX_SOURCE_SIZE, 1, MAX_FILE_SIZE, fp);
	source_size += read;
	fclose(fp);
	printf("la source functions_v9.6.cl est %llu %llu\n", read, source_size);

	//-------------
	fp = fopen("kernel_specific.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel_specific.cl\n");
		exit(1);
	}
	read = fread_s(source_str + source_size, MAX_SOURCE_SIZE, 1, MAX_FILE_SIZE, fp);
	source_size += read;
	fclose(fp);
	printf("la source kernel_specific.cl est %llu %llu\n", read, source_size);

	//-------------
	fp = fopen("..\\kernel_max_v9.6.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel_max_v9.6.cl\n");
		exit(1);
	}
	read = fread_s(source_str + source_size, MAX_SOURCE_SIZE, 1, MAX_FILE_SIZE, fp);
	source_size += read;
	fclose(fp);
	printf("la source kernel_max_v9.6.cl est %llu %llu\n", read, source_size);

	//-------------
	fp = fopen("..\\kernel_le_force_v9.6.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel_le_force_v9.6.cl\n");
		exit(1);
	}
	read = fread_s(source_str + source_size, MAX_SOURCE_SIZE, 1, MAX_FILE_SIZE, fp);
	source_size += read;
	fclose(fp);
	printf("la source kernel_le_force_v9.6.cl est %llu %llu\n", read, source_size);

	//-------------
	fp = fopen("..\\kernel_poss_contact_v9.6.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel_poss_contact_v9.6.cl\n");
		exit(1);
	}
	read = fread_s(source_str + source_size, MAX_SOURCE_SIZE, 1, MAX_FILE_SIZE, fp);
	source_size += read;
	fclose(fp);
	printf("la source kernel_poss_contact_v9.6.cl est %llu %llu\n", read, source_size);

	//-------------
	fp = fopen("..\\kernel_contact_v9.6.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel_contact_v9.6.cl\n");
		exit(1);
	}
	read = fread_s(source_str + source_size, MAX_SOURCE_SIZE, 1, MAX_FILE_SIZE, fp);
	source_size += read;
	fclose(fp);
	printf("la source kernel_contact_v9.6.cl est %llu %llu\n", read, source_size);

	//-------------
	fp = fopen("..\\kernel_integrate_shift_v9.6.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel_integrate_shift_v9.6.cl\n");
		exit(1);
	}
	read = fread_s(source_str + source_size, MAX_SOURCE_SIZE, 1, MAX_FILE_SIZE, fp);
	source_size += read;
	fclose(fp);
	printf("la source kernel_integrate_shift_v9.6.cl est %llu %llu\n", read, source_size);

	//-------------- Create a program from the kernel source
	program = clCreateProgramWithSource(context, 1, (const char**)&source_str, (const size_t*)&source_size, &ret);
	checkError(ret, "clCreateProgramWithSource");

	ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
	checkError(ret, "clBuildProgram");

	char* build_log; size_t log_size;
	ret = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
	checkError(ret, "clGetProgramBuildInfo");
	build_log = (char*)malloc((log_size + 1));

	// Second call to get the log
	ret = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL);
	checkError(ret, "clGetProgramBuildInfo");
	build_log[log_size] = '\0';
	printf("--- Build log ---\n ");
	fprintf(stderr, "%s\n", build_log);
	free(build_log);


	//------------- clCreateKernel kernel_Calculate_le
	kernel_Calculate_le = clCreateKernel(program, "kernel_Calculate_le", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_k0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_l0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_lt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_ltm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_e_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_thetatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_dissip_visc_global_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Calculate_le, j++, sizeof(cl_mem), (void*)&fiber_dissip_visc_stretch_devPtr); checkError(ret, "clSetKernelArg");

	//------------- clCreateKernel kernel_BendingForce
	kernel_BendingForce = clCreateKernel(program, "kernel_BendingForce", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_b_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_l0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_kappa1_bar_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_kappa2_bar_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_e_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_BendingForce, j++, sizeof(cl_mem), (void*)&fiber_m1_devPtr); checkError(ret, "clSetKernelArg");

	//------------- clCreateKernel twist_Force
	kernel_TwistForce = clCreateKernel(program, "kernel_TwistForce", &ret);
	checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_c_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_l0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_e_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_m1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_TwistForce, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");

	//------------- clCreateKernel kernel_Refresh_and_Store_Poss_Contact
	kernel_Refresh_and_Store_Poss_Contact = clCreateKernel(program, "kernel_Refresh_and_Store_Poss_Contact", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&oldPossContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&possContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&oldPossContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&possContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&old_cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&old_cont_deltatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&cont_deltatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&old_cont_utm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Refresh_and_Store_Poss_Contact, j++, sizeof(cl_mem), (void*)&cont_utm_devPtr); checkError(ret, "clSetKernelArg");

	//------------- clCreateKernel kernel_Poss_Contact
	kernel_Poss_Contact = clCreateKernel(program, "kernel_Poss_Contact", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&fiber_radius_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&n_possCont_in_wg_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&possContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&possContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&iShiftCont_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, 3 * (WG2D + 1) * sizeof(cl_float), NULL); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Poss_Contact, j++, 3 * (WG2D + 1) * sizeof(cl_float), NULL); checkError(ret, "clSetKernelArg");

	//------------- clCreateKernel kernel_Align_Poss_Contact
	kernel_Align_Poss_Contact = clCreateKernel(program, "kernel_Align_Poss_Contact", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&n_possCont_in_wg_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&oldPossContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&possContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&oldPossContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&possContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&iShiftCont_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&old_cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&old_cont_deltatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&cont_deltatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&old_cont_utm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Align_Poss_Contact, j++, sizeof(cl_mem), (void*)&cont_utm_devPtr); checkError(ret, "clSetKernelArg");

	//-------------------------------- kernel_contact
	kernel_Contact = clCreateKernel(program, "kernel_Contact", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&fiber_radius_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&possContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&possContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&iShiftCont_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_s1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_s2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_deltat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_deltatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_Ct_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_nt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_ut_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Contact, j++, sizeof(cl_mem), (void*)&cont_utm_devPtr); checkError(ret, "clSetKernelArg");

	//-------------------------------- kernel_RemoveDoubleContact
	kernel_RemoveDoubleContact = clCreateKernel(program, "kernel_RemoveDoubleContact", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_RemoveDoubleContact, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_RemoveDoubleContact, j++, sizeof(cl_mem), (void*)&possContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_RemoveDoubleContact, j++, sizeof(cl_mem), (void*)&possContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_RemoveDoubleContact, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_RemoveDoubleContact, j++, sizeof(cl_mem), (void*)&cont_deltat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_RemoveDoubleContact, j++, sizeof(cl_mem), (void*)&cont_Ct_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_CalculateContactForce
	kernel_CalculateContactForce = clCreateKernel(program, "kernel_CalculateContactForce", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_radius_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_l0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_thetatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&fiber_e_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&possContSeg1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&possContSeg2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&iShiftCont_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_s1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_s2_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_deltat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_deltatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_Ct_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_nt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_ut_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_utm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_fc_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&nForce_per_Segment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&force_per_Segment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&nMoment_per_Segment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&moment_per_Segment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_dissip_ft_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_dissip_fn_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_W1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_CalculateContactForce, j++, sizeof(cl_mem), (void*)&cont_W2_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------- kernel_AddContactForce
	kernel_AddContactForce = clCreateKernel(program, "kernel_AddContactForce", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_AddContactForce, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_AddContactForce, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_AddContactForce, j++, sizeof(cl_mem), (void*)&nForce_per_Segment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_AddContactForce, j++, sizeof(cl_mem), (void*)&force_per_Segment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_AddContactForce, j++, sizeof(cl_mem), (void*)&nMoment_per_Segment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_AddContactForce, j++, sizeof(cl_mem), (void*)&moment_per_Segment_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_ShiftContact
	kernel_ShiftContact = clCreateKernel(program, "kernel_ShiftContact", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_ShiftContact, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_ShiftContact, j++, sizeof(cl_mem), (void*)&cont_deltat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_ShiftContact, j++, sizeof(cl_mem), (void*)&cont_deltatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_ShiftContact, j++, sizeof(cl_mem), (void*)&cont_ut_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_ShiftContact, j++, sizeof(cl_mem), (void*)&cont_utm_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Integrate_and_Shift
	kernel_Integrate_and_Shift = clCreateKernel(program, "kernel_Integrate_and_Shift", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_mass_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_j0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_thetatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Integrate_and_Shift, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Compute_m1bar
	kernel_Compute_m1bar = clCreateKernel(program, "kernel_Compute_m1bar", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_l0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_mass_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_j0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_lt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_e_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_thetatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_m1_bar_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_m1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1bar, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Compute_m1
	kernel_Compute_m1 = clCreateKernel(program, "kernel_Compute_m1", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_l0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_mass_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_j0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_lt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_e_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_thetatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_m1_bar_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_m1_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Compute_m1, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Max_Displacement
	kernel_Max_Displacement_Step1 = clCreateKernel(program, "kernel_Max_Displacement_Step1", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Max_Displacement_Step1, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Displacement_Step1, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Displacement_Step1, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Displacement_Step1, j++, sizeof(cl_mem), (void*)&reduction_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Displacement_Step1, j++, (2 * SIZE_REDUCTION) * sizeof(cl_float), NULL); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Max_Displacement_Step2
	kernel_Max_Displacement_Step2 = clCreateKernel(program, "kernel_Max_Displacement_Step2", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Max_Displacement_Step2, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Displacement_Step2, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Displacement_Step2, j++, sizeof(cl_mem), (void*)&reduction_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Count_Contact_Step1
	kernel_Count_Contact_Step1 = clCreateKernel(program, "kernel_Count_Contact_Step1", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Count_Contact_Step1, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Count_Contact_Step1, j++, sizeof(cl_mem), (void*)&partialCount_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Count_Contact_Step1, j++, (2 * SIZE_REDUCTION) * sizeof(cl_uint), NULL); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Count_Contact_Step2
	kernel_Count_Contact_Step2 = clCreateKernel(program, "kernel_Count_Contact_Step2", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Count_Contact_Step2, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Count_Contact_Step2, j++, sizeof(cl_mem), (void*)&partialCount_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Max_Delta_Step1
	kernel_Max_Delta_Step1 = clCreateKernel(program, "kernel_Max_Delta_Step1", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Max_Delta_Step1, j++, sizeof(cl_mem), (void*)&cont_exist_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Delta_Step1, j++, sizeof(cl_mem), (void*)&cont_deltat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Delta_Step1, j++, sizeof(cl_mem), (void*)&reduction_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Delta_Step1, j++, (2 * SIZE_REDUCTION) * sizeof(cl_float), NULL); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Max_Delta_Step2
	kernel_Max_Delta_Step2 = clCreateKernel(program, "kernel_Max_Delta_Step2", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Max_Delta_Step2, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Delta_Step2, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Max_Delta_Step2, j++, sizeof(cl_mem), (void*)&reduction_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Unwarp
	kernel_Unwarp = clCreateKernel(program, "kernel_Unwarp", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Unwarp, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Unwarp, j++, sizeof(cl_mem), (void*)&fiber_status_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Unwarp, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Unwarp, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");


	//------------------------ kernel_Dissipation_Step1
	kernel_Dissipation_Step1 = clCreateKernel(program, "kernel_Dissipation_Step1", &ret); checkError(ret, "clCreateKernel");
	j = 1;	//############### l'argument 0 est definie dans Execute_Kernel_Dissipation
	ret = clSetKernelArg(kernel_Dissipation_Step1, j++, sizeof(cl_mem), (void*)&reduction_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Dissipation_Step1, j++, (2 * SIZE_REDUCTION) * sizeof(cl_float), NULL); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Dissipation_Step2
	kernel_Dissipation_Step2 = clCreateKernel(program, "kernel_Dissipation_Step2", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Dissipation_Step2, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Dissipation_Step2, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Dissipation_Step2, j++, sizeof(cl_mem), (void*)&reduction_devPtr); checkError(ret, "clSetKernelArg");

	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//------------------------ kernel_Specific_OneTime
	kernel_Specific_OneTime = clCreateKernel(program, "kernel_Specific_OneTime", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Specific_OneTime, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_OneTime, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Specific_Position
	kernel_Specific_Position = clCreateKernel(program, "kernel_Specific_Position", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_radius_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_l0_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_flag_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_e_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_thetatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_m1_bar_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Position, j++, sizeof(cl_mem), (void*)&fiber_m1_devPtr); checkError(ret, "clSetKernelArg");

	//------------------------ kernel_Specific_Force
	kernel_Specific_Force = clCreateKernel(program, "kernel_Specific_Force", &ret); checkError(ret, "clCreateKernel");
	j = 0;
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&param_UInt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&param_Double_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_n_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&iFiberFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&iFromSegment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_flag_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_xtm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_xt_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_thetatm_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_thetat_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_f_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_moment_devPtr); checkError(ret, "clSetKernelArg");
	ret = clSetKernelArg(kernel_Specific_Force, j++, sizeof(cl_mem), (void*)&fiber_Wop_devPtr); checkError(ret, "clSetKernelArg");
}

//-------------------------------------------------------
void Execute_Kernel_Calculate_le(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG };

	//---------------------------
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Calculate_le, 1, NULL, global_Size, local_Size, 0, NULL, NULL);checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_Bending_Force(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG / 4 };

	//---------------------------;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_BendingForce, 1, NULL, global_Size, local_Size, 0, NULL, NULL);checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_Twist_Force(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG / 2 };

	//---------------------------;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_TwistForce, 1, NULL, global_Size, local_Size, 0, NULL, NULL);checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_Specific_Force(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG };

	//---------------------------;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Specific_Force, 1, NULL, global_Size, local_Size, 0, NULL, NULL);checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_PossibleContact(parameter* parameterPtr) {

	int			        nContactMax = parameterPtr->nContactMax;
	unsigned int		iGroup, maxPossibleContacts = 0, iGroupMaxPossibleContacts;
	size_t              global_Size[1];
	size_t              local_Size[1];
	cl_int              ret;
	cl_double			integratedDispacement, epsStar;

	//--------------- test si il faut appeler le kernel
	if (parameterPtr->flag_Reset_Contact == 1) {
		integratedDispacement = 1.e10;    //set integrated displacement to a large value pour forcer le reset
		ret = clEnqueueWriteBuffer(command_queue, param_UInt_devPtr, CL_TRUE, 39 * sizeof(cl_uint), sizeof(cl_uint), &(parameterPtr->flag_Reset_Contact), 0, NULL, NULL); 
		checkError(ret, "clEnqueueWriteBuffer");
		parameterPtr->flag_Reset_Contact = 0;
	}
	else {		//-- read integratedDispacement
		ret = clEnqueueReadBuffer(command_queue, param_Double_devPtr, CL_TRUE, 32 * sizeof(cl_double), sizeof(cl_double), &integratedDispacement, 0, NULL, NULL);
		checkError(ret, "clEnqueueReadBuffer");
	}

	if (integratedDispacement < parameterPtr->epsStar)  //si deplacement integré trop petit, ne rien faire
		return;

	//--------------- reset le deplacement integré
	integratedDispacement = 0.;
	ret = clEnqueueWriteBuffer(command_queue, param_Double_devPtr, CL_TRUE, 32 * sizeof(cl_double), sizeof(cl_double), &integratedDispacement, 0, NULL, NULL);
	checkError(ret, "clEnqueueWriteBuffer");

	//---------------- store and reset possible contact  
	global_Size[0] = parameterPtr->nContactMax;
	local_Size[0] = WG;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Refresh_and_Store_Poss_Contact, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	//--------------- reset liste de potentiel contact
	for (iGroup = 0; iGroup < nContactMax / WG; iGroup++)
		n_possCont_in_wg_hostPtr[iGroup] = 0;
	ret = clEnqueueWriteBuffer(command_queue, n_possCont_in_wg_devPtr, CL_TRUE, 0, nContactMax / WG * sizeof(cl_uint), n_possCont_in_wg_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueWriteBuffer");

	//--------------- execute PossibleContact
	size_t	global_Size3D[3] = {parameterPtr->nSegment, parameterPtr->nSegment, 1};
	size_t	local_Size3D[3] = {WG2D, WG2D, 1};

	if (parameterPtr->periodic == PERIODIC_NO)
		global_Size3D[2] = 1;
	else if (parameterPtr->periodic == PERIODIC_XY)
		global_Size3D[2] = 9;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Poss_Contact, 3, NULL, global_Size3D, local_Size3D, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	//---------------- lit le nombre de potentiel contact par group et determine le nombre max de potentiel contact	
	ret = clEnqueueReadBuffer(command_queue, n_possCont_in_wg_devPtr, CL_TRUE, 0, nContactMax / WG * sizeof(cl_uint), n_possCont_in_wg_hostPtr, 0, NULL, NULL);
	checkError(ret, "clEnqueueReadBuffer");

	maxPossibleContacts = 0;
	iGroupMaxPossibleContacts = 0;
	for (iGroup = 0; iGroup < nContactMax / WG; iGroup++) {
		if (n_possCont_in_wg_hostPtr[iGroup] > maxPossibleContacts) {
			maxPossibleContacts = n_possCont_in_wg_hostPtr[iGroup];
			iGroupMaxPossibleContacts = iGroup;
		}
	}

	printf("Execute_Kernel_PossibleContact : iter=%u maxPossibleContacts = %d/%d iGroup=%u eps=%e\n", 
		parameterPtr->iter, maxPossibleContacts, WG, iGroupMaxPossibleContacts, parameterPtr->epsStar);

		//--------------- ajuste epsStar
	if (maxPossibleContacts > 3 * WG / 4)
		parameterPtr->epsStar -= 0.02;
	else
		parameterPtr->epsStar += 0.01;

	if (parameterPtr->epsStar < 0.01)
		parameterPtr->epsStar = 0.01;

	epsStar = parameterPtr->epsStar;
	ret = clEnqueueWriteBuffer(command_queue, param_Double_devPtr, CL_TRUE, 30 * sizeof(cl_double), sizeof(cl_double), &epsStar, 0, NULL, NULL);
	checkError(ret, "clEnqueueWriteBuffer");

	//--------------- execute kernel_Align_Poss_Contact
	global_Size[0] = parameterPtr->nContactMax / WG;
	local_Size[0] = 1;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Align_Poss_Contact, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_Contact(parameter* parameterPtr) {

	size_t global_Size[1] = { parameterPtr->nContactMax };
	size_t local_Size[1] = { WG / 2 };
	cl_int			ret;

	//--------------- execute kernel_Contact
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Contact, 1, NULL, global_Size, local_Size, 0, NULL, NULL);checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_RemoveDoubleContact(parameter* parameterPtr) {

	size_t global_Size[1] = { parameterPtr->nContactMax };
	size_t local_Size[1] = { WG };
	cl_int			ret;

	//--------------- execute kernel_Contact
	ret = clEnqueueNDRangeKernel(command_queue, kernel_RemoveDoubleContact, 1, NULL, global_Size, local_Size, 0, NULL, NULL);checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_CalculateContactForce(parameter* parameterPtr) {

	size_t global_Size[1] = { parameterPtr->nContactMax };
	size_t local_Size[1] = { WG / 2 };		//####### if WG -> overload 
	cl_int			ret;
	int             k;

	//--------------- reset liste de contacts
	for (k = 0; k < parameterPtr->nSegment; k++) {
		nForce_per_Segment_hostPtr[k] = 0;
		nMoment_per_Segment_hostPtr[k] = 0;
	}
	ret = clEnqueueWriteBuffer(command_queue, nForce_per_Segment_devPtr, CL_TRUE, 0, parameterPtr->nSegment * sizeof(cl_uint), nForce_per_Segment_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");
	ret = clEnqueueWriteBuffer(command_queue, nMoment_per_Segment_devPtr, CL_TRUE, 0, parameterPtr->nSegment * sizeof(cl_uint), nMoment_per_Segment_hostPtr, 0, NULL, NULL);checkError(ret, "clEnqueueWriteBuffer");

	//--------------- execute kernel_Contact
	ret = clEnqueueNDRangeKernel(command_queue, kernel_CalculateContactForce, 1, NULL, global_Size, local_Size, 0, NULL, NULL);checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_AddContactForce(parameter* parameterPtr) {

	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG };
	cl_int			ret;

	ret = clEnqueueNDRangeKernel(command_queue, kernel_AddContactForce, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_ShiftContact(parameter* parameterPtr) {

	size_t global_Size[1] = { parameterPtr->nContactMax };
	size_t local_Size[1] = { WG };
	cl_int			ret;

	ret = clEnqueueNDRangeKernel(command_queue, kernel_ShiftContact, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_kernel_Integrate_and_Shift(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG };

	//---------------------------;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Integrate_and_Shift, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_kernel_Compute_m1bar(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG };

	//---------------------------;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Compute_m1bar, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_kernel_Compute_m1(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1] = { parameterPtr->nSegment };
	size_t local_Size[1] = { WG };

	//---------------------------;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Compute_m1, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_kernel_Specific_Position(parameter* parameterPtr) {

	cl_int			ret;
	size_t global_Size[1];
	size_t local_Size[1];

	//---------------------------
	global_Size[0] = 1;
	local_Size[0] = 1;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Specific_OneTime, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	global_Size[0] = parameterPtr->nSegment;
	local_Size[0] = WG;
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Specific_Position, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
void Execute_Kernel_Max_Displacement(parameter* parameterPtr) {

	cl_int              ret;
	size_t              global_Size[1] = { parameterPtr->nSegment };
	size_t              local_Size[1] = { SIZE_REDUCTION };

	//--------------- execute kernel_Max_Displacement
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Max_Displacement_Step1, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	global_Size[0] = 1;
	local_Size[0] = 1;

	ret = clEnqueueNDRangeKernel(command_queue, kernel_Max_Displacement_Step2, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");
}

//-------------------------------------------------------
unsigned int Execute_Kernel_NContact_From_Device(parameter* parameterPtr) {

	cl_int              ret;
	size_t              global_Size[1] = { parameterPtr->nContactMax };
	size_t              local_Size[1] = { SIZE_REDUCTION };

	//--------------- execute kernel_Max_Displacement
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Count_Contact_Step1, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	global_Size[0] = 1;
	local_Size[0] = 1;

	ret = clEnqueueNDRangeKernel(command_queue, kernel_Count_Contact_Step2, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	return Read_One_UInt_Param_dev2host(7);
}

//-------------------------------------------------------
float Execute_Kernel_Max_Delta(parameter* parameterPtr) {

	cl_int              ret;
	size_t              global_Size[1] = { parameterPtr->nContactMax };
	size_t              local_Size[1] = { SIZE_REDUCTION };

	//--------------- execute kernel_Max_Displacement
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Max_Delta_Step1, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	global_Size[0] = 1;
	local_Size[0] = 1;

	ret = clEnqueueNDRangeKernel(command_queue, kernel_Max_Delta_Step2, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	return Read_One_Double_Param_dev2host(81);
}

//-------------------------------------------------------
void Execute_Kernel_Unwarp(parameter* parameterPtr) {

	size_t global_Size[1] = { parameterPtr->nFiber };
	size_t local_Size[1] = { 1 };
	cl_int			ret;

	//--------------- execute kernel_Max_Displacement
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Unwarp, 1, NULL, global_Size, local_Size, 0, NULL, NULL);
	checkError(ret, "clEnqueueNDRangeKernel");

	return;
}

//-------------------------------------------------------
void Execute_Kernel_Dissipation(parameter* parameterPtr) {

	Execute_Kernel_Sum_One_Array_In_Slot(fiber_dissip_visc_global_devPtr, parameterPtr->nSegment, 33);
	Execute_Kernel_Sum_One_Array_In_Slot(fiber_dissip_visc_stretch_devPtr, parameterPtr->nSegment, 34);
	Execute_Kernel_Sum_One_Array_In_Slot(cont_dissip_ft_devPtr, parameterPtr->nContactMax, 35);
	Execute_Kernel_Sum_One_Array_In_Slot(cont_dissip_fn_devPtr, parameterPtr->nContactMax, 36);
	Execute_Kernel_Sum_One_Array_In_Slot(fiber_Wop_devPtr, parameterPtr->nSegment, 37);
}


//-------------------------------------------------------
void Execute_Kernel_Sum_One_Array_In_Slot(cl_mem dataBuffer, size_t size, unsigned int slot) {

	cl_int              ret;
	size_t              global_Size[1], local_Size[1];

	global_Size[0] = size;
	local_Size[0] = SIZE_REDUCTION;
	ret = clSetKernelArg(kernel_Dissipation_Step1, 0, sizeof(cl_mem), (void*)&dataBuffer); checkError(ret, "clSetKernelArg");
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Dissipation_Step1, 1, NULL, global_Size, local_Size, 0, NULL, NULL); checkError(ret, "clEnqueueNDRangeKernel");
	
	global_Size[0] = 1;
	local_Size[0] = 1;
	Write_One_UInt_Param_host2dev(NPARAM - 2, size);
	Write_One_UInt_Param_host2dev(NPARAM - 1, slot);
	ret = clEnqueueNDRangeKernel(command_queue, kernel_Dissipation_Step2, 1, NULL, global_Size, local_Size, 0, NULL, NULL); checkError(ret, "clEnqueueNDRangeKernel");
}