#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <CL/cl.h>

//------------------------ types 
typedef struct parameter {
	unsigned int		iter, line;
	unsigned int		periodic, periodicVisu;
	unsigned int		nFiber, nSegment, nContact, nContactMax;
	unsigned int		iter0,iter1,low_High_Eta,iter3, rRate, low_High_Mu,iter6,iter7,iter8,iterStop;	//20-29
	unsigned int		flag0,flag1,flag2,flag3,flag4,flag5,flag6,flag7,flag_Phase,flag_Reset_Contact;	//30-39
	unsigned int		NN;	//40
	
	double				dt,kn,kt,lambda,lambda_internal,lambda_contact_n, lambda_contact_t;
	double				mu, muTarget, lambdaTarget;	//7-10-11
	double				traction, pressure, pressureTarget;	//15-16-17
	double				ls, bending, c;
	double				d, dTarget, r, rTarget, R, RTarget;	//18-19-20-21-22-23
	double				h, hTarget, dhdt;	//24-25-26
	double				RHelix;	//param_Double[27]
	double				pitch;	//param_Double[28]
	double				epsStar;	//param_Double[30]
	double				fix, fx, fy, fz;	//40-41-42-43
	double				sc, thetac;//47-48
	double				lx, ly, lz;	//50-51-52
	double				lxTarget, lyTarget, lzTarget;	//53-54-55
	double				nCellX, nCellY;	//56-57
	double				bound, boundX, boundY, boundZ; //60-61-62-63
	double				alpha, alphaX, alphaY, alphaZ;	//64-65-66-67
	double				theta, beta;	//68-69
	double				delta;	//70
}parameter;

typedef struct fiber {
	unsigned int		n, status;
	double				radius, masse, l0, k0, bending, j, c;
	double				** xt, ** xtm, ** thetat, ** thetatm, ** e, ** f;
	double				** moment, ** m1_bar, **m1;
	unsigned int		* flag;
	double				* lt, * ltm, * kappa1_bar, * kappa2_bar;
}fiber;

typedef struct contact {
	unsigned int		type1, type2, fiber1, fiber2, node1, node2, exist;
	double				s1, s2, deltat, deltatm, W1, W2;
	double				C_t[3], n_t[3], u_t[3], u_tm[3], fc[3];
	double				shift[3];
}contact;

//-----------------fiber.c------------------------------------------
void AllocateOneFiberLib(fiber* theFiberPtr, unsigned int iFiber, unsigned int npt);
 
//---------------- energy.c
double KineticEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr);
double BendingEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr);
double TwistEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr);
double TractionEnergyLib(fiber* theFiberPtr, parameter* theParameterPtr);
double Elastic_Normal_Contact_Energy(contact* theContactPtr, parameter* theParameterPtr);
double Elastic_Tangential_Contact_Energy(contact* theContactPtr, parameter* theParameterPtr);
void InternalForcesLib(fiber theFiber, unsigned int i, int direction, double* force);
void InternalTractionForcesLib(fiber theFiber, unsigned int i, int direction, double* force);
void InternalBendingForcesLib(fiber theFiber, unsigned int i, int direction, double* force);
void InternalBendingMomentLib(fiber theFiber, unsigned int i, int direction, double* moment);

//--------------- geometry.c
void Calculate_Relaxed_M1_Bar(fiber* theFiberPtr);
void Calculate_M1_From_M1Bar_And_Theta(fiber* theFiberPtr);
void CalculateLeLib(fiber* theFiberPtr);
double KappaLib(fiber theFiber, unsigned int i, int direction);
double	ZMaxLib(fiber* theFiberPtr, parameter* theParameterPtr);
double	ZMinLib(fiber* theFiberPtr, parameter* theParameterPtr);

//--------------- print.c
void PrintPositionsLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int index);
void PrintForcesLib(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int index);
void Print_leLib(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int index);
void Print_thetaLib(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int index);
double PrintBendingEnergyPerFiber(fiber* theFiberPtr, parameter* theParameterPtr);

//--------------- read_save_config.c
void SaveConfigLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int n);
void ReadConfigLib(fiber** theFiberPtr, parameter* theParameterPtr, char* pathName, int n);
void SaveReducedConfigLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int n);
int ReadReducedConfigLib(fiber** theFiberPtr, parameter* theParameterPtr, char* pathName, int n);
void PrintReducedConfigLib(fiber* theFiberPtr, parameter* theParameterPtr, char* pathName, int n);

//--------------- utilities.c
unsigned int NSegment(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int iFiber, unsigned int i);

//---------------- math.c ------------------------
double			MyRand(void);
double			MyRandWithSeed(long seed);
double			DotProduct(double* v1, double* v2);
double			CosTheta(double* v1, double* v2);
double			Distance(double* v1, double* v2);
double			Norm(double* v);
double			Normalize(double *v);
void			CrossProduct(double* v1, double* v2, double* cross);
void			VectorMultiply(double* v, double scale);
int				IMin(int a, int b);
int				IMax(int a, int b);
double			DMin(double a, double b);
double			DMax(double a, double b);
double			DAbs(double a);
double			DSign(double a) ;
double			Angle(double* v1, double* v2);
void			SetRandomOriention(double* e, double* m1);
void Rotate_Point_Around_Axis(double* x, double* pointOnAxis, double* direction, double angle);
void Shift_Point(double* x, double* shift);

//----------------- device_init.c
void			InitDevice(void);
void			Create_Kernels(void);
void			Execute_Kernel_Calculate_le(parameter* parameterPtr);
void			Execute_Kernel_Bending_Force(parameter* parameterPtr);
void			Execute_Kernel_Twist_Force(parameter* parameterPtr);
void			Execute_Kernel_Specific_Force(parameter* parameterPtr);
void            Execute_Kernel_PossibleContact(parameter* parameterPtr);
void            Execute_Kernel_Contact(parameter* parameterPtr);
void			Execute_Kernel_RemoveDoubleContact(parameter* parameterPtr);
void			Execute_Kernel_CalculateContactForce(parameter* parameterPtr);
void			Execute_Kernel_AddContactForce(parameter* parameterPtr);
void			Execute_Kernel_ShiftContact(parameter* parameterPtr);
void			Execute_kernel_Integrate_and_Shift(parameter* parameterPtr);
void			Execute_kernel_Compute_m1bar(parameter* parameterPtr);
void			Execute_kernel_Compute_m1(parameter* parameterPtr);
void			Execute_kernel_Specific_Position(parameter* parameterPtr);
void			Execute_Kernel_Max_Displacement(parameter* parameterPtr);
float			Execute_Kernel_Max_Delta(parameter* parameterPtr);
unsigned int	Execute_Kernel_NContact_From_Device(parameter* parameterPtr);
void			Execute_Kernel_Unwarp(parameter* parameterPtr);
void			Execute_Kernel_Dissipation(parameter* parameterPtr);

//----------------- device_host_device.c
void			Create_Host_Ptrs(void);
void			Create_Device_Ptrs(void);

void            Copy_Utils_host2dev(void);
void			Copy_Param_host2dev(void);
void			Copy_Param_dev2host(void);

void			Write_One_UInt_Param_host2dev(unsigned int slot, unsigned int value);
void			Write_One_Double_Param_host2dev(unsigned int slot, double value);
unsigned int	Read_One_UInt_Param_dev2host(unsigned int slot);
double			Read_One_Double_Param_dev2host(unsigned int slot);

void			Copy_Fiber_host2dev(void);
void			Copy_Fiber_dev2host(void);
void			Copy_FiberPosition_dev2host(void);
void			Copy_Contact_dev2host(void);
void			Print_Possible_dev2host(unsigned int flag);
unsigned int	NPossContact_dev(void);
unsigned int	NContact_dev(void);

void			ResetDissipationBuffers(void);
void			ResetW12Buffers(void);

//----------------- device_error_code.c
void			check_error(cl_int err, const char* operation, char* filename, int line);
#define checkError(E, S) check_error(E, S, __FILE__, __LINE__)

//----------------- times.c
double MyTimeLib(void);