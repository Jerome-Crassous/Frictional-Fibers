#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//---------------------------------------------------------
#include <windows.h>
#include <winbase.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <CL/cl.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

#define GN_CELL_X			16
#define GN_CELL_Y			16

#define NLOOP_RAW			284
#define NLOOP				120
#define SCALE				(119./283.)

#define LC_LW_RAW			70.
#define LC_LW				(LC_LW_RAW*SCALE) //29.63

//---------------- core.c
void			DoOneIteration(void);

//---------------- packing.c
void			SetFibersParameters(void);
void			SetFibersInitialPositions(void);