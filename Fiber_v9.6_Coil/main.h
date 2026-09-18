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

//---------------- core.c
void			DoOneIteration(void);

//---------------- coil.c
void			SetFibersParameters(void);
void			SetFibersInitialPositions(void);