#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

//-----------------------------------------------------
double MyTimeLib(void) {//retourne le temps en seconde

	double                  times;
	static unsigned int				inited = 0;
	static LARGE_INTEGER	lpFrequency, lpPerformanceCountBegin;
	LARGE_INTEGER			lpPerformanceCount;

	if (!inited) {
		QueryPerformanceFrequency(&lpFrequency);
		QueryPerformanceCounter(&lpPerformanceCountBegin);
		inited = 1;
	}

	QueryPerformanceCounter(&lpPerformanceCount);

	times = (double)(lpPerformanceCount.QuadPart - lpPerformanceCountBegin.QuadPart);
	times = times / ((double)lpFrequency.QuadPart);
	times *= 1;
	return(times);
}