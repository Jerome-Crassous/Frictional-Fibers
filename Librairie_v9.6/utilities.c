#include <stdlib.h>
#include <stdio.h>
#include "..\fiberLib_OpenCL_v9.6.h"
#include "..\fiberLib_Common_Macros_v9.6.h"

#define N_FIBER_MAX		1000000

//-----------------------------------------------------
unsigned int NSegment(fiber* theFiberPtr, parameter* theParameterPtr, unsigned int iFiber, unsigned int i) {

	unsigned int				iFiberLoc;
	static int					inited = 0;
	static unsigned int			*shift;

	if (!inited) {
		if(theParameterPtr->nFiber> N_FIBER_MAX)
			printf("##### theParameterPtr->nFiber> N_FIBER_MAX in NSegment\n");

		shift = (unsigned int*)malloc(N_FIBER_MAX * sizeof(unsigned int));

		for (iFiberLoc = 0; iFiberLoc < theParameterPtr->nFiber; iFiberLoc++) {
			if (iFiberLoc == 0)
				shift[iFiberLoc] = 0;
			else
				shift[iFiberLoc] = shift[iFiberLoc - 1] + theFiberPtr[iFiberLoc - 1].n;
		}
		inited = 1;
	}
	if ((shift[iFiber] + i) >= theParameterPtr->nSegment)
		printf("##### warning nSegment\n");

	return shift[iFiber] + i;
}